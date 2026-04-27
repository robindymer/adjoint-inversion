using Plots
gr() # GR backend

# --- CONFIGURATION ---
data_dir = "./data/marmousi2"
dx = 1.25
dz = 1.25

# Decimation factor (Integers only)
# 1 = Full resolution (Slow)
# 4 = Moderate speed (Good for final checks)
# 10 = Fast (Good for quick debugging)
downsample_factor = 4 

# --- HELPER FUNCTIONS ---

function get_dimensions(dir_path)
    dim_file = joinpath(dir_path, "dimensions.txt")
    if isfile(dim_file)
        dims = split(read(dim_file, String))
        return parse(Int, dims[1]), parse(Int, dims[2])
    else
        println("Warning: dimensions.txt not found. Using defaults.")
        return 13601, 2801
    end
end

function read_seismic_bin(filename, nx, nz)
    # Read full binary first (Fastest sequential read)
    data = Array{Float64}(undef, nz, nx)
    open(filename, "r") do io
        read!(io, data)
    end
    return data
end

# --- MAIN EXECUTION ---

# 1. Get Original Dimensions
nx_orig, nz_orig = get_dimensions(data_dir)
println("Original Dimensions: $nx_orig (x) × $nz_orig (z)")

# 2. Load Data (Full Resolution)
println("Loading binaries...")
# Note: We load full data first to ensure binary alignment is correct, 
# then we slice it. Reading strided bytes from disk is harder/slower.
lam_full = read_seismic_bin(joinpath(data_dir, "lambda.bin"), nx_orig, nz_orig)
mu_full  = read_seismic_bin(joinpath(data_dir, "mu.bin"), nx_orig, nz_orig)
rho_full = read_seismic_bin(joinpath(data_dir, "rho.bin"), nx_orig, nz_orig)

# 3. Downsample
if downsample_factor > 1
    println("Downsampling by factor of $downsample_factor...")
    
    # Select indices: 1, 1+step, 1+2*step, ...
    # Julia ranges are [start:step:stop]
    idx_z = 1:downsample_factor:nz_orig
    idx_x = 1:downsample_factor:nx_orig
    
    # Slice the arrays
    lam = lam_full[idx_z, idx_x]
    mu  = mu_full[idx_z, idx_x]
    rho = rho_full[idx_z, idx_x]
    
    # Clear full arrays from memory (optional, helps garbage collector)
    lam_full = nothing
    mu_full = nothing
    rho_full = nothing
    GC.gc() # Force garbage collection
    
    # Update dimensions and spacing for the plot
    nx, nz = size(lam, 2), size(lam, 1)
    dx_plot = dx * downsample_factor
    dz_plot = dz * downsample_factor
else
    lam, mu, rho = lam_full, mu_full, rho_full
    nx, nz = nx_orig, nz_orig
    dx_plot, dz_plot = dx, dz
end

println("Plotting Grid: $nx (x) × $nz (z)")

# 4. Plotting
println("Generating plots...")

# Define physical axes ranges (in meters) using the NEW spacing
xs = range(0, stop=(nx-1)*dx_plot, length=nx)
zs = range(0, stop=(nz-1)*dz_plot, length=nz)

plot_settings = (
    yflip = true, 
    xlabel = "Distance (m)", 
    ylabel = "Depth (m)", 
    aspect_ratio = :equal
)

# Plot 1: Lambda
p1 = heatmap(xs, zs, lam, 
    title = "Lambda (λ) - Decimated $(downsample_factor)x", 
    c = :viridis, 
    clabel = "Pa"; 
    plot_settings...
)

# Plot 2: Mu
p2 = heatmap(xs, zs, mu, 
    title = "Shear Modulus (μ)", 
    c = :plasma, 
    clabel = "Pa"; 
    plot_settings...
)

# Plot 3: Density
p3 = heatmap(xs, zs, rho, 
    title = "Density (ρ)", 
    c = :cividis, 
    clabel = "kg/m³"; 
    plot_settings...
)

# Combine
final_plot = plot(p1, p2, p3, layout = (3, 1), size = (1000, 1200))

# Display
display(final_plot)

# Save
output_file = joinpath(data_dir, "marmousi_plots_$(downsample_factor)x.png")
savefig(output_file)
println("Saved plot to $output_file")