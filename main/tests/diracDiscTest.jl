using LinearAlgebra
using Plots

include("../types.jl")
# include("grid.jl")
include("../stencils.jl")
include("../diracDisc.jl")

function test_dirac_1D()
    println("=" ^ 60)
    println("Testing 1D Dirac Delta Discretization")
    println("=" ^ 60)
    
    # Create 1D grid
    nx = 101
    x_min, x_max = -1.0, 1.0
    x = range(x_min, x_max, length=nx)
    h = x[2] - x[1]
    
    # Create norm matrix (identity for simplicity, or use SBP norm)
    H = Diagonal(h * ones(nx))
    quadrature = diag(H)
    
    # Source location
    x_s = 0.3
    
    # Orders
    m_order = 3  # moment conditions
    s_order = 3  # smoothness conditions
    
    # Compute delta function
    delta = diracDiscr1D(x_s, collect(x), m_order, s_order, quadrature)
    
    # Check moment conditions
    println("\nMoment condition checks:")
    for k in 0:m_order-1
        moment = sum(delta .* (collect(x) .- x_s).^k .* h)
        expected = k == 0 ? 1.0 : 0.0
        println("  ∫ δ(x) (x - x_s)^$k dx = $moment (expected: $expected)")
    end
    
    # Plot
    p1 = plot(x, delta, 
              xlabel="x", ylabel="δ(x)", 
              title="1D Dirac Delta Discretization\nSource at x = $x_s",
              label="Delta function",
              linewidth=2,
              marker=:circle,
              markersize=3)
    vline!([x_s], label="Source location", linestyle=:dash, linewidth=2)
    
    savefig(p1, "dirac_1D_test.png")
    println("\nPlot saved to dirac_1D_test.png")
    
    return delta
end

function test_dirac_2D()
    println("\n" * "=" ^ 60)
    println("Testing 2D Dirac Delta Discretization")
    println("=" ^ 60)
    
    # Create 2D grid
    nx, ny = 51, 51
    x_min, x_max = -1.0, 1.0
    y_min, y_max = -1.0, 1.0
    
    x = range(x_min, x_max, length=nx)
    y = range(y_min, y_max, length=ny)
    hx = x[2] - x[1]
    hy = y[2] - y[1]
    
    # Create Grid struct (adjust based on your Grid definition)
    g = Grid(x_min, x_max, y_min, y_max, nx, ny, hx, hy, collect(x), collect(y))
    # Create norm matrices
    quadrature_x = diag(Diagonal(hx * ones(nx)))
    quadrature_y = diag(Diagonal(hy * ones(ny)))
    quadrature = [quadrature_x, quadrature_y]
    
    # Source location
    x_s = [0.25, -0.3]
    
    # Orders
    m_order = 3
    s_order = 3
    
    # Compute 2D delta function
    delta_2D = diracDiscr(g, x_s, m_order, s_order, quadrature)
    
    # Reshape to 2D for visualization
    delta_matrix = reshape(delta_2D, nx, ny)
    
    # Check integral (should be 1)
    integral = sum(delta_2D) * hx * hy
    println("\n∫∫ δ(x,y) dx dy = $integral (expected: 1.0)")
    
    # Create meshgrid for plotting
    X = repeat(collect(x), 1, ny)
    Y = repeat(collect(y)', nx, 1)
    
    # Plot 1: Surface plot
    p1 = surface(x, y, delta_matrix',
                 xlabel="x", ylabel="y", zlabel="δ(x,y)",
                 title="2D Dirac Delta - Surface Plot",
                 camera=(45, 30),
                 colorbar=true)
    
    # Plot 2: Contour plot
    p2 = contourf(x, y, delta_matrix',
                  xlabel="x", ylabel="y",
                  title="2D Dirac Delta - Contour Plot",
                  colorbar=true,
                  levels=20)
    scatter!([x_s[1]], [x_s[2]], 
             marker=:star, markersize=10, 
             color=:red, label="Source")
    
    # Plot 3: Heatmap
    p3 = heatmap(x, y, delta_matrix',
                 xlabel="x", ylabel="y",
                 title="2D Dirac Delta - Heatmap",
                 colorbar=true,
                 aspect_ratio=1)
    scatter!([x_s[1]], [x_s[2]], 
             marker=:star, markersize=10, 
             color=:red, label="Source")
    
    # Plot 4: Cross-sections
    ix_s = argmin(abs.(collect(x) .- x_s[1]))
    iy_s = argmin(abs.(collect(y) .- x_s[2]))
    
    p4 = plot(x, delta_matrix[:, iy_s],
              xlabel="x", ylabel="δ(x, y_s)",
              label="Cross-section at y = $(x_s[2])",
              linewidth=2,
              marker=:circle)
    plot!(y, delta_matrix[ix_s, :],
          xlabel="y", ylabel="δ(x_s, y)",
          label="Cross-section at x = $(x_s[1])",
          linewidth=2,
          marker=:square)
    title!("Cross-sections through source point")
    
    # Combine plots
    p_combined = plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 1000))
    
    savefig(p_combined, "dirac_2D_test.png")
    println("Plot saved to dirac_2D_test.png")
    
    return delta_2D, delta_matrix
end

function test_convergence_1D()
    println("\n" * "=" ^ 60)
    println("Testing Convergence of 1D Dirac Delta")
    println("=" ^ 60)
    
    x_s = 0.3
    m_order = 4
    s_order = 4
    
    # Test different grid resolutions
    n_values = [21, 41, 81, 161, 321]
    errors = Float64[]
    
    println("\nGrid points | Integral error")
    println("-" ^ 35)
    
    for nx in n_values
        x = range(-1.0, 1.0, length=nx)
        h = x[2] - x[1]
        H = Diagonal(h * ones(nx))
        
        delta = diracDiscr1D(x_s, collect(x), m_order, s_order, Matrix(H))
        
        # Check if integral equals 1
        integral = sum(delta) * h
        error = abs(integral - 1.0)
        push!(errors, error)
        
        println("$nx          | $error")
    end
    
    # Plot convergence
    p = plot(n_values, errors,
             xlabel="Number of grid points",
             ylabel="Absolute error in ∫δ(x)dx",
             title="Convergence of Dirac Delta Discretization",
             label="Error",
             marker=:circle,
             markersize=6,
             linewidth=2,
             yscale=:log10,
             xscale=:log10)
    
    savefig(p, "dirac_convergence_test.png")
    println("\nConvergence plot saved to dirac_convergence_test.png")
end

function test_different_orders()
    println("\n" * "=" ^ 60)
    println("Testing Different Order Combinations")
    println("=" ^ 60)
    
    nx = 101
    x = range(-1.0, 1.0, length=nx)
    h = x[2] - x[1]
    H = Diagonal(h * ones(nx))
    x_s = 0.2
    
    orders = [(2, 2), (3, 3), (4, 4), (5, 5)]
    
    p = plot(xlabel="x", ylabel="δ(x)",
             title="Dirac Delta with Different Orders",
             legend=:topright)
    
    for (m_ord, s_ord) in orders
        delta = diracDiscr1D(x_s, collect(x), m_ord, s_ord, Matrix(H))
        plot!(x, delta, label="m=$m_ord, s=$s_ord", linewidth=2)
    end
    
    vline!([x_s], label="Source", linestyle=:dash, color=:black, linewidth=2)
    
    savefig(p, "dirac_orders_test.png")
    println("\nOrder comparison plot saved to dirac_orders_test.png")
end

# Run all tests
function run_all_tests()
    println("\n" * "█" ^ 60)
    println("Running Dirac Delta Discretization Tests")
    println("█" ^ 60)
    
    test_dirac_1D()
    test_dirac_2D()
    # test_convergence_1D()
    # test_different_orders()
    
    println("\n" * "█" ^ 60)
    println("All tests completed!")
    println("█" ^ 60)
end

# Run tests if executed directly
# if abspath(PROGRAM_FILE) == @__FILE__
run_all_tests()
# end