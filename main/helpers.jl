include("types.jl")

# TODO: Update for variable materials
# compute strains and total energy
function compute_energy(v1, v2, u1, u2, grid::Grid, stencils::Stencils, material::Material)
    D, λ, μ = material.D, material.λ, material.μ
    # Unpack grid
    hx, hy = grid.hx, grid.hy

    A = hx*hy
    # velocities
    Ekin = 0.5 * sum(v1.^2 .+ v2.^2) * A * D

    # strains
    u1_x = apply_D1x(u1, grid, stencils)
    u1_y = apply_D1y(u1, grid, stencils)
    u2_x = apply_D1x(u2, grid, stencils)
    u2_y = apply_D1y(u2, grid, stencils)

    exx = u1_x
    eyy = u2_y
    exy = 0.5 .* (u2_x .+ u1_y)

    # potential energy density per point:
    # W = 0.5*λ*(exx + eyy)^2 + μ*(exx^2 + eyy^2 + 2*exy^2)
    tr = exx .+ eyy
    W = 0.5*λ .* (tr.^2) .+ μ .* (exx.^2 .+ eyy.^2 .+ 2.0 .* exy.^2)
    Epot = sum(W) * A

    return Ekin + Epot, Ekin, Epot
end
function write_solution(v1, v2, u1, u2, filename, grid::Grid)
    open(filename, "w") do io
        for j in 1:grid.my
            for i in 1:grid.mx
                idx = (i-1)*grid.my + j
                println(io, "$(grid.x_values[i]) $(grid.y_values[j]) $(u1[idx]) $(u2[idx]) $(v1[idx]) $(v2[idx])")
            end
            println(io, "")
        end
    end
end

function read_solution(filename, grid::Grid)
    v1 = zeros(grid.mx*grid.my)
    v2 = zeros(grid.mx*grid.my)
    u1 = zeros(grid.mx*grid.my)
    u2 = zeros(grid.mx*grid.my)
    open(filename, "r") do io
        i = 1
        j = 1
        for line in eachline(io)
            if isempty(line)
                i = 1
                j += 1
                continue
            end
            vals = split(line)
            idx = (i-1)*grid.my + j
            u1[idx] = parse(Float64, vals[3])
            u2[idx] = parse(Float64, vals[4])
            v1[idx] = parse(Float64, vals[5])
            v2[idx] = parse(Float64, vals[6])
            i += 1
        end
    end
    return v1, v2, u1, u2
end

function save_error_map(u_num, u_exact_vec, grid::Grid, filename::String)
    Nx, Ny = grid.mx, grid.my
    y_vals = grid.y_values
    err = zeros(Nx, Ny)

    # build explicit error array
    for i in 1:Nx
        for j in 1:Ny
            k = (j-1)*Nx + i
            err[i,j] = abs(u_num[k] - u_exact_vec[k])
        end
    end

    # transpose so x runs horizontally, y vertically
    heatmap(grid.x_values, y_vals, err',
            xlabel="x", ylabel="y", title="Error map",
            colorbar_title="|error|")

    savefig(filename)
    println("Saved error map to $filename")
end

# TODO: Incorporate u2
function compute_L2_error(u1_num, u2_num, t, grid::Grid, stencils::Stencils, material::Material)
    result = 0.0
    # u1_exact = zeros(length(u1_num))
    # TODO: Probably not use closures for H in periodic case
    if Params.USE_CUDA && CUDA.functional()
        u1_num_cpu = Array(u1_num)
        u2_num_cpu = Array(u2_num)
    else
        u1_num_cpu = u1_num
        u2_num_cpu = u2_num
    end
    for j in 1:grid.my
        for i in 1:grid.mx
            e = u1_num_cpu[idx(i, j, grid.mx, grid.my)] - mms.u1(grid.x_values[i], grid.y_values[j], t, material.ω)
            result += e * stencils.Hx[i] * stencils.Hy[j] * e
            # u1_exact[idx(i, j, grid.mx, grid.my)] = mms.u1(grid.x_values[i], grid.y_values[j], t, material.ω)
        end
    end

    L2_err_u1 = sqrt(result)

    # Save error map
    # save_error_map(u1_num, u1_exact, grid, "error_map_u1_m$(grid.mx)_$(Params.SCHEMA).png", staggered)

    return L2_err_u1
end

function compare_solutions(v1_num, v2_num, u1_num, u2_num, v1_ref, v2_ref, u1_ref, u2_ref, grid::Grid)
    err_v1 = v1_num .- v1_ref
    err_v2 = v2_num .- v2_ref
    err_u1 = u1_num .- u1_ref
    err_u2 = u2_num .- u2_ref
    # NOTE: Naive error, but fine if we only care about existance of difference in solutions
    A = grid.hx * grid.hy
    L2_err_v = sqrt( sum(err_v1.^2 .+ err_v2.^2) * A )
    L2_err_u = sqrt( sum(err_u1.^2 .+ err_u2.^2) * A )
    return L2_err_v, L2_err_u
end

function compute_divergence(u1, u2, grid::Grid, stencils::Stencils, material::Material)
    # Unpack grid
    mx, my = grid.mx, grid.my

    div = zeros(mx*my)

    # Compute divergence at inner points
    u1_x, u1_y, u2_x, u2_y, _, _, _, _, _, _ = compute_derivatives(u1, u2, grid, stencils, material)
    for j in 1:my
        for i in 1:mx
            div[idx(i, j, mx, my)] = u1_x[idx(i, j, mx, my)] + u2_y[idx(i, j, mx, my)]
        end
    end

    return div
end

function compute_curl(u1, u2, grid::Grid, stencils::Stencils, material::Material)
    # Unpack grid
    mx, my = grid.mx, grid.my
    # Unpack stencils
    d1x_l, d1x_r = stencils.d1x_l, stencils.d1x_r
    d1y_l, d1y_r = stencils.d1y_l, stencils.d1y_r
    # Unpack material
    λ, μ = material.λ, material.μ

    curl = zeros(mx*my)

    # Compute curl at inner points
    u1_x, u1_y, u2_x, u2_y, _, _, _, _, _, _ = compute_derivatives(u1, u2, grid, stencils, material)
    for j in 1:my
        for i in 1:mx
            curl[idx(i, j, mx, my)] = u2_x[idx(i, j, mx, my)] - u1_y[idx(i, j, mx, my)]
        end
    end

    return curl
end

function plot_gif(anim, x_values, y_values, u1, u2, t, ω, mx, my)
    # Numerical solution
    p1 = surface(x_values, y_values, u1; 
            zlim=(-1, 1), 
            xlabel="x", ylabel="y", zlabel="u1", 
            title="3D view t=$(round(t, digits=2))",
            camera=(30, 30))   # oblique view

    p2 = surface(x_values, y_values, u1;
                zlim=(-1, 1),
                xlabel="x", ylabel="y", zlabel="u1",
                title="Top-down view",
                camera=(0, 90))   # from above
    
    if Params.USE_MMS
        # Analytical solution
        u1_exact = zeros(length(x_values)*length(y_values))
        # u2_exact = zeros(length(x_values)*length(y_values))
        for i in 1:mx
            for j in 1:my
                u1_exact[idx(i, j, mx, my)] = mms.u1(x_values[i], y_values[j], t, ω)
                # u2_exact[idx(i, j, mx, my)] = mms.u2(x_values[i], y_values[j], t, ω)
            end
        end
        p3 = surface(x_values, y_values, u1_exact;
                    zlim=(-1, 1),
                    xlabel="x", ylabel="y", zlabel="u1",
                    title="Exact solution t=$(round(t, digits=2))",
                    camera=(30, 30))   # oblique view
        p4 = surface(x_values, y_values, u1_exact;
                    zlim=(-1, 1),
                    xlabel="x", ylabel="y", zlabel="u1",
                    title="Exact solution Top-down view",
                    camera=(0, 90))   # from above
        
        plt = plot(p1, p2, p3, p4, layout=(2,2), size=(2000,1600))
    else
        plt = plot(p1, p2, layout=(1,2), size=(2000,800))
    end

    frame(anim, plt)
end

function plot_simulation(v1, v2, grid, t, steps, anim, ω, u1=nothing, u2=nothing, folder_name=nothing, pvd=nothing)
    if Params.SAVE_GIF || Params.SAVE_VTK
        # If we run CUDA, transfer data back to CPU for plotting
        if Params.USE_CUDA && CUDA.functional()
            v1_cpu = Array(v1)
            v2_cpu = Array(v2)
            u1_cpu = Array(u1)
            u2_cpu = Array(u2)
        else
            v1_cpu = v1
            v2_cpu = v2
            u1_cpu = u1
            u2_cpu = u2
        end
    end
    if Params.SAVE_GIF
        plot_gif(anim, grid.x_values, grid.y_values, u1_cpu, u2_cpu, t, ω, grid.mx, grid.my)
    end
    if Params.SAVE_VTK
        # div_u = compute_divergence(u1, u2, grid, stencils, material)
        # curl_u = compute_curl(u1, u2, grid, stencils, material)
        # println("Max curl: $(maximum(abs.(curl_u))) at t=$t")
        vtkfile = vtk_grid("$folder_name/field_$(lpad(steps,3,'0'))", grid.x_values, grid.y_values)

        vtk_point_data(vtkfile, u1_cpu, "u1")
        vtk_point_data(vtkfile, u2_cpu, "u2")
        vtk_point_data(vtkfile, v1_cpu, "v1")
        vtk_point_data(vtkfile, v2_cpu, "v2")

        # vtk_point_data(vtkfile, div_u, "divergence")
        # vtk_point_data(vtkfile, curl_u, "curl")
        pvd[t] = vtkfile
        vtk_save(vtkfile)
    end
end

function get_populated_material(λ_tag::String, μ_tag::String, D_tag::String, mx::Int, my::Int)
    x_l = Params.x_l
    x_r = Params.x_r
    y_l = Params.y_l
    y_r = Params.y_r
    λ = zeros(mx*my)
    μ = zeros(mx*my)
    D = zeros(mx*my)

    if Params.PERIODIC
        hx = (x_r - x_l) / mx
        hy = (y_r - y_l) / my
    else
        hx = (x_r - x_l) / (mx - 1)
        hy = (y_r - y_l) / (my - 1)
    end

    for i in 1:mx
        for j in 1:my
            x = x_l + (i - 1) * hx
            y = y_l + (j - 1) * hy

            λ[idx(i,j,mx,my)] = MaterialFunctions.λ_func(λ_tag, x, y)
            μ[idx(i,j,mx,my)] = MaterialFunctions.μ_func(μ_tag, x, y)
            D[idx(i,j,mx,my)] = MaterialFunctions.D_func(D_tag, x, y)
        end
    end

    return λ, μ, D
end

function get_populated_marmousi_material(data_dir::String, mx::Int, my::Int, downsample_factor::Int)
    # 1. Read dimensions
    dim_file = joinpath(data_dir, "dimensions.txt")
    dims = parse.(Int, split(read(dim_file, String)))
    nx_orig, nz_orig = dims[1], dims[2]
    
    println("Loading Marmousi data from $data_dir")
    println("  Original size: $nx_orig x $nz_orig")
    
    # 2. Read binary data
    function read_marmousi_bin(filename, nx, nz)
        data = Array{Float64}(undef, nz, nx)
        open(filename, "r") do io
            read!(io, data)
        end
        return data
    end
    
    lam_full = read_marmousi_bin(joinpath(data_dir, "lambda.bin"), nx_orig, nz_orig)
    mu_full = read_marmousi_bin(joinpath(data_dir, "mu.bin"), nx_orig, nz_orig)
    D_full = read_marmousi_bin(joinpath(data_dir, "rho.bin"), nx_orig, nz_orig)
    
    # 3. Downsample to target grid
    # Calculate downsampling indices to match target mx, my
    idx_x = round.(Int, range(1, nx_orig, length=mx))
    idx_z = round.(Int, range(1, nz_orig, length=my))
    
    lam_2d = lam_full[idx_z, idx_x]
    mu_2d = mu_full[idx_z, idx_x]
    D_2d = D_full[idx_z, idx_x]
    
    println("  Downsampled to: $mx x $my")
    
    # 4. Convert to 1D arrays with simulation indexing
    λ = zeros(mx*my)
    μ = zeros(mx*my)
    D = zeros(mx*my)
    
    for i in 1:mx
        for j in 1:my
            λ[idx(i, j, mx, my)] = lam_2d[j, i]  # Note: 2D array is [z, x] but we need [x, y]
            μ[idx(i, j, mx, my)] = mu_2d[j, i]
            D[idx(i, j, mx, my)] = D_2d[j, i]
        end
    end
    
    println("  Lambda range: $(minimum(λ)) to $(maximum(λ)) Pa")
    println("  Mu range: $(minimum(μ)) to $(maximum(μ)) Pa")
    println("  Density range: $(minimum(D)) to $(maximum(D)) kg/m³")
    
    return λ, μ, D
end