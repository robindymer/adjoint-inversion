function setup_mms!(mx, my, x_values, y_values, u1_0, u2_0, v1_0, v2_0, mms, ω)
    for i in 1:mx
        for j in 1:my
            u1_0[idx(i, j, mx, my)] = mms.u1(x_values[i], y_values[j], 0.0, ω)
            u2_0[idx(i, j, mx, my)] = mms.u2(x_values[i], y_values[j], 0.0, ω)
            v1_0[idx(i, j, mx, my)] = mms.u1_t(x_values[i], y_values[j], 0.0, ω)
            v2_0[idx(i, j, mx, my)] = mms.u2_t(x_values[i], y_values[j], 0.0, ω)
        end
    end
end

function simulate(x_l, x_r, y_l, y_r, mx, my, D, λ, μ, ω, T, CFL, SBP_ORDER)
    # First, set up the grid and stencils
    x_width = x_r - x_l
    y_width = y_r - y_l
    if Params.PERIODIC
        hx = x_width / mx
        hy = y_width / my
        x_values = LinRange(x_l, x_r - hx, mx)
        y_values = LinRange(y_l, y_r - hy, my)
    else
        hx = x_width / (mx - 1)
        hy = y_width / (my - 1)
        x_values = LinRange(x_l, x_r, mx)
        y_values = LinRange(y_l, y_r, my)
    end


    grid = Grid(x_l, x_r, y_l, y_r, mx, my, hx, hy, x_values, y_values)
    stencils = compute_stencils(mx, my, hx, hy, SBP_ORDER)
    material = Material(D, λ, μ, ω)
    bcs = BoundaryConditions(Params.BC_WEST, Params.BC_EAST, Params.BC_SOUTH, Params.BC_NORTH)

    anim = Animation()
    # Simulation parameters
    c = sqrt(maximum(λ) + 2 * maximum(μ)) / sqrt(minimum(D))  # max wave speed
    dt = CFL * min(hx, hy) / c
    # Convert dt to float32 
    dt = Params.PRECISION(dt)
    println("dt: $dt")
    println("Maximum speed: $c")
    t = 0.0
    # Plotting parameters
    n_step = Params.output_interval
    steps = 0

    v1_0 = zeros(mx * my)
    v2_0 = zeros(mx * my)
    u1_0 = zeros(mx * my)
    u2_0 = zeros(mx * my)

    if Params.USE_MMS
        println("Using MMS test")
        setup_mms!(mx, my, x_values, y_values, u1_0, u2_0, v1_0, v2_0, mms, ω)
    else
        println("Using Gaussian initial condition")
        # Using Gaussian initial condition
        for i in 1:mx
            for j in 1:my
                u1_0[idx(i, j, mx, my)] = initialcondition(x_values[i], y_values[j])[1]
                u2_0[idx(i, j, mx, my)] = initialcondition(x_values[i], y_values[j])[2]
            end
        end
    end

    u1 = nothing
    u2 = nothing

    v1 = copy(v1_0)
    v2 = copy(v2_0)
    u1 = copy(u1_0)
    u2 = copy(u2_0)

    # Initialize VTK variables (will remain nothing if SAVE_VTK is false)
    folder_name = nothing
    pvd = nothing

    if Params.SAVE_VTK
        # Save VTK files pvd file
        folder_name = Params.OUTPUT_PATH
        mkpath(folder_name)
        pvd = paraview_collection(Params.OUTPUT_NAME)
        # Write simulation details to text file in folder
        open("$folder_name/simulation.txt", "w") do f
            write(f, Params.format_parameters(dt=dt, hx=hx, hy=hy))
        end
    end

    derivatives = Derivatives(zeros(mx * my), zeros(mx * my), zeros(mx * my), zeros(mx * my),
        zeros(mx * my), zeros(mx * my), zeros(mx * my), zeros(mx * my),
        zeros(mx * my), zeros(mx * my), zeros(mx * my), zeros(mx * my),
        zeros(mx * my), zeros(mx * my))
    iD = 1.0 ./ D

    device_data = nothing
    rk_cache = nothing
    D2_stencils = nothing

    if Params.SCHEMA in (:narrow, :mixed, :mixed_upwind_μ, :mixed_upwind_λ)
        # Precompute D2 stencils
        D2x_λ_∂Ω_W_vec, D2x_λ_∂Ω_E_vec, D2x_μ_∂Ω_W_vec, 
        D2x_μ_∂Ω_E_vec, D2x_λ_vec, D2x_μ_vec, D2y_λ_∂Ω_N_vec, D2y_λ_∂Ω_S_vec, 
        D2y_μ_∂Ω_N_vec, D2y_μ_∂Ω_S_vec, D2y_μ_vec, D2y_λ_vec = precompute_D2_∂Ω(mx, my, grid, stencils, material)

        D2_stencils = D2Stencils(D2x_λ_∂Ω_W_vec, D2x_λ_∂Ω_E_vec, D2x_μ_∂Ω_W_vec, 
        D2x_μ_∂Ω_E_vec, D2x_λ_vec, D2x_μ_vec, D2y_λ_∂Ω_N_vec, D2y_λ_∂Ω_S_vec, 
        D2y_μ_∂Ω_N_vec, D2y_μ_∂Ω_S_vec, D2y_μ_vec, D2y_λ_vec)
    end
    
    if Params.USE_CUDA && CUDA.functional()
        
        USES_UPWIND = Params.SCHEMA in [:upwind, :mixed_upwind_μ, :mixed_upwind_λ]
        if Params.PULSE
            quadrature = [stencils.Hx, stencils.Hy]
            delta_2D = diracDiscr(grid, Params.X_S, Params.M_ORDER, Params.S_ORDER, quadrature)
        end
        to_matrix(v) = CuMatrix(reduce(hcat, v))'
        device_data = DeviceData(
            CuArray(derivatives.u1_x),
            CuArray(derivatives.u1_y),
            CuArray(derivatives.u2_x),
            CuArray(derivatives.u2_y),
            CuArray(derivatives.u2_yλx),
            CuArray(derivatives.u2_xμy),
            CuArray(derivatives.u1_xλy),
            CuArray(derivatives.u1_yμx),
            CuArray(derivatives.u1_xλx),
            CuArray(derivatives.u1_xμx),
            CuArray(derivatives.u1_yμy),
            CuArray(derivatives.u2_yλy),
            CuArray(derivatives.u2_yμy),
            CuArray(derivatives.u2_xμx),
            CuArray(v1),
            CuArray(v2),
            CuArray(u1),
            CuArray(u2),
            CuArray(material.λ),
            CuArray(material.μ),
            CuArray(material.D),
            CuArray(iD),
            USES_UPWIND ? nothing : CuArray(stencils.D1x),
            USES_UPWIND ? nothing : CuArray(stencils.D1y),
            USES_UPWIND ? CuArray(stencils.Dpx) : nothing,
            USES_UPWIND ? CuArray(stencils.Dpy) : nothing,
            USES_UPWIND ? CuArray(stencils.Dmx) : nothing,
            USES_UPWIND ? CuArray(stencils.Dmy) : nothing,
            Params.PULSE ? CuArray(delta_2D) : nothing,
            Params.PERIODIC || USES_UPWIND ? nothing : to_matrix(stencils.D1x_∂Ω_W),
            Params.PERIODIC || USES_UPWIND ? nothing : to_matrix(stencils.D1x_∂Ω_E),
            Params.PERIODIC || USES_UPWIND ? nothing : to_matrix(stencils.D1y_∂Ω_N),
            Params.PERIODIC || USES_UPWIND ? nothing : to_matrix(stencils.D1y_∂Ω_S),
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dpx_∂Ω_W) : nothing,
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dpx_∂Ω_E) : nothing,
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dpy_∂Ω_N) : nothing,
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dpy_∂Ω_S) : nothing,
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dmx_∂Ω_W) : nothing,
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dmx_∂Ω_E) : nothing,
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dmy_∂Ω_N) : nothing,
            !Params.PERIODIC && USES_UPWIND ? to_matrix(stencils.Dmy_∂Ω_S) : nothing,
            Params.PERIODIC ? nothing : CuArray(stencils.d1x_l),
            Params.PERIODIC ? nothing : CuArray(stencils.d1x_r),
            Params.PERIODIC ? nothing : CuArray(stencils.d1y_l),
            Params.PERIODIC ? nothing : CuArray(stencils.d1y_r),
            Params.PERIODIC || USES_UPWIND ? nothing : CuArray(stencils.d1x_l_wide),
            Params.PERIODIC || USES_UPWIND ? nothing : CuArray(stencils.d1x_r_wide),
            Params.PERIODIC || USES_UPWIND ? nothing : CuArray(stencils.d1y_l_wide),
            Params.PERIODIC || USES_UPWIND ? nothing : CuArray(stencils.d1y_r_wide),
            D2_stencils === nothing ? nothing : CuArray(D2x_λ_∂Ω_W_vec), D2_stencils === nothing ? nothing : CuArray(D2x_λ_∂Ω_E_vec), 
            D2_stencils === nothing ? nothing : CuArray(D2x_μ_∂Ω_W_vec), D2_stencils === nothing ? nothing : CuArray(D2x_μ_∂Ω_E_vec), 
            D2_stencils === nothing ? nothing : CuMatrix(D2x_λ_vec), D2_stencils === nothing ? nothing : CuMatrix(D2x_μ_vec), 
            D2_stencils === nothing ? nothing : CuArray(D2y_λ_∂Ω_N_vec), D2_stencils === nothing ? nothing : CuArray(D2y_λ_∂Ω_S_vec), 
            D2_stencils === nothing ? nothing : CuArray(D2y_μ_∂Ω_N_vec), D2_stencils === nothing ? nothing : CuArray(D2y_μ_∂Ω_S_vec), 
            D2_stencils === nothing ? nothing : CuMatrix(D2y_μ_vec), D2_stencils === nothing ? nothing : CuMatrix(D2y_λ_vec)
        )

        # Ensure inputs are on GPU
        v1 = device_data.v1
        v2 = device_data.v2
        u1 = device_data.u1
        u2 = device_data.u2

        z = CUDA.zeros(Params.PRECISION, mx*my)
        rhs_cache = RHS_Cache(
            CUDA.zeros(Params.PRECISION, mx*my), CUDA.zeros(Params.PRECISION, mx*my), 
            CUDA.zeros(Params.PRECISION, mx*my), CUDA.zeros(Params.PRECISION, mx*my)
        )
        rk_cache = RK4_Cache(
            copy(z), copy(z), copy(z), copy(z), # k1
            copy(z), copy(z), copy(z), copy(z), # k2
            copy(z), copy(z), copy(z), copy(z), # k3
            copy(z), copy(z), copy(z), copy(z), # k4
            rhs_cache
        )
    end

    if Params.TEST_SYMMETRY
        # Test SBP properties
        # checkSymmetry(grid, stencils, material, derivatives, device_data)
        testSymmetry(grid, stencils, material, derivatives, bcs, device_data, D2_stencils, rk_cache)
        error("Stopping after symmetry check")
    end

    while t < T
        if t + dt > T
            dt = T - t
        end
        if Params.DEBUG_CUDA
            RK4(v1, v2, u1, u2, t, dt, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rk_cache)
            CUDA.synchronize()
            data = CUDA.@profile RK4(v1, v2, u1, u2, t, dt, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rk_cache)
            display(data)
        else
            v1, v2, u1, u2 = RK4(v1, v2, u1, u2, t, dt, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rk_cache)
        end

        if (steps % n_step == 0)
            plot_simulation(v1, v2, grid, t, steps, anim, ω, u1, u2, folder_name, pvd)
            println("t=$t")
        end

        if Params.DEBUG_CUDA
            break
        end
        t += dt
        steps += 1
    end
    if Params.SAVE_GIF
        plot_simulation(v1, v2, grid, t, steps, anim, ω, u1, u2, folder_name, pvd)
        gif(anim, "$(Params.ASSETS_PATH)/wave2D_$(Params.SCHEMA)_$(Params.USE_CUDA ? "GPU" : "CPU").gif", fps=10)
    end
    if Params.SAVE_VTK
        plot_simulation(v1, v2, grid, t, steps, anim, ω, u1, u2, folder_name, pvd)
        vtk_save(pvd)
    end
    # Compare with reference solution. First read wave2D_ref.out
    # v1_ref, v2_ref, u1_ref, u2_ref = read_solution("wave2D_ref_$(x_l)-$(x_r)x$(y_l)-$(y_r)_mx$(mx)-my$(my)_T$(T).out", grid)
    # L2_err_v, L2_err_u = compare_solutions(v1, v2, u1, u2, v1_ref, v2_ref, u1_ref, u2_ref, grid)
    # println("L2 error in velocity: $L2_err_v")
    # println("L2 error in displacement: $L2_err_u")

    # Write new reference solution
    # write_solution(v1, v2, u1, u2, "wave2D_ref_$(x_l)-$(x_r)x$(y_l)-$(y_r)_mx$(mx)-my$(my)_T$(T).out", grid)

    if Params.USE_MMS
        L2_err = compute_L2_error(u1, u2, t, grid, stencils, material, false)
        return L2_err
    else
        return nothing
    end
end
