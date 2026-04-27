# Get (v_t, u_t) from (v, u)
function RHS(v1, v2, u1, u2, t, grid::Grid, stencils::Stencils, material::Material,
    bcs::BoundaryConditions, derivatives::Derivatives, device_data::Union{DeviceData,Nothing},
    D2_stencils::Union{D2Stencils, Nothing}, cache::Union{RHS_Cache, Nothing})

    # Unpack grid
    mx, my = grid.mx, grid.my
    x_values, y_values = grid.x_values, grid.y_values
    quad_term_x, quad_term_y = Params.PRECISION(stencils.quad_term_x), Params.PRECISION(stencils.quad_term_y)
    # Unpack material
    D, λ, μ, ω = material.D, material.λ, material.μ, material.ω
    # Inverse density
    iD = 1.0 ./ D

    if Params.USE_CUDA && CUDA.functional()
        device_data.u1 .= u1
        device_data.u2 .= u2
        device_data.v1 .= v1
        device_data.v2 .= v2

        v1_t, v2_t, u1_t, u2_t = cache.v1_t, cache.v2_t, cache.u1_t, cache.u2_t
    else
        # v_t = (λ+2μ)u_xx
        v1_t = zeros(mx * my)
        v2_t = zeros(mx * my)
    end

    # ---- For constant material parameters ----
    # u1_x, u1_y, u2_x, u2_y, u1_xx, u1_yy, u2_xx, u2_yy, u2_xy, u1_yx = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material)
    # v1_t .= iD .* ((λ+2μ) .* u1_xx .+ μ .* u1_yy .+ (λ+μ) .* u2_xy)
    # v2_t .= iD .* (μ .* u2_xx .+ (λ+2μ) .* u2_yy .+ (λ+μ) .* u1_yx)
    # ---- Variable material parameters, narrow D2 on non-mixed derivatives ----
    # NOTE: For the periodic case, we can get all derivatives in one go easily, since we don't have any closures
    if Params.PERIODIC
        u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy, u1_xλy, u1_yμx = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)
        # Check if derivatives are on GPU (if compute_derivatives returned GPU arrays)
        if Params.USE_CUDA && CUDA.functional()
            # Use GPU arrays for computation
            v1_t .= device_data.iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
            v2_t .= device_data.iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
        else
            # Use CPU arrays
            v1_t .= iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
            v2_t .= iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
        end
    else
        # TODO: Support for upwind and mixed-upwind with non-periodic BC
        if Params.SCHEMA == :narrow
            # Profile.clear()
            u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy, u1_xλy, u1_yμx = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)
            # pprof()
            if Params.USE_CUDA && CUDA.functional()
                v1_t .= device_data.iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
                v2_t .= device_data.iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
            else
                v1_t .= iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
                v2_t .= iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
            end
        elseif Params.SCHEMA == :wide
            # ---- For variable material parameters ----
            u1_x, u1_y, u2_x, u2_y, _, _, _, _, _, _ = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)

            # TODO: Fix real solution for the shared memory problem, copy not optimal
            if Params.USE_CUDA && CUDA.functional()
                u1_x = copy(u1_x)
                u1_y = copy(u1_y)
                u2_x = copy(u2_x)
                u2_y = copy(u2_y)
            end

            if Params.USE_CUDA && CUDA.functional()
                λ_local = device_data.λ
                μ_local = device_data.μ
            else
                λ_local = λ
                μ_local = μ
            end
            sigma_11 = λ_local .* (u1_x .+ u2_y) .+ 2.0f0 .* μ_local .* u1_x
            sigma_22 = λ_local .* (u1_x .+ u2_y) .+ 2.0f0 .* μ_local .* u2_y
            sigma_12 = μ_local .* (u1_y .+ u2_x)

            # TODO: Minimize computational cost
            sigma_11_x, sigma_11_y, _, _, _, _, _, _, _, _ = compute_derivatives(sigma_11, sigma_11, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)
            if Params.USE_CUDA && CUDA.functional()
                sigma_11_x = copy(sigma_11_x)
                sigma_11_y = copy(sigma_11_y)
            end
            sigma_12_x, sigma_12_y, _, _, _, _, _, _, _, _ = compute_derivatives(sigma_12, sigma_12, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)
            if Params.USE_CUDA && CUDA.functional()
                sigma_12_x = copy(sigma_12_x)
                sigma_12_y = copy(sigma_12_y)
            end
            sigma_22_x, sigma_22_y, _, _, _, _, _, _, _, _ = compute_derivatives(sigma_22, sigma_22, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)

            if Params.USE_CUDA && CUDA.functional()
                v1_t .= device_data.iD .* (sigma_11_x .+ sigma_12_y)
                v2_t .= device_data.iD .* (sigma_12_x .+ sigma_22_y)
            else
                v1_t .= iD .* (sigma_11_x .+ sigma_12_y)
                v2_t .= iD .* (sigma_12_x .+ sigma_22_y)
            end
        elseif Params.SCHEMA == :mixed
            # Wide for u1_xλx, u2_yλy, narrow on rest (except mixed derivatives u2_yλx, u2_xμy, u1_xλy, u1_yμx)
            # Narrow on u1_xμx, u1_yμy, u2_yμy, u2_xμx

            if Params.USE_CUDA && CUDA.functional()
                λ_local = device_data.λ
                μ_local = device_data.μ
            else
                λ_local = λ
                μ_local = μ
            end

            u1_x, u1_y, u2_x, u2_y, _, u1_xμx, u1_yμy, _, u2_yμy, u2_xμx, u2_yλx, u2_xμy, u1_xλy, u1_yμx = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)
            if Params.USE_CUDA && CUDA.functional()
                u1_x = copy(u1_x)
                u1_y = copy(u1_y)
                u2_x = copy(u2_x)
                u2_y = copy(u2_y)
                u1_xμx = copy(u1_xμx)
                u1_yμy = copy(u1_yμy)
                u2_yμy = copy(u2_yμy)
                u2_xμx = copy(u2_xμx)
                u2_yλx = copy(u2_yλx)
                u2_xμy = copy(u2_xμy)
                u1_xλy = copy(u1_xλy)
                u1_yμx = copy(u1_yμx)
            end
            # NOTE: :wide used here
            u1_xλx, _, _, u2_yλy, _, _, _, _, _, _ = compute_derivatives(λ_local .* u1_x, λ_local .* u2_y, :wide, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)

            if Params.USE_CUDA && CUDA.functional()
                v1_t .= device_data.iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
                v2_t .= device_data.iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
            else
                v1_t .= iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
                v2_t .= iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
            end
        elseif Params.SCHEMA == :upwind
            # ---- For variable material parameters ----
            u1_x, u1_y, u2_x, u2_y, _, _, _, _, _, _ = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils; plus_or_minus=:minus)

            if Params.USE_CUDA && CUDA.functional()
                u1_x = copy(u1_x)
                u1_y = copy(u1_y)
                u2_x = copy(u2_x)
                u2_y = copy(u2_y)
            end

            if Params.USE_CUDA && CUDA.functional()
                λ_local = device_data.λ
                μ_local = device_data.μ
            else
                λ_local = λ
                μ_local = μ
            end

            sigma_11 = λ_local .* (u1_x .+ u2_y) .+ 2.0f0 .* μ_local .* u1_x
            sigma_22 = λ_local .* (u1_x .+ u2_y) .+ 2.0f0 .* μ_local .* u2_y
            sigma_12 = μ_local .* (u1_y .+ u2_x)

            # TODO: Minimize computational cost
            sigma_11_x, sigma_11_y, _, _, _, _, _, _, _, _ = compute_derivatives(sigma_11, sigma_11, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils; plus_or_minus=:plus)
            if Params.USE_CUDA && CUDA.functional()
                sigma_11_x = copy(sigma_11_x)
                sigma_11_y = copy(sigma_11_y)
            end
            sigma_12_x, sigma_12_y, _, _, _, _, _, _, _, _ = compute_derivatives(sigma_12, sigma_12, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils; plus_or_minus=:plus)
            if Params.USE_CUDA && CUDA.functional()
                sigma_12_x = copy(sigma_12_x)
                sigma_12_y = copy(sigma_12_y)
            end
            sigma_22_x, sigma_22_y, _, _, _, _, _, _, _, _ = compute_derivatives(sigma_22, sigma_22, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils; plus_or_minus=:plus)

            if Params.USE_CUDA && CUDA.functional()
                v1_t .= device_data.iD .* (sigma_11_x .+ sigma_12_y)
                v2_t .= device_data.iD .* (sigma_12_x .+ sigma_22_y)
            else
                v1_t .= iD .* (sigma_11_x .+ sigma_12_y)
                v2_t .= iD .* (sigma_12_x .+ sigma_22_y)
            end
        elseif Params.SCHEMA == :mixed_upwind_μ
            # # Wide (where D1 = D+ + D- / 2) for u1_xλx, u2_yλy, u2_yλx, u1_xλy
            # # Narrow on u1_xμx, u1_yμy, u2_yμy, u2_xμx
            # # Upwind on u2_xμy, u1_yμx

            # # Narrow
            # u1_x, u1_y, u2_x, u2_y, _, u1_xμx, u1_yμy, _, u2_yμy, u2_xμx, _, _, _, _ = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data)
            # # Wide
            # u1_xλx, _, _, u2_yλy, _, _, _, _, _, _, _, _, _, _ = compute_derivatives(λ .* u1_x, λ .* u2_y, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data; skip_narrow=true)
            # u1_xλy, _, _, u2_yλx, _, _, _, _, _, _, _, _, _, _ = compute_derivatives(λ .* u1_y, λ .* u2_x, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data; skip_narrow=true)
            # # Upwind
            # u1_x_u, u1_y_u, u2_x_u, u2_y_u, _, _, _, _, _, _ = compute_derivatives(u1, u2, :upwind, Params.PERIODIC, grid, stencils, material, derivatives, device_data; plus_or_minus=:minus)
            # u2_xμy, _, _, u1_yμx, _, _, _, _, _, _ = compute_derivatives(μ .* u2_y_u, μ .* u1_x_u, :upwind, Params.PERIODIC, grid, stencils, material, derivatives, device_data; plus_or_minus=:plus)

            # v1_t .= iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
            # v2_t .= iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
            error("mixed_upwind_μ not implemented yet")
        elseif Params.SCHEMA == :mixed_upwind_λ
            # Wide (where D1 = D+ + D- / 2) for u2_xμy, u1_yμx
            # Narrow on u1_xμx, u1_yμy, u2_yμy, u2_xμx
            # Upwind on u1_xλx, u2_yλy, u2_yλx, u1_xλy
            # Narrow (wide on u2_xμy, u1_yμx)
            u1_x, u1_y, u2_x, u2_y, _, u1_xμx, u1_yμy, _, u2_yμy, u2_xμx, _, u2_xμy, _, u1_yμx = compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils)
            # Upwind
            u1_x_u, u1_y_u, u2_x_u, u2_y_u, _, _, _, _, _, _ = compute_derivatives(u1, u2, :upwind, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils; plus_or_minus=:minus)
            u1_xλx, _, u2_yλx, _, _, _, _, _, _, _ = compute_derivatives(λ .* u1_x_u, λ .* u2_y_u, :upwind, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils; plus_or_minus=:plus)
            _, u1_xλy, _, u2_yλy, _, _, _, _, _, _ = compute_derivatives(λ .* u1_x_u, λ .* u2_y_u, :upwind, Params.PERIODIC, grid, stencils, material, derivatives, device_data, D2_stencils; plus_or_minus=:plus)

            v1_t .= iD .* (u1_xλx .+ 2.0f0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
            v2_t .= iD .* (u2_xμx .+ u2_yλy .+ 2.0f0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
            
            # Might be useful
            u1_y_narrow, u2_y_upwind = u1_y, u2_y_u  # For west/east boundaries
            u1_x_upwind, u2_x_narrow = u1_x_u, u2_x  # For south/north boundaries
        else
            error("Unknown Params.SCHEMA type")
        end
    end

    if Params.PULSE
        # Build dirac distribution function
        if Params.USE_CUDA && CUDA.functional()
            v1_t .+= device_data.delta_2D * sourceTime(t, Params.T0, Params.SIGMA)
            v2_t .+= device_data.delta_2D * sourceTime(t, Params.T0, Params.SIGMA)
        else
            quadrature = [stencils.Hx, stencils.Hy]
            delta_2D = diracDiscr(grid, Params.X_S, Params.M_ORDER, Params.S_ORDER, quadrature)
            v1_t += sourceTime(t, Params.T0, Params.SIGMA) * delta_2D
            v2_t += sourceTime(t, Params.T0, Params.SIGMA) * delta_2D
        end
    end

    if Params.USE_MMS
        # Add body force from MMS
        if Params.USE_CUDA && CUDA.functional()
            fx_vec = zeros(mx * my)
            fy_vec = zeros(mx * my)
            for j in 1:my
                for i in 1:mx
                    fx_vec[idx(i, j, mx, my)] = mms.fx(x_values[i], y_values[j], t, ω)
                    fy_vec[idx(i, j, mx, my)] = mms.fy(x_values[i], y_values[j], t, ω)
                end
            end

            # NOTE: Slow!
            v1_t .+= CuArray(fx_vec)
            v2_t .+= CuArray(fy_vec)
        else
            for j in 1:my
                for i in 1:mx
                    v1_t[idx(i, j, mx, my)] += mms.fx(x_values[i], y_values[j], t, ω)
                    v2_t[idx(i, j, mx, my)] += mms.fy(x_values[i], y_values[j], t, ω)
                end
            end
        end
    end

    if !Params.PERIODIC
        if Params.USE_CUDA && CUDA.functional()
            d1x_l, d1x_r = device_data.d1x_l, device_data.d1x_r
            d1y_l, d1y_r = device_data.d1y_l, device_data.d1y_r
        else
            if Params.SCHEMA in (:mixed_upwind_λ, :mixed_upwind_μ)
                d1x_l, d1x_r = stencils.d1x_l_narrow, stencils.d1x_r_narrow
                d1y_l, d1y_r = stencils.d1y_l_narrow, stencils.d1y_r_narrow
                # println("d1x_l is narrow")
            else
                d1x_l, d1x_r = stencils.d1x_l, stencils.d1x_r
                d1y_l, d1y_r = stencils.d1y_l, stencils.d1y_r
            end
        end
        if Params.SCHEMA == :upwind
            # No narrow derivatives for pure upwind
            d1x_l_wide, d1x_r_wide = d1x_l, d1x_r
            d1y_l_wide, d1y_r_wide = d1y_l, d1y_r
        elseif Params.SCHEMA in (:mixed_upwind_λ, :mixed_upwind_μ)
            d1x_l_wide, d1x_r_wide = stencils.d1x_l_upwind, stencils.d1x_r_upwind
            d1y_l_wide, d1y_r_wide = stencils.d1y_l_upwind, stencils.d1y_r_upwind
            # println("d1x_l_wide is upwind")
        # elseif Params.SCHEMA in (:mixed_upwind_λ, :mixed_upwind_μ)
        #     # # For mixed, we utilize D1 = D+ + D- / 2, therefore use upwind norm
        #     # d1x_l_wide, d1x_r_wide = stencils.d1x_l_central, stencils.d1x_r_central
        #     # d1y_l_wide, d1y_r_wide = stencils.d1y_l_central, stencils.d1y_r_central
        #     error("mixed_upwind_λ and mixed_upwind_μ not implemented yet")
        else
            if Params.USE_CUDA && CUDA.functional()
                d1x_l_wide, d1x_r_wide = device_data.d1x_l_wide, device_data.d1x_r_wide
                d1y_l_wide, d1y_r_wide = device_data.d1y_l_wide, device_data.d1y_r_wide
            else
                d1x_l_wide, d1x_r_wide = stencils.d1x_l_wide, stencils.d1x_r_wide
                d1y_l_wide, d1y_r_wide = stencils.d1y_l_wide, stencils.d1y_r_wide
            end
        end

        if Params.SCHEMA == :wide
            # Override narrow derivatives with wide ones
            d1x_l, d1x_r = d1x_l_wide, d1x_r_wide
            d1y_l, d1y_r = d1y_l_wide, d1y_r_wide
        end

        # Neumann BC SAT terms
        tau_W1 = -1.0f0
        tau_W2 = -1.0f0
        tau_E1 = -1.0f0
        tau_E2 = -1.0f0
        tau_S1 = -1.0f0
        tau_S2 = -1.0f0
        tau_N1 = -1.0f0
        tau_N2 = -1.0f0

        d_radius = max(length(d1x_l), length(d1x_r), length(d1y_l), length(d1y_r),
            length(d1x_l_wide), length(d1x_r_wide), length(d1y_l_wide), length(d1y_r_wide))

        if Params.USE_CUDA && CUDA.functional()
            threads = 256
            blocks = cld(my, threads)
            schema_is_mixed = (Params.SCHEMA == :mixed || Params.SCHEMA == :mixed_upwind_λ)
        end

        # Can look into changing these around
        # u1_y_for_west_east = u1_y
        # u2_y_for_west_east = u2_y
        # u1_x_for_south_north = u1_x
        # u2_x_for_south_north = u2_x
        # E.g...
        if Params.SCHEMA == :mixed_upwind_λ
            u1_y_for_west_east = u1_y      # narrow (correct for μ)
            u2_y_for_west_east = u2_y_u    # upwind (needed for λ at west/east)
            u1_x_for_south_north = u1_x_u  # upwind (needed for λ at south/north)
            u2_x_for_south_north = u2_x    # narrow (correct for μ)
        else
            u1_y_for_west_east = u1_y
            u2_y_for_west_east = u2_y
            u1_x_for_south_north = u1_x
            u2_x_for_south_north = u2_x
        end

        if bcs.west == :traction_free
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_traction_free_west_kernel!(
                    v1_t, v2_t, u1, u2, u1_y, u2_y, d1x_l, d1x_l_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.iD, quad_term_x,
                    tau_W1, tau_W2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_traction_free_west!(Params.SCHEMA, v1_t, v2_t, u1, u2, u1_y_for_west_east, u2_y_for_west_east, mx, my,
                    d1x_l, d1x_l_wide, d_radius, λ, μ, iD, quad_term_x, tau_W1, tau_W2)
            end
        elseif bcs.west == :non_reflecting
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_non_reflecting_west_kernel!(
                    v1_t, v2_t, v1, v2, u1, u2, u1_y, u2_y, d1x_l, d1x_l_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.D, device_data.iD, quad_term_x,
                    tau_W1, tau_W2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_non_reflecting_west!(Params.SCHEMA, v1_t, v2_t, v1, v2, u1, u2, u1_y_for_west_east, u2_y_for_west_east, mx, my,
                    d1x_l, d1x_l_wide, d_radius, λ, μ, D, iD, quad_term_x, tau_W1, tau_W2)
            end
        else
            error("Unknown BC type for west boundary")
        end

        if bcs.east == :traction_free
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_traction_free_east_kernel!(
                    v1_t, v2_t, u1, u2, u1_y, u2_y, d1x_r, d1x_r_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.iD, quad_term_x,
                    tau_E1, tau_E2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_traction_free_east!(Params.SCHEMA, v1_t, v2_t, u1, u2, u1_y_for_west_east, u2_y_for_west_east, mx, my,
                    d1x_r, d1x_r_wide, d_radius, λ, μ, iD, quad_term_x, tau_E1, tau_E2)
            end
        elseif bcs.east == :non_reflecting
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_non_reflecting_east_kernel!(
                    v1_t, v2_t, v1, v2, u1, u2, u1_y, u2_y, d1x_r, d1x_r_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.D, device_data.iD, quad_term_x,
                    tau_E1, tau_E2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_non_reflecting_east!(Params.SCHEMA, v1_t, v2_t, v1, v2, u1, u2, u1_y_for_west_east, u2_y_for_west_east, mx, my,
                    d1x_r, d1x_r_wide, d_radius, λ, μ, D, iD, quad_term_x, tau_E1, tau_E2)
            end
        else
            error("Unknown BC type for east boundary")
        end

        if bcs.south == :traction_free
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_traction_free_south_kernel!(
                    v1_t, v2_t, u1, u2, u1_x, u2_x, d1y_l, d1y_l_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.iD, quad_term_y,
                    tau_S1, tau_S2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_traction_free_south!(Params.SCHEMA, v1_t, v2_t, u1, u2, u1_x_for_south_north, u2_x_for_south_north, mx, my,
                    d1y_l, d1y_l_wide, d_radius, λ, μ, iD, quad_term_y, tau_S1, tau_S2)
            end
        elseif bcs.south == :non_reflecting
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_non_reflecting_south_kernel!(
                    v1_t, v2_t, v1, v2, u1, u2, u1_x, u2_x, d1y_l, d1y_l_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.D, device_data.iD, quad_term_y,
                    tau_S1, tau_S2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_non_reflecting_south!(Params.SCHEMA, v1_t, v2_t, v1, v2, u1, u2, u1_x_for_south_north, u2_x_for_south_north, mx, my,
                    d1y_l, d1y_l_wide, d_radius, λ, μ, D, iD, quad_term_y, tau_S1, tau_S2)
            end
        else
            error("Unknown BC type for south boundary")
        end

        if bcs.north == :traction_free
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_traction_free_north_kernel!(
                    v1_t, v2_t, u1, u2, u1_x, u2_x, d1y_r, d1y_r_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.iD, quad_term_y,
                    tau_N1, tau_N2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_traction_free_north!(Params.SCHEMA, v1_t, v2_t, u1, u2, u1_x_for_south_north, u2_x_for_south_north, mx, my,
                    d1y_r, d1y_r_wide, d_radius, λ, μ, iD, quad_term_y, tau_N1, tau_N2)
            end
        elseif bcs.north == :non_reflecting
            if Params.USE_CUDA && CUDA.functional()
                @cuda threads = threads blocks = blocks apply_non_reflecting_north_kernel!(
                    v1_t, v2_t, v1, v2, u1, u2, u1_x, u2_x, d1y_r, d1y_r_wide, d_radius,
                    device_data.λ, device_data.μ, device_data.D, device_data.iD, quad_term_y,
                    tau_N1, tau_N2, mx, my, schema_is_mixed
                )
                # CUDA.synchronize()
            else
                apply_non_reflecting_north!(Params.SCHEMA, v1_t, v2_t, v1, v2, u1, u2, u1_x_for_south_north, u2_x_for_south_north, mx, my,
                    d1y_r, d1y_r_wide, d_radius, λ, μ, D, iD, quad_term_y, tau_N1, tau_N2)
            end
        else
            error("Unknown BC type for north boundary")
        end

    end

    if Params.USE_CUDA && CUDA.functional()
        u1_t .= v1
        u2_t .= v2
    else
        u1_t = copy(v1)
        u2_t = copy(v2)
    end

    return v1_t, v2_t, u1_t, u2_t
end
