@inline function apply_D_periodic!(SCHEMA, u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy,
                                   u1_xλy, u1_yμx, u1, u2, i, j, mx, my, grid, stencils, material, rowλ_x, rowμ_x, rowλ_y, rowμ_y)
    Ω_radius = stencils.Ω_radius # k max for interior
    # interior_offset = stencils.SBP_ORDER == 2 ? 1 : 2 # Or floor(Ω_radius/2)
    interior_offset = Ω_radius ÷ 2

    if SCHEMA == :narrow
        # D2x_∂Ω only sees x direction, D2y_∂Ω only y direction
        extract_row_x!(rowλ_x, material.λ, j, mx, my)
        extract_row_x!(rowμ_x, material.μ, j, mx, my)
        extract_col_y!(rowλ_y, material.λ, i, mx, my)
        extract_col_y!(rowμ_y, material.μ, i, mx, my)


        # Interior D2
        D2x_λ = D2_inner_periodic(rowλ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2x_μ = D2_inner_periodic(rowμ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2y_μ = D2_inner_periodic(rowμ_y, j, my, grid.hy, stencils.SBP_ORDER)
        D2y_λ = D2_inner_periodic(rowλ_y, j, my, grid.hy, stencils.SBP_ORDER)
        for k in 1:Ω_radius
            wrapped_x = mod1(i - interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - interior_offset - 1 + k, my)
            # -- D1 --
            u1_x[idx(i, j, mx, my)] += stencils.D1x[k] * u1[idx(wrapped_x, j, mx, my)]
            u2_x[idx(i, j, mx, my)] += stencils.D1x[k] * u2[idx(wrapped_x, j, mx, my)]

            u1_y[idx(i, j, mx, my)] += stencils.D1y[k] * u1[idx(i, wrapped_y, mx, my)]
            u2_y[idx(i, j, mx, my)] += stencils.D1y[k] * u2[idx(i, wrapped_y, mx, my)]

            # -- D2, narrow --
            u1_xλx[idx(i, j, mx, my)] += D2x_λ[k] * u1[idx(wrapped_x, j, mx, my)]
            u1_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u1[idx(wrapped_x, j, mx, my)]
            u2_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u2[idx(wrapped_x, j, mx, my)]

            u1_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u1[idx(i, wrapped_y, mx, my)]
            u2_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u2[idx(i, wrapped_y, mx, my)]
            u2_yλy[idx(i, j, mx, my)] += D2y_λ[k] * u2[idx(i, wrapped_y, mx, my)]
        end
    elseif SCHEMA == :wide
        # -- D1 --
        for k in 1:Ω_radius
            wrapped_x = mod1(i - interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - interior_offset - 1 + k, my)
            # -- D1 --
            u1_x[idx(i, j, mx, my)] += stencils.D1x[k] * u1[idx(wrapped_x, j, mx, my)]
            u2_x[idx(i, j, mx, my)] += stencils.D1x[k] * u2[idx(wrapped_x, j, mx, my)]

            u1_y[idx(i, j, mx, my)] += stencils.D1y[k] * u1[idx(i, wrapped_y, mx, my)]
            u2_y[idx(i, j, mx, my)] += stencils.D1y[k] * u2[idx(i, wrapped_y, mx, my)]
        end
        # -- D2, wide (apply D1 twice with material coefficients in between) --
        # x-direction: u1_xλx = ∂x(λ ∂x u1), u1_xμx = ∂x(μ ∂x u1), u2_xμx = ∂x(μ ∂x u2)
        for kx in 1:Ω_radius
            wrapped_kx = mod1(i - interior_offset - 1 + kx, mx)

            # Compute inner derivative: ∂x u at position wrapped_kx
            inner_u1_x = 0.0
            inner_u2_x = 0.0
            @inbounds for m in 1:Ω_radius
                wrapped_m = mod1(wrapped_kx - interior_offset - 1 + m, mx)
                inner_u1_x += stencils.D1x[m] * u1[idx(wrapped_m, j, mx, my)]
                inner_u2_x += stencils.D1x[m] * u2[idx(wrapped_m, j, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.D1x[kx]
            λkx = material.λ[idx(wrapped_kx, j, mx, my)]
            μkx = material.μ[idx(wrapped_kx, j, mx, my)]

            u1_xλx[idx(i, j, mx, my)] += w * λkx * inner_u1_x
            u1_xμx[idx(i, j, mx, my)] += w * μkx * inner_u1_x
            u2_xμx[idx(i, j, mx, my)] += w * μkx * inner_u2_x
        end

        # y-direction: u1_yμy = ∂y(μ ∂y u1), u2_yμy = ∂y(μ ∂y u2), u2_yλy = ∂y(λ ∂y u2)
        for ky in 1:Ω_radius
            wrapped_ky = mod1(j - interior_offset - 1 + ky, my)

            # Compute inner derivative: ∂y u at position wrapped_ky
            inner_u1_y = 0.0
            inner_u2_y = 0.0
            @inbounds for n in 1:Ω_radius
                wrapped_n = mod1(wrapped_ky - interior_offset - 1 + n, my)
                inner_u1_y += stencils.D1y[n] * u1[idx(i, wrapped_n, mx, my)]
                inner_u2_y += stencils.D1y[n] * u2[idx(i, wrapped_n, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.D1y[ky]
            λky = material.λ[idx(i, wrapped_ky, mx, my)]
            μky = material.μ[idx(i, wrapped_ky, mx, my)]

            u1_yμy[idx(i, j, mx, my)] += w * μky * inner_u1_y
            u2_yμy[idx(i, j, mx, my)] += w * μky * inner_u2_y
            u2_yλy[idx(i, j, mx, my)] += w * λky * inner_u2_y
        end
    elseif SCHEMA == :mixed
        # D2x_∂Ω only sees x direction, D2y_∂Ω only y direction
        extract_row_x!(rowλ_x, material.λ, j, mx, my)
        extract_row_x!(rowμ_x, material.μ, j, mx, my)
        extract_col_y!(rowλ_y, material.λ, i, mx, my)
        extract_col_y!(rowμ_y, material.μ, i, mx, my)


        # Interior D2
        D2x_λ = D2_inner_periodic(rowλ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2x_μ = D2_inner_periodic(rowμ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2y_μ = D2_inner_periodic(rowμ_y, j, my, grid.hy, stencils.SBP_ORDER)
        D2y_λ = D2_inner_periodic(rowλ_y, j, my, grid.hy, stencils.SBP_ORDER)
        # -- D2, mixed i.e. wide for λ and narrow for μ --
        for k in 1:Ω_radius
            wrapped_x = mod1(i - interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - interior_offset - 1 + k, my)
            # -- D1 --
            u1_x[idx(i, j, mx, my)] += stencils.D1x[k] * u1[idx(wrapped_x, j, mx, my)]
            u2_x[idx(i, j, mx, my)] += stencils.D1x[k] * u2[idx(wrapped_x, j, mx, my)]

            u1_y[idx(i, j, mx, my)] += stencils.D1y[k] * u1[idx(i, wrapped_y, mx, my)]
            u2_y[idx(i, j, mx, my)] += stencils.D1y[k] * u2[idx(i, wrapped_y, mx, my)]

            # -- D2, narrow --
            u1_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u1[idx(wrapped_x, j, mx, my)]
            u2_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u2[idx(wrapped_x, j, mx, my)]

            u1_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u1[idx(i, wrapped_y, mx, my)]
            u2_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u2[idx(i, wrapped_y, mx, my)]
        end
        for kx in 1:Ω_radius
            wrapped_kx = mod1(i - interior_offset - 1 + kx, mx)

            # Compute inner derivative: ∂x u at position wrapped_kx
            inner_u1_x = 0.0
            @inbounds for m in 1:Ω_radius
                wrapped_m = mod1(wrapped_kx - interior_offset - 1 + m, mx)
                inner_u1_x += stencils.D1x[m] * u1[idx(wrapped_m, j, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.D1x[kx]
            λkx = material.λ[idx(wrapped_kx, j, mx, my)]

            u1_xλx[idx(i, j, mx, my)] += w * λkx * inner_u1_x
        end

        for ky in 1:Ω_radius
            wrapped_ky = mod1(j - interior_offset - 1 + ky, my)

            # Compute inner derivative: ∂y u at position wrapped_ky
            inner_u2_y = 0.0
            @inbounds for n in 1:Ω_radius
                wrapped_n = mod1(wrapped_ky - interior_offset - 1 + n, my)
                inner_u2_y += stencils.D1y[n] * u2[idx(i, wrapped_n, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.D1y[ky]
            λky = material.λ[idx(i, wrapped_ky, mx, my)]

            u2_yλy[idx(i, j, mx, my)] += w * λky * inner_u2_y
        end
    elseif SCHEMA == :upwind
        # -- D1 with centered stencil (average of Dp and Dm) --
        for k in 1:Ω_radius
            wrapped_x = mod1(i - interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - interior_offset - 1 + k, my)
            # Use average of upwind operators to get centered derivative
            D1x_centered = 0.5 * (stencils.Dpx[k] + stencils.Dmx[k])
            D1y_centered = 0.5 * (stencils.Dpy[k] + stencils.Dmy[k])

            u1_x[idx(i, j, mx, my)] += D1x_centered * u1[idx(wrapped_x, j, mx, my)]
            u2_x[idx(i, j, mx, my)] += D1x_centered * u2[idx(wrapped_x, j, mx, my)]

            u1_y[idx(i, j, mx, my)] += D1y_centered * u1[idx(i, wrapped_y, mx, my)]
            u2_y[idx(i, j, mx, my)] += D1y_centered * u2[idx(i, wrapped_y, mx, my)]
        end

        # -- D2, wide: Use Dm(c * Dp(u)) for SBP property --
        # x-direction: u1_xλx = ∂x(λ ∂x u1), u1_xμx = ∂x(μ ∂x u1), u2_xμx = ∂x(μ ∂x u2)
        for kx in 1:Ω_radius
            wrapped_kx = mod1(i - interior_offset - 1 + kx, mx)

            # Compute inner derivative: ∂x u at position wrapped_kx
            inner_u1_x = 0.0
            inner_u2_x = 0.0
            @inbounds for m in 1:Ω_radius
                wrapped_m = mod1(wrapped_kx - interior_offset - 1 + m, mx)
                inner_u1_x += stencils.Dmx[m] * u1[idx(wrapped_m, j, mx, my)]
                inner_u2_x += stencils.Dmx[m] * u2[idx(wrapped_m, j, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.Dpx[kx]
            λkx = material.λ[idx(wrapped_kx, j, mx, my)]
            μkx = material.μ[idx(wrapped_kx, j, mx, my)]

            u1_xλx[idx(i, j, mx, my)] += w * λkx * inner_u1_x
            u1_xμx[idx(i, j, mx, my)] += w * μkx * inner_u1_x
            u2_xμx[idx(i, j, mx, my)] += w * μkx * inner_u2_x
        end

        # y-direction: u1_yμy = ∂y(μ ∂y u1), u2_yμy = ∂y(μ ∂y u2), u2_yλy = ∂y(λ ∂y u2)
        for ky in 1:Ω_radius
            wrapped_ky = mod1(j - interior_offset - 1 + ky, my)

            # Compute inner derivative: ∂y u at position wrapped_ky
            inner_u1_y = 0.0
            inner_u2_y = 0.0
            @inbounds for n in 1:Ω_radius
                wrapped_n = mod1(wrapped_ky - interior_offset - 1 + n, my)
                inner_u1_y += stencils.Dmy[n] * u1[idx(i, wrapped_n, mx, my)]
                inner_u2_y += stencils.Dmy[n] * u2[idx(i, wrapped_n, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.Dpy[ky]
            λky = material.λ[idx(i, wrapped_ky, mx, my)]
            μky = material.μ[idx(i, wrapped_ky, mx, my)]

            u1_yμy[idx(i, j, mx, my)] += w * μky * inner_u1_y
            u2_yμy[idx(i, j, mx, my)] += w * μky * inner_u2_y
            u2_yλy[idx(i, j, mx, my)] += w * λky * inner_u2_y
        end
    elseif SCHEMA == :mixed_upwind_μ
        narrow_Ω_radius = stencils.narrow_Ω_radius
        narrow_interior_offset = narrow_Ω_radius ÷ 2
        # D2x_∂Ω only sees x direction, D2y_∂Ω only y direction
        extract_row_x!(rowλ_x, material.λ, j, mx, my)
        extract_row_x!(rowμ_x, material.μ, j, mx, my)
        extract_col_y!(rowλ_y, material.λ, i, mx, my)
        extract_col_y!(rowμ_y, material.μ, i, mx, my)


        # Interior D2
        D2x_λ = D2_inner_periodic(rowλ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2x_μ = D2_inner_periodic(rowμ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2y_μ = D2_inner_periodic(rowμ_y, j, my, grid.hy, stencils.SBP_ORDER)
        D2y_λ = D2_inner_periodic(rowλ_y, j, my, grid.hy, stencils.SBP_ORDER)
        # -- D2, mixed i.e. wide for λ and narrow for μ --
        for k in 1:Ω_radius
            wrapped_x = mod1(i - interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - interior_offset - 1 + k, my)
            # -- D1 --
            D1x_centered = 0.5 * (stencils.Dpx[k] + stencils.Dmx[k])
            D1y_centered = 0.5 * (stencils.Dpy[k] + stencils.Dmy[k])
            # NOTE: Can do it like this instead
            # D1x = [1/12, -2/3, 0, 2/3, -1/12]*1/grid.hx
            # D1y = [1/12, -2/3, 0, 2/3, -1/12]*1/grid.hy
            # D1x_centered = D1x[k]
            # D1y_centered = D1y[k]

            u1_x[idx(i, j, mx, my)] += D1x_centered * u1[idx(wrapped_x, j, mx, my)]
            u2_x[idx(i, j, mx, my)] += D1x_centered * u2[idx(wrapped_x, j, mx, my)]

            u1_y[idx(i, j, mx, my)] += D1y_centered * u1[idx(i, wrapped_y, mx, my)]
            u2_y[idx(i, j, mx, my)] += D1y_centered * u2[idx(i, wrapped_y, mx, my)]
        end
        for k in 1:narrow_Ω_radius
            wrapped_x = mod1(i - narrow_interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - narrow_interior_offset - 1 + k, my)

            # -- D2, narrow --
            u1_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u1[idx(wrapped_x, j, mx, my)]
            u2_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u2[idx(wrapped_x, j, mx, my)]

            u1_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u1[idx(i, wrapped_y, mx, my)]
            u2_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u2[idx(i, wrapped_y, mx, my)]
        end

        # Compute centered stencils from upwind operators (once, outside loops)
        D1x_centered = 0.5 .* (stencils.Dpx .+ stencils.Dmx)
        D1y_centered = 0.5 .* (stencils.Dpy .+ stencils.Dmy)

        # Wide derivatives with centered stencils
        for kx in 1:Ω_radius
            wrapped_kx = mod1(i - interior_offset - 1 + kx, mx)

            # Compute inner derivative: ∂x u at position wrapped_kx
            inner_u1_x = 0.0
            @inbounds for m in 1:Ω_radius
                wrapped_m = mod1(wrapped_kx - interior_offset - 1 + m, mx)
                inner_u1_x += D1x_centered[m] * u1[idx(wrapped_m, j, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = D1x_centered[kx]
            λkx = material.λ[idx(wrapped_kx, j, mx, my)]

            u1_xλx[idx(i, j, mx, my)] += w * λkx * inner_u1_x
        end

        for ky in 1:Ω_radius
            wrapped_ky = mod1(j - interior_offset - 1 + ky, my)

            # Compute inner derivative: ∂y u at position wrapped_ky
            inner_u2_y = 0.0
            @inbounds for n in 1:Ω_radius
                wrapped_n = mod1(wrapped_ky - interior_offset - 1 + n, my)
                inner_u2_y += D1y_centered[n] * u2[idx(i, wrapped_n, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = D1y_centered[ky]
            λky = material.λ[idx(i, wrapped_ky, mx, my)]

            u2_yλy[idx(i, j, mx, my)] += w * λky * inner_u2_y
        end
    elseif SCHEMA == :mixed_upwind_λ
        narrow_Ω_radius = stencils.narrow_Ω_radius
        narrow_interior_offset = narrow_Ω_radius ÷ 2
        # D2x_∂Ω only sees x direction, D2y_∂Ω only y direction
        extract_row_x!(rowλ_x, material.λ, j, mx, my)
        extract_row_x!(rowμ_x, material.μ, j, mx, my)
        extract_col_y!(rowλ_y, material.λ, i, mx, my)
        extract_col_y!(rowμ_y, material.μ, i, mx, my)

        # Interior D2
        D2x_λ = D2_inner_periodic(rowλ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2x_μ = D2_inner_periodic(rowμ_x, i, mx, grid.hx, stencils.SBP_ORDER)
        D2y_μ = D2_inner_periodic(rowμ_y, j, my, grid.hy, stencils.SBP_ORDER)
        D2y_λ = D2_inner_periodic(rowλ_y, j, my, grid.hy, stencils.SBP_ORDER)
        # -- D2, mixed i.e. wide for λ and narrow for μ --
        for k in 1:Ω_radius
            wrapped_x = mod1(i - interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - interior_offset - 1 + k, my)
            # -- D1 --
            D1x_centered = 0.5 * (stencils.Dpx[k] + stencils.Dmx[k])
            D1y_centered = 0.5 * (stencils.Dpy[k] + stencils.Dmy[k])

            u1_x[idx(i, j, mx, my)] += D1x_centered * u1[idx(wrapped_x, j, mx, my)]
            u2_x[idx(i, j, mx, my)] += D1x_centered * u2[idx(wrapped_x, j, mx, my)]

            u1_y[idx(i, j, mx, my)] += D1y_centered * u1[idx(i, wrapped_y, mx, my)]
            u2_y[idx(i, j, mx, my)] += D1y_centered * u2[idx(i, wrapped_y, mx, my)]
        end
        for k in 1:narrow_Ω_radius
            wrapped_x = mod1(i - narrow_interior_offset - 1 + k, mx)
            wrapped_y = mod1(j - narrow_interior_offset - 1 + k, my)

            # -- D2, narrow --
            u1_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u1[idx(wrapped_x, j, mx, my)]
            u2_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u2[idx(wrapped_x, j, mx, my)]

            u1_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u1[idx(i, wrapped_y, mx, my)]
            u2_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u2[idx(i, wrapped_y, mx, my)]
        end

        for kx in 1:Ω_radius
            wrapped_kx = mod1(i - interior_offset - 1 + kx, mx)

            # Compute inner derivative: ∂x u at position wrapped_kx
            inner_u1_x = 0.0
            inner_u2_x = 0.0
            @inbounds for m in 1:Ω_radius
                wrapped_m = mod1(wrapped_kx - interior_offset - 1 + m, mx)
                inner_u1_x += stencils.Dmx[m] * u1[idx(wrapped_m, j, mx, my)]
                inner_u2_x += stencils.Dmx[m] * u2[idx(wrapped_m, j, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.Dpx[kx]
            λkx = material.λ[idx(wrapped_kx, j, mx, my)]
            μkx = material.μ[idx(wrapped_kx, j, mx, my)]

            u1_xλx[idx(i, j, mx, my)] += w * λkx * inner_u1_x
            # u1_xμx[idx(i, j, mx, my)] += w * μkx * inner_u1_x
            # u2_xμx[idx(i, j, mx, my)] += w * μkx * inner_u2_x
        end

        # y-direction: u1_yμy = ∂y(μ ∂y u1), u2_yμy = ∂y(μ ∂y u2), u2_yλy = ∂y(λ ∂y u2)
        for ky in 1:Ω_radius
            wrapped_ky = mod1(j - interior_offset - 1 + ky, my)

            # Compute inner derivative: ∂y u at position wrapped_ky
            inner_u1_y = 0.0
            inner_u2_y = 0.0
            @inbounds for n in 1:Ω_radius
                wrapped_n = mod1(wrapped_ky - interior_offset - 1 + n, my)
                inner_u1_y += stencils.Dmy[n] * u1[idx(i, wrapped_n, mx, my)]
                inner_u2_y += stencils.Dmy[n] * u2[idx(i, wrapped_n, mx, my)]
            end

            # Apply outer derivative with material coefficients
            w = stencils.Dpy[ky]
            λky = material.λ[idx(i, wrapped_ky, mx, my)]
            μky = material.μ[idx(i, wrapped_ky, mx, my)]

            # u1_yμy[idx(i, j, mx, my)] += w * μky * inner_u1_y
            # u2_yμy[idx(i, j, mx, my)] += w * μky * inner_u2_y
            u2_yλy[idx(i, j, mx, my)] += w * λky * inner_u2_y
        end
    else
        error("Schema $SCHEMA not supported.")
    end

    if SCHEMA == :upwind
        # Mixed derivatives: Use CENTERED approximations (Dp + Dm)/2 for symmetry
        for kx in 1:Ω_radius
            for ky in 1:Ω_radius
                # Periodic wrapping
                wrapped_x = mod1(i - interior_offset - 1 + kx, mx)
                wrapped_y = mod1(j - interior_offset - 1 + ky, my)

                u2_yλx[idx(i, j, mx, my)] += stencils.Dpx[kx] * material.λ[idx(wrapped_x, j, mx, my)] * stencils.Dmy[ky] * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u2_xμy[idx(i, j, mx, my)] += stencils.Dpy[ky] * material.μ[idx(i, wrapped_y, mx, my)] * stencils.Dmx[kx] * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u1_xλy[idx(i, j, mx, my)] += stencils.Dpy[ky] * material.λ[idx(i, wrapped_y, mx, my)] * stencils.Dmx[kx] * u1[idx(wrapped_x, wrapped_y, mx, my)]
                u1_yμx[idx(i, j, mx, my)] += stencils.Dpx[kx] * material.μ[idx(wrapped_x, j, mx, my)] * stencils.Dmy[ky] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            end
        end
    elseif SCHEMA == :mixed_upwind_μ
        for kx in 1:Ω_radius
            for ky in 1:Ω_radius
                wrapped_x = mod1(i - interior_offset - 1 + kx, mx)
                wrapped_y = mod1(j - interior_offset - 1 + ky, my)

                # Centered stencils from averaging upwind operators
                D1x_centered = 0.5 * (stencils.Dpx[kx] + stencils.Dmx[kx])
                D1y_centered = 0.5 * (stencils.Dpy[ky] + stencils.Dmy[ky])

                u2_yλx[idx(i, j, mx, my)] += D1x_centered * material.λ[idx(wrapped_x, j, mx, my)] * D1y_centered * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u2_xμy[idx(i, j, mx, my)] += stencils.Dpy[ky] * material.μ[idx(i, wrapped_y, mx, my)] * stencils.Dmx[kx] * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u1_xλy[idx(i, j, mx, my)] += D1y_centered * material.λ[idx(i, wrapped_y, mx, my)] * D1x_centered * u1[idx(wrapped_x, wrapped_y, mx, my)]
                u1_yμx[idx(i, j, mx, my)] += stencils.Dpx[kx] * material.μ[idx(wrapped_x, j, mx, my)] * stencils.Dmy[ky] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            end
        end
    elseif SCHEMA == :mixed_upwind_λ
        for kx in 1:Ω_radius
            for ky in 1:Ω_radius
                wrapped_x = mod1(i - interior_offset - 1 + kx, mx)
                wrapped_y = mod1(j - interior_offset - 1 + ky, my)

                # Centered stencils from averaging upwind operators
                D1x_centered = 0.5 * (stencils.Dpx[kx] + stencils.Dmx[kx])
                D1y_centered = 0.5 * (stencils.Dpy[ky] + stencils.Dmy[ky])

                u2_yλx[idx(i, j, mx, my)] += stencils.Dpx[kx] * material.λ[idx(wrapped_x, j, mx, my)] * stencils.Dmy[ky] * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u2_xμy[idx(i, j, mx, my)] += D1y_centered * material.μ[idx(i, wrapped_y, mx, my)] * D1x_centered * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u1_xλy[idx(i, j, mx, my)] += stencils.Dpy[ky] * material.λ[idx(i, wrapped_y, mx, my)] * stencils.Dmx[kx] * u1[idx(wrapped_x, wrapped_y, mx, my)]
                u1_yμx[idx(i, j, mx, my)] += D1x_centered * material.μ[idx(wrapped_x, j, mx, my)] * D1y_centered * u1[idx(wrapped_x, wrapped_y, mx, my)]
            end
        end
    else
        for kx in 1:Ω_radius
            for ky in 1:Ω_radius
                # Periodic wrapping
                wrapped_x = mod1(i - interior_offset - 1 + kx, mx)
                wrapped_y = mod1(j - interior_offset - 1 + ky, my)
                u2_yλx[idx(i, j, mx, my)] += stencils.D1x[kx] * material.λ[idx(wrapped_x, j, mx, my)] * stencils.D1y[ky] * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u2_xμy[idx(i, j, mx, my)] += stencils.D1y[ky] * material.μ[idx(i, wrapped_y, mx, my)] * stencils.D1x[kx] * u2[idx(wrapped_x, wrapped_y, mx, my)]
                u1_xλy[idx(i, j, mx, my)] += stencils.D1y[ky] * material.λ[idx(i, wrapped_y, mx, my)] * stencils.D1x[kx] * u1[idx(wrapped_x, wrapped_y, mx, my)]
                u1_yμx[idx(i, j, mx, my)] += stencils.D1x[kx] * material.μ[idx(wrapped_x, j, mx, my)] * stencils.D1y[ky] * u1[idx(wrapped_x, wrapped_y, mx, my)]
            end
        end
    end
end