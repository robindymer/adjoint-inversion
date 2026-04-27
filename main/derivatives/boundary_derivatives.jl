@inline function apply_D1_boundary!(u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx, u1, u2, i, j, mx, my, stencils, material, plus_or_minus,
                            D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S, D1x, D1y)
    λ = material.λ
    μ = material.μ

    ∂Ω_radius = stencils.∂Ω_radius # k max for boundary
    Ω_radius = stencils.Ω_radius # k max for interior
    # interior_offset = stencils.SBP_ORDER == 2 ? 1 : 2 # Or floor(Ω_radius/2)
    interior_offset = Ω_radius ÷ 2
    interior_radius_y = my - 2 * stencils.BW_D1
    interior_radius_x = mx - 2 * stencils.BW_D1
    if i <= stencils.BW_D1
        if j <= stencils.BW_D1
            # NW corner!
            # Dx, Dxx, Dy and Dyy at boundary
            for k in 1:∂Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x_∂Ω_W[i][k] * u1[idx(k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_∂Ω_W[i][k] * u2[idx(k, j, mx, my)]

                u1_y[idx(i, j, mx, my)] += D1y_∂Ω_N[j][k] * u1[idx(i, k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_∂Ω_N[j][k] * u2[idx(i, k, mx, my)]
            end
            for kx in 1:∂Ω_radius
                for ky in 1:∂Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x_∂Ω_W[i][kx] * λ[idx(kx, j, mx, my)] * D1y_∂Ω_N[j][ky] * u2[idx(kx, ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y_∂Ω_N[j][ky] * μ[idx(i, ky, mx, my)] * D1x_∂Ω_W[i][kx] * u2[idx(kx, ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y_∂Ω_N[j][ky] * λ[idx(i, ky, mx, my)] * D1x_∂Ω_W[i][kx] * u1[idx(kx, ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x_∂Ω_W[i][kx] * μ[idx(kx, j, mx, my)] * D1y_∂Ω_N[j][ky] * u1[idx(kx, ky, mx, my)]
                end
            end
        elseif j >= my - stencils.BW_D1 + 1
            # SW corner!
            for k in 1:∂Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x_∂Ω_W[i][k] * u1[idx(k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_∂Ω_W[i][k] * u2[idx(k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_y[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][k] * u1[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
            end
            for kx in 1:∂Ω_radius
                for ky in 1:∂Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x_∂Ω_W[i][kx] * λ[idx(kx, j, mx, my)] * D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * u2[idx(kx, my - ∂Ω_radius + ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * μ[idx(i, my - ∂Ω_radius + ky, mx, my)] * D1x_∂Ω_W[i][kx] * u2[idx(kx, my - ∂Ω_radius + ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * λ[idx(i, my - ∂Ω_radius + ky, mx, my)] * D1x_∂Ω_W[i][kx] * u1[idx(kx, my - ∂Ω_radius + ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x_∂Ω_W[i][kx] * μ[idx(kx, j, mx, my)] * D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * u1[idx(kx, my - ∂Ω_radius + ky, mx, my)]
                end
            end
        else
            # W edge!
            for k in 1:∂Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x_∂Ω_W[i][k] * u1[idx(k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_∂Ω_W[i][k] * u2[idx(k, j, mx, my)]
            end
            # interior y-derivatives
            for k in 1:Ω_radius
                u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                # println("i = $i, j = $j, k = $k, idx(i, j - interior_offset - 1 + k, mx, my) = $(idx(i, j - interior_offset - 1 + k, mx, my))")
            end
            for kx in 1:∂Ω_radius
                for ky in 1:Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x_∂Ω_W[i][kx] * λ[idx(kx, j, mx, my)] * D1y[ky] * u2[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y[ky] * μ[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_∂Ω_W[i][kx] * u2[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y[ky] * λ[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_∂Ω_W[i][kx] * u1[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x_∂Ω_W[i][kx] * μ[idx(kx, j, mx, my)] * D1y[ky] * u1[idx(kx, j - interior_offset - 1 + ky, mx, my)]
                end
            end
        end
    elseif i >= mx - stencils.BW_D1 + 1
        if j <= stencils.BW_D1
            # NE corner
            for k in 1:∂Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][k] * u2[idx(mx - ∂Ω_radius + k, j, mx, my)]

                u1_y[idx(i, j, mx, my)] += D1y_∂Ω_N[j][k] * u1[idx(i, k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_∂Ω_N[j][k] * u2[idx(i, k, mx, my)]
            end
            for kx in 1:∂Ω_radius
                for ky in 1:∂Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * λ[idx(mx - ∂Ω_radius + kx, j, mx, my)] * D1y_∂Ω_N[j][ky] * u2[idx(mx - ∂Ω_radius + kx, ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y_∂Ω_N[j][ky] * μ[idx(i, ky, mx, my)] * D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * u2[idx(mx - ∂Ω_radius + kx, ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y_∂Ω_N[j][ky] * λ[idx(i, ky, mx, my)] * D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * u1[idx(mx - ∂Ω_radius + kx, ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * μ[idx(mx - ∂Ω_radius + kx, j, mx, my)] * D1y_∂Ω_N[j][ky] * u1[idx(mx - ∂Ω_radius + kx, ky, mx, my)]
                end
            end
        elseif j >= my - stencils.BW_D1 + 1
            # SE corner
            for k in 1:∂Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][k] * u2[idx(mx - ∂Ω_radius + k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_y[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][k] * u1[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
            end
            for kx in 1:∂Ω_radius
                for ky in 1:∂Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * λ[idx(mx - ∂Ω_radius + kx, j, mx, my)] * D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * u2[idx(mx - ∂Ω_radius + kx, my - ∂Ω_radius + ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * μ[idx(i, my - ∂Ω_radius + ky, mx, my)] * D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * u2[idx(mx - ∂Ω_radius + kx, my - ∂Ω_radius + ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * λ[idx(i, my - ∂Ω_radius + ky, mx, my)] * D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * u1[idx(mx - ∂Ω_radius + kx, my - ∂Ω_radius + ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * μ[idx(mx - ∂Ω_radius + kx, j, mx, my)] * D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * u1[idx(mx - ∂Ω_radius + kx, my - ∂Ω_radius + ky, mx, my)]
                end
            end
        else
            # E edge
            for k in 1:∂Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][k] * u2[idx(mx - ∂Ω_radius + k, j, mx, my)]
            end
            # interior y-derivatives
            for k in 1:Ω_radius
                u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
            end
            for kx in 1:∂Ω_radius
                for ky in 1:Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * λ[idx(mx - ∂Ω_radius + kx, j, mx, my)] * D1y[ky] * u2[idx(mx - ∂Ω_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y[ky] * μ[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * u2[idx(mx - ∂Ω_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y[ky] * λ[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * u1[idx(mx - ∂Ω_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x_∂Ω_E[i - (stencils.BW_D1 + interior_radius_x)][kx] * μ[idx(mx - ∂Ω_radius + kx, j, mx, my)] * D1y[ky] * u1[idx(mx - ∂Ω_radius + kx, j - interior_offset - 1 + ky, mx, my)]
                end
            end
        end
    else
        if j <= stencils.BW_D1
            # N edge
            for k in 1:∂Ω_radius
                u1_y[idx(i, j, mx, my)] += D1y_∂Ω_N[j][k] * u1[idx(i, k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_∂Ω_N[j][k] * u2[idx(i, k, mx, my)]
            end
            # interior x-derivatives
            for k in 1:Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
            # Mixed derivatives
            for kx in 1:Ω_radius
                for ky in 1:∂Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x[kx] * λ[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_∂Ω_N[j][ky] * u2[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y_∂Ω_N[j][ky] * μ[idx(i, ky, mx, my)] * D1x[kx] * u2[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y_∂Ω_N[j][ky] * λ[idx(i, ky, mx, my)] * D1x[kx] * u1[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x[kx] * μ[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_∂Ω_N[j][ky] * u1[idx(i - interior_offset - 1 + kx, ky, mx, my)]
                end
            end
        elseif j >= my - stencils.BW_D1 + 1
            # S edge
            for k in 1:∂Ω_radius
                # Since stencils already are mirrored, we go from down to up
                u1_y[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][k] * u1[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_y[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
            end
            # interior x-derivatives
            for k in 1:Ω_radius
                u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
            # Mixed derivatives
            for kx in 1:Ω_radius
                for ky in 1:∂Ω_radius
                    u2_yλx[idx(i, j, mx, my)] += D1x[kx] * λ[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * u2[idx(i - interior_offset - 1 + kx, my - ∂Ω_radius + ky, mx, my)]
                    u2_xμy[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * μ[idx(i, my - ∂Ω_radius + ky, mx, my)] * D1x[kx] * u2[idx(i - interior_offset - 1 + kx, my - ∂Ω_radius + ky, mx, my)]
                    u1_xλy[idx(i, j, mx, my)] += D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * λ[idx(i, my - ∂Ω_radius + ky, mx, my)] * D1x[kx] * u1[idx(i - interior_offset - 1 + kx, my - ∂Ω_radius + ky, mx, my)]
                    u1_yμx[idx(i, j, mx, my)] += D1x[kx] * μ[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y_∂Ω_S[j - (stencils.BW_D1 + interior_radius_y)][ky] * u1[idx(i - interior_offset - 1 + kx, my - ∂Ω_radius + ky, mx, my)]
                end
            end
        end
    end
end

@inline function apply_D2_narrow_boundary!(u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u1, u2, i, j, mx, my, 
    Ω_radius, stencils::Stencils, D2_stencils::Union{D2Stencils, Nothing})
    # NOTE: Saving this since if material parameters change in time, we need this!
    # extract_row_x!(rowλ_x, material.λ, j, mx, my)
    # extract_row_x!(rowμ_x, material.μ, j, mx, my)
    # extract_col_y!(rowλ_y, material.λ, i, mx, my)
    # extract_col_y!(rowμ_y, material.μ, i, mx, my)

    # D2x_λ_∂Ω_W = D2_∂Ω_l(rowλ_x, grid.hx, stencils.SBP_ORDER)
    # D2x_λ_∂Ω_E = D2_∂Ω_r(rowλ_x, grid.hx, stencils.SBP_ORDER)
    # D2x_μ_∂Ω_W = D2_∂Ω_l(rowμ_x, grid.hx, stencils.SBP_ORDER)
    # D2x_μ_∂Ω_E = D2_∂Ω_r(rowμ_x, grid.hx, stencils.SBP_ORDER)

    # D2y_λ_∂Ω_N = D2_∂Ω_l(rowλ_y, grid.hy, stencils.SBP_ORDER)
    # D2y_λ_∂Ω_S = D2_∂Ω_r(rowλ_y, grid.hy, stencils.SBP_ORDER)
    # D2y_μ_∂Ω_N = D2_∂Ω_l(rowμ_y, grid.hy, stencils.SBP_ORDER)
    # D2y_μ_∂Ω_S = D2_∂Ω_r(rowμ_y, grid.hy, stencils.SBP_ORDER)

    # if (i > stencils.BW_D2) && (i < mx - stencils.BW_D2 + 1)
    #     D2x_λ = D2_inner(rowλ_x, i, grid.hx, stencils.SBP_ORDER)
    #     D2x_μ = D2_inner(rowμ_x, i, grid.hx, stencils.SBP_ORDER)
    # end
    # if (j > stencils.BW_D2) && (j < my - stencils.BW_D2 + 1)
    #     D2y_μ = D2_inner(rowμ_y, j, grid.hy, stencils.SBP_ORDER)
    #     D2y_λ = D2_inner(rowλ_y, j, grid.hy, stencils.SBP_ORDER)
    # end

    D2x_λ_∂Ω_W = D2_stencils.D2x_λ_∂Ω_W_vec[j]
    D2x_λ_∂Ω_E = D2_stencils.D2x_λ_∂Ω_E_vec[j]
    D2x_μ_∂Ω_W = D2_stencils.D2x_μ_∂Ω_W_vec[j]
    D2x_μ_∂Ω_E = D2_stencils.D2x_μ_∂Ω_E_vec[j]
    
    D2y_λ_∂Ω_N = D2_stencils.D2y_λ_∂Ω_N_vec[i]
    D2y_λ_∂Ω_S = D2_stencils.D2y_λ_∂Ω_S_vec[i]
    D2y_μ_∂Ω_N = D2_stencils.D2y_μ_∂Ω_N_vec[i]
    D2y_μ_∂Ω_S = D2_stencils.D2y_μ_∂Ω_S_vec[i]

    D2x_λ = D2_stencils.D2x_λ_vec[i, j]
    D2x_μ = D2_stencils.D2x_μ_vec[i, j]

    D2y_μ = D2_stencils.D2y_μ_vec[i, j]
    D2y_λ = D2_stencils.D2y_λ_vec[i, j]

    ∂Ω_radius = stencils.∂Ω_radius # k max for boundary
    # interior_offset = stencils.SBP_ORDER == 2 ? 1 : 2 # Or floor(Ω_radius/2)
    interior_offset = Ω_radius ÷ 2
    interior_radius_y = my - 2 * stencils.BW_D2
    interior_radius_x = mx - 2 * stencils.BW_D2

    if i <= stencils.BW_D2
        if j <= stencils.BW_D2
            # NW corner!
            # Dx, Dxx, Dy and Dyy at boundary
            for k in 1:∂Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ_∂Ω_W[i, k] * u1[idx(k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_W[i, k] * u1[idx(k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_W[i, k] * u2[idx(k, j, mx, my)]

                u1_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_N[j, k] * u1[idx(i, k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_N[j, k] * u2[idx(i, k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ_∂Ω_N[j, k] * u2[idx(i, k, mx, my)]
            end
        elseif j >= my - stencils.BW_D2 + 1
            # SW corner!
            for k in 1:∂Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ_∂Ω_W[i, k] * u1[idx(k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_W[i, k] * u1[idx(k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_W[i, k] * u2[idx(k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u1[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
            end
        else
            # W edge!
            for k in 1:∂Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ_∂Ω_W[i, k] * u1[idx(k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_W[i, k] * u1[idx(k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_W[i, k] * u2[idx(k, j, mx, my)]
            end
            # interior y-derivatives
            for k in 1:Ω_radius
                u1_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                # println("i = $i, j = $j, k = $k, idx(i, j - interior_offset - 1 + k, mx, my) = $(idx(i, j - interior_offset - 1 + k, mx, my))")
            end
        end
    elseif i >= mx - stencils.BW_D2 + 1
        if j <= stencils.BW_D2
            # NE corner
            for k in 1:∂Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u2[idx(mx - ∂Ω_radius + k, j, mx, my)]

                u1_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_N[j, k] * u1[idx(i, k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_N[j, k] * u2[idx(i, k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ_∂Ω_N[j, k] * u2[idx(i, k, mx, my)]
            end
        elseif j >= my - stencils.BW_D2 + 1
            # SE corner
            for k in 1:∂Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u2[idx(mx - ∂Ω_radius + k, j, mx, my)]

                # Since stencils already are mirrored, we go from down to up
                u1_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u1[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
            end
        else
            # E edge
            for k in 1:∂Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u1[idx(mx - ∂Ω_radius + k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ_∂Ω_E[i - (stencils.BW_D2 + interior_radius_x), k] * u2[idx(mx - ∂Ω_radius + k, j, mx, my)]
            end
            # interior y-derivatives
            for k in 1:Ω_radius
                u1_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
            end
        end
    else
        if j <= stencils.BW_D2
            # N edge
            for k in 1:∂Ω_radius
                u1_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_N[j, k] * u1[idx(i, k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_N[j, k] * u2[idx(i, k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ_∂Ω_N[j, k] * u2[idx(i, k, mx, my)]
            end
            # interior x-derivatives
            for k in 1:Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
        elseif j >= my - stencils.BW_D2 + 1
            # S edge
            for k in 1:∂Ω_radius
                # Since stencils already are mirrored, we go from down to up
                u1_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u1[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_yμy[idx(i, j, mx, my)] += D2y_μ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
                u2_yλy[idx(i, j, mx, my)] += D2y_λ_∂Ω_S[j - (stencils.BW_D2 + interior_radius_y), k] * u2[idx(i, my - ∂Ω_radius + k, mx, my)]
            end
            # interior x-derivatives
            for k in 1:Ω_radius
                u1_xλx[idx(i, j, mx, my)] += D2x_λ[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u1_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
                u2_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]
            end
        end
    end
end