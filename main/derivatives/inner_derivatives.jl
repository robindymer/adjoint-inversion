@inline function apply_D1_inner!(u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                                u1, u2, i, j, mx, my, stencils, material, D1x, D1y)
    λ = material.λ
    μ = material.μ

    Ω_radius = stencils.Ω_radius # k max for interior
    # interior_offset = stencils.SBP_ORDER == 2 ? 1 : 2 # Or floor(Ω_radius/2)
    interior_offset = Ω_radius ÷ 2
    # Interior points
    for k in 1:Ω_radius
        u1_x[idx(i, j, mx, my)] += D1x[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u2_x[idx(i, j, mx, my)] += D1x[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]

        u1_y[idx(i, j, mx, my)] += D1y[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_y[idx(i, j, mx, my)] += D1y[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
    end
    for kx in 1:Ω_radius
        for ky in 1:Ω_radius
            u2_yλx[idx(i, j, mx, my)] += D1x[kx] * λ[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y[ky] * u2[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
            u2_xμy[idx(i, j, mx, my)] += D1y[ky] * μ[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x[kx] * u2[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
            u1_xλy[idx(i, j, mx, my)] += D1y[ky] * λ[idx(i, j - interior_offset - 1 + ky, mx, my)] * D1x[kx] * u1[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
            u1_yμx[idx(i, j, mx, my)] += D1x[kx] * μ[idx(i - interior_offset - 1 + kx, j, mx, my)] * D1y[ky] * u1[idx(i - interior_offset - 1 + kx, j - interior_offset - 1 + ky, mx, my)]
        end
    end
end

@inline function apply_D2_narrow_inner!(u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u1, u2, i, j, 
    mx, my, Ω_radius, D2_stencils::D2Stencils)
    # NOTE: Saving this since if material parameters change in time, we need this!
    # extract_row_x!(rowλ_x, material.λ, j, mx, my)
    # extract_row_x!(rowμ_x, material.μ, j, mx, my)
    # extract_col_y!(rowλ_y, material.λ, i, mx, my)
    # extract_col_y!(rowμ_y, material.μ, i, mx, my)

    # D2x_λ = D2_inner(rowλ_x, i, grid.hx, stencils.SBP_ORDER)
    # D2x_μ = D2_inner(rowμ_x, i, grid.hx, stencils.SBP_ORDER)
    # D2y_μ = D2_inner(rowμ_y, j, grid.hy, stencils.SBP_ORDER)
    # D2y_λ = D2_inner(rowλ_y, j, grid.hy, stencils.SBP_ORDER)
    D2x_λ = D2_stencils.D2x_λ_vec[i, j]
    D2x_μ = D2_stencils.D2x_μ_vec[i, j]
    D2y_μ = D2_stencils.D2y_μ_vec[i, j]
    D2y_λ = D2_stencils.D2y_λ_vec[i, j]

    interior_offset = Ω_radius ÷ 2

    for k in 1:Ω_radius
        u1_xλx[idx(i, j, mx, my)] += D2x_λ[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u1_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u1[idx(i - interior_offset - 1 + k, j, mx, my)]
        u2_xμx[idx(i, j, mx, my)] += D2x_μ[k] * u2[idx(i - interior_offset - 1 + k, j, mx, my)]

        u1_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u1[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_yμy[idx(i, j, mx, my)] += D2y_μ[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
        u2_yλy[idx(i, j, mx, my)] += D2y_λ[k] * u2[idx(i, j - interior_offset - 1 + k, mx, my)]
    end
end