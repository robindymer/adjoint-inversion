###########################
#### TRACTION-FREE BC ####
##########################

## Not staggered ##
###################

function apply_traction_free_west!(schema, v1_t, v2_t, u1, u2, u1_y, u2_y, mx, my, 
                                   d1x_l, d1x_l_wide, d_radius, λ, μ, iD, quad_term_x, tau_W1, tau_W2)
    # Left boundary (x = x_l)
    for j in 1:my        
        # -(λ+2μ)u1_x - λu2_y = 0
        # -μ(u2_x + u1_y) = 0
        ux1 = 0.0
        ux2 = 0.0
        ux1_wide = 0.0
        for k in 1:d_radius
            ux1 += d1x_l[k]*u1[idx(k, j, mx, my)]
            ux2 += d1x_l[k]*u2[idx(k, j, mx, my)]
            ux1_wide += d1x_l_wide[k]*u1[idx(k, j, mx, my)] # NOTE: Cleaner, but perfomance issue?
        end
            
        u1y = u1_y[idx(1, j, mx, my)]
        u2y = u2_y[idx(1, j, mx, my)]

        λ_idx = λ[idx(1, j, mx, my)]
        μ_idx = μ[idx(1, j, mx, my)]

        if schema == :mixed
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        elseif schema == :mixed_upwind_λ
            # For this case, wide = upwind and rest is narrow (centered)
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        else
            sigma_11 = (λ_idx + 2.0*μ_idx)*ux1 + λ_idx*u2y
        end
        sigma_12 = μ_idx*(u1y + ux2)

        v1_t[idx(1, j, mx, my)] += tau_W1 * quad_term_x * iD[idx(1, j, mx, my)] * ( -sigma_11 )
        v2_t[idx(1, j, mx, my)] += tau_W2 * quad_term_x * iD[idx(1, j, mx, my)] * ( -sigma_12 )
    end
end

function apply_traction_free_east!(schema, v1_t, v2_t, u1, u2, u1_y, u2_y, mx, my, 
                                   d1x_r, d1x_r_wide, d_radius, λ, μ, iD, quad_term_x, tau_E1, tau_E2)
    # Right boundary (x = x_r)
    for j in 1:my
        # (λ+2μ)u1_x + λu2_y = 0
        # μ(u2_x + u1_y) = 0
        ux1 = 0.0
        ux2 = 0.0
        ux1_wide = 0.0
        for k in 1:d_radius
            ux1 += d1x_r[k]*u1[idx(mx - d_radius + k, j, mx, my)]
            ux2 += d1x_r[k]*u2[idx(mx - d_radius + k, j, mx, my)]
            ux1_wide += d1x_r_wide[k]*u1[idx(mx - d_radius + k, j, mx, my)]
        end

        u1y = u1_y[idx(mx, j, mx, my)]
        u2y = u2_y[idx(mx, j, mx, my)]

        λ_idx = λ[idx(mx, j, mx, my)]
        μ_idx = μ[idx(mx, j, mx, my)]

        if schema == :mixed
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        elseif schema == :mixed_upwind_λ
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        else
            sigma_11 = (λ_idx + 2.0*μ_idx)*ux1 + λ_idx*u2y
        end
        sigma_12 = μ_idx*(u1y + ux2)

        v1_t[idx(mx, j, mx, my)] += tau_E1 * quad_term_x * iD[idx(mx, j, mx, my)] * ( sigma_11 )
        v2_t[idx(mx, j, mx, my)] += tau_E2 * quad_term_x * iD[idx(mx, j, mx, my)] * ( sigma_12 )
    end
end

function apply_traction_free_south!(schema, v1_t, v2_t, u1, u2, u1_x, u2_x, mx, my, 
                                    d1y_l, d1y_l_wide, d_radius, λ, μ, iD, quad_term_y, tau_S1, tau_S2)
    # Bottom boundary (y = y_l)
    for i in 1:mx
        # -μ(u2_x + u1_y) = 0
        # -λu1_x - (λ+2μ)u2_y = 0
        uy1 = 0.0
        uy2 = 0.0
        uy2_wide = 0.0
        for k in 1:d_radius
            uy1 += d1y_l[k]*u1[idx(i, k, mx, my)]
            uy2 += d1y_l[k]*u2[idx(i, k, mx, my)]
            uy2_wide += d1y_l_wide[k]*u2[idx(i, k, mx, my)]
        end

        u1x = u1_x[idx(i, 1, mx, my)]
        u2x = u2_x[idx(i, 1, mx, my)]

        λ_idx = λ[idx(i, 1, mx, my)]
        μ_idx = μ[idx(i, 1, mx, my)]

        sigma_12 = μ_idx*(uy1 + u2x)
        if schema == :mixed
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        elseif schema == :mixed_upwind_λ
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        else
            sigma_22 = λ_idx*u1x + (λ_idx+2μ_idx)*uy2
        end

        v1_t[idx(i, 1, mx, my)] += tau_S1 * quad_term_y * iD[idx(i, 1, mx, my)] * ( -sigma_12 )
        v2_t[idx(i, 1, mx, my)] += tau_S2 * quad_term_y * iD[idx(i, 1, mx, my)] * ( -sigma_22 )
    end
end

function apply_traction_free_north!(schema, v1_t, v2_t, u1, u2, u1_x, u2_x, mx, my, 
                                    d1y_r, d1y_r_wide, d_radius, λ, μ, iD, quad_term_y, tau_N1, tau_N2)
    # Top boundary (y = y_r)
    for i in 1:mx
        # μ(u2_x + u1_y) = 0
        # λu1_x + (λ+2μ)u2_y = 0
        uy1 = 0.0
        uy2 = 0.0
        uy2_wide = 0.0
        for k in 1:d_radius
            uy1 += d1y_r[k]*u1[idx(i, my - d_radius + k, mx, my)]
            uy2 += d1y_r[k]*u2[idx(i, my - d_radius + k, mx, my)]
            uy2_wide += d1y_r_wide[k]*u2[idx(i, my - d_radius + k, mx, my)]
        end

        u1x = u1_x[idx(i, my, mx, my)]
        u2x = u2_x[idx(i, my, mx, my)]

        λ_idx = λ[idx(i, my, mx, my)]
        μ_idx = μ[idx(i, my, mx, my)]

        sigma_12 = μ_idx*(uy1 + u2x)
        if schema == :mixed
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        elseif schema == :mixed_upwind_λ
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        else
            sigma_22 = λ_idx*u1x + (λ_idx+2μ_idx)*uy2
        end

        v1_t[idx(i, my, mx, my)] += tau_N1 * quad_term_y * iD[idx(i, my, mx, my)] * ( sigma_12 )
        v2_t[idx(i, my, mx, my)] += tau_N2 * quad_term_y * iD[idx(i, my, mx, my)] * ( sigma_22 )
    end
end

###########################
#### NON-REFLECTING BC ####
###########################

function apply_non_reflecting_west!(schema, v1_t, v2_t, v1, v2, u1, u2, u1_y, u2_y, mx, my, 
                                    d1x_l, d1x_l_wide, d_radius, λ, μ, D, iD, quad_term_x, tau_W1, tau_W2)
    # Left boundary (x = x_l) - waves propagating in -x direction
    for j in 1:my        
        # Compute stress components at boundary using displacement derivatives
        ux1 = 0.0
        ux2 = 0.0
        ux1_wide = 0.0
        for k in 1:d_radius
            ux1 += d1x_l[k]*u1[idx(k, j, mx, my)]
            ux2 += d1x_l[k]*u2[idx(k, j, mx, my)]
            ux1_wide += d1x_l_wide[k]*u1[idx(k, j, mx, my)]
        end
        
        u1y = u1_y[idx(1, j, mx, my)]
        u2y = u2_y[idx(1, j, mx, my)]

        λ_idx = λ[idx(1, j, mx, my)]
        μ_idx = μ[idx(1, j, mx, my)]
        rho_idx = D[idx(1, j, mx, my)]
        c_p = sqrt((λ_idx + 2.0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        # sigma_11 = (λ+2μ)*u1_x + λ*u2_y
        # sigma_12 = μ*(u1_y + u2_x)
        # Absorbing BC: v1 = sigma_11/(ρ*c_p), v2 = sigma_12/(ρ*c_s)
        if schema == :mixed
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        elseif schema == :mixed_upwind_λ
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        else
            sigma_11 = (λ_idx + 2.0*μ_idx)*ux1 + λ_idx*u2y
        end
        sigma_12 = μ_idx*(u1y + ux2)
        
        v1_t[idx(1, j, mx, my)] += tau_W1 * quad_term_x * iD[idx(1, j, mx, my)] * (v1[idx(1, j, mx, my)]*rho_idx*c_p - sigma_11)
        v2_t[idx(1, j, mx, my)] += tau_W2 * quad_term_x * iD[idx(1, j, mx, my)] * (v2[idx(1, j, mx, my)]*rho_idx*c_s - sigma_12)
    end
end

function apply_non_reflecting_east!(schema, v1_t, v2_t, v1, v2, u1, u2, u1_y, u2_y, mx, my, 
                                   d1x_r, d1x_r_wide, d_radius, λ, μ, D, iD, quad_term_x, tau_E1, tau_E2)
    # Right boundary (x = x_r) - waves propagating in +x direction  
    for j in 1:my
        ux1 = 0.0
        ux2 = 0.0
        ux1_wide = 0.0
        for k in 1:d_radius
            ux1 += d1x_r[k]*u1[idx(mx - d_radius + k, j, mx, my)]
            ux2 += d1x_r[k]*u2[idx(mx - d_radius + k, j, mx, my)]
            ux1_wide += d1x_r_wide[k]*u1[idx(mx - d_radius + k, j, mx, my)]
        end

        u1y = u1_y[idx(mx, j, mx, my)]
        u2y = u2_y[idx(mx, j, mx, my)]

        λ_idx = λ[idx(mx, j, mx, my)]
        μ_idx = μ[idx(mx, j, mx, my)]
        rho_idx = D[idx(mx, j, mx, my)]
        c_p = sqrt((λ_idx + 2.0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        if schema == :mixed
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        elseif schema == :mixed_upwind_λ
            sigma_11 = λ_idx*ux1_wide + 2.0*μ_idx*ux1 + λ_idx*u2y
        else
            sigma_11 = (λ_idx + 2.0*μ_idx)*ux1 + λ_idx*u2y
        end
        sigma_12 = μ_idx*(u1y + ux2)
        
        v1_t[idx(mx, j, mx, my)] += tau_E1 * quad_term_x * iD[idx(mx, j, mx, my)] * (v1[idx(mx, j, mx, my)]*rho_idx*c_p + sigma_11)
        v2_t[idx(mx, j, mx, my)] += tau_E2 * quad_term_x * iD[idx(mx, j, mx, my)] * (v2[idx(mx, j, mx, my)]*rho_idx*c_s + sigma_12)
    end
end

function apply_non_reflecting_south!(schema, v1_t, v2_t, v1, v2, u1, u2, u1_x, u2_x, mx, my, 
                                    d1y_l, d1y_l_wide, d_radius, λ, μ, D, iD, quad_term_y, tau_S1, tau_S2)
    # Bottom boundary (y = y_l) - waves propagating in -y direction
    for i in 1:mx
        uy1 = 0.0
        uy2 = 0.0
        uy2_wide = 0.0
        for k in 1:d_radius
            uy1 += d1y_l[k]*u1[idx(i, k, mx, my)]
            uy2 += d1y_l[k]*u2[idx(i, k, mx, my)]
            uy2_wide += d1y_l_wide[k]*u2[idx(i, k, mx, my)]
        end

        u1x = u1_x[idx(i, 1, mx, my)]
        u2x = u2_x[idx(i, 1, mx, my)]

        λ_idx = λ[idx(i, 1, mx, my)]
        μ_idx = μ[idx(i, 1, mx, my)]
        rho_idx = D[idx(i, 1, mx, my)]
        c_p = sqrt((λ_idx + 2.0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        # sigma_22 = λ*u1_x + (λ+2μ)*u2_y
        # sigma_12 = μ*(u1_y + u2_x)
        sigma_12 = μ_idx*(uy1 + u2x)
        if schema == :mixed
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        elseif schema == :mixed_upwind_λ
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        else
            sigma_22 = λ_idx*u1x + (λ_idx+2μ_idx)*uy2
        end
        
        v2_t[idx(i, 1, mx, my)] += tau_S1 * quad_term_y * iD[idx(i, 1, mx, my)] * (v2[idx(i, 1, mx, my)]*rho_idx*c_p - sigma_22)
        v1_t[idx(i, 1, mx, my)] += tau_S2 * quad_term_y * iD[idx(i, 1, mx, my)] * (v1[idx(i, 1, mx, my)]*rho_idx*c_s - sigma_12)
    end
end

function apply_non_reflecting_north!(schema, v1_t, v2_t, v1, v2, u1, u2, u1_x, u2_x, mx, my, 
                                    d1y_r, d1y_r_wide, d_radius, λ, μ, D, iD, quad_term_y, tau_N1, tau_N2)
    # Top boundary (y = y_r) - waves propagating in +y direction
    for i in 1:mx
        uy1 = 0.0
        uy2 = 0.0
        uy2_wide = 0.0
        for k in 1:d_radius
            uy1 += d1y_r[k]*u1[idx(i, my - d_radius + k, mx, my)]
            uy2 += d1y_r[k]*u2[idx(i, my - d_radius + k, mx, my)]
            uy2_wide += d1y_r_wide[k]*u2[idx(i, my - d_radius + k, mx, my)]
        end

        u1x = u1_x[idx(i, my, mx, my)]
        u2x = u2_x[idx(i, my, mx, my)]

        λ_idx = λ[idx(i, my, mx, my)]
        μ_idx = μ[idx(i, my, mx, my)]
        rho_idx = D[idx(i, my, mx, my)]
        c_p = sqrt((λ_idx + 2.0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        sigma_12 = μ_idx*(uy1 + u2x)
        if schema == :mixed
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        elseif schema == :mixed_upwind_λ
            sigma_22 = λ_idx*u1x + λ_idx*uy2_wide+2*μ_idx*uy2
        else
            sigma_22 = λ_idx*u1x + (λ_idx+2μ_idx)*uy2
        end
        
        v2_t[idx(i, my, mx, my)] += tau_N1 * quad_term_y * iD[idx(i, my, mx, my)] * (v2[idx(i, my, mx, my)]*rho_idx*c_p + sigma_22)
        v1_t[idx(i, my, mx, my)] += tau_N2 * quad_term_y * iD[idx(i, my, mx, my)] * (v1[idx(i, my, mx, my)]*rho_idx*c_s + sigma_12)
    end
end