###########################
#### TRACTION-FREE BC ####
##########################

function apply_traction_free_west_kernel!(v1_t, v2_t, u1, u2, u1_y, u2_y, 
                                          d1x_l, d1x_l_wide, d_radius, λ, μ, iD, 
                                          quad_term_x, tau_W1, tau_W2, mx, my, schema_is_mixed)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    
    if j <= my
        # Compute derivatives
        ux1 = 0.0f0
        ux2 = 0.0f0
        ux1_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(k, j, mx, my)
            ux1 += d1x_l[k] * u1[idx_k]
            ux2 += d1x_l[k] * u2[idx_k]
            ux1_wide += d1x_l_wide[k] * u1[idx_k]
        end
        
        idx_bdry = idx(1, j, mx, my)
        u1y = u1_y[idx_bdry]
        u2y = u2_y[idx_bdry]
        
        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        
        # Compute stress components
        sigma_11 = schema_is_mixed ? λ_idx*ux1_wide + 2.0f0*μ_idx*ux1 + λ_idx*u2y : (λ_idx + 2.0f0*μ_idx)*ux1 + λ_idx*u2y
        sigma_12 = μ_idx*(u1y + ux2)
        
        # Update velocity time derivatives
        iD_val = iD[idx_bdry]
        v1_t[idx_bdry] += tau_W1 * quad_term_x * iD_val * (-sigma_11)
        v2_t[idx_bdry] += tau_W2 * quad_term_x * iD_val * (-sigma_12)
    end
    
    return nothing
end

function apply_traction_free_east_kernel!(v1_t, v2_t, u1, u2, u1_y, u2_y, 
                                         d1x_r, d1x_r_wide, d_radius, λ, μ, iD, 
                                         quad_term_x, tau_E1, tau_E2, mx, my, schema_is_mixed)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    
    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    if j <= my
        # Compute derivatives
        ux1 = 0.0f0
        ux2 = 0.0f0
        ux1_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(mx - d_radius + k, j, mx, my)
            ux1 += d1x_r[k] * u1[idx_k]
            ux2 += d1x_r[k] * u2[idx_k]
            ux1_wide += d1x_r_wide[k] * u1[idx_k]
        end
        
        idx_bdry = idx(mx, j, mx, my)
        u1y = u1_y[idx_bdry]
        u2y = u2_y[idx_bdry]
        
        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        
        # Compute stress components
        sigma_11 = schema_is_mixed ? λ_idx*ux1_wide + 2.0f0*μ_idx*ux1 + λ_idx*u2y : (λ_idx + 2.0f0*μ_idx)*ux1 + λ_idx*u2y
        sigma_12 = μ_idx*(u1y + ux2)
        
        # Update velocity time derivatives
        iD_val = iD[idx_bdry]
        v1_t[idx_bdry] += tau_E1 * quad_term_x * iD_val * (sigma_11)
        v2_t[idx_bdry] += tau_E2 * quad_term_x * iD_val * (sigma_12)
    end
    
    return nothing
end

function apply_traction_free_south_kernel!(v1_t, v2_t, u1, u2, u1_x, u2_x, 
                                          d1y_l, d1y_l_wide, d_radius, λ, μ, iD, 
                                          quad_term_y, tau_S1, tau_S2, mx, my, schema_is_mixed)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    
    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    if i <= mx
        # Compute derivatives
        uy1 = 0.0f0
        uy2 = 0.0f0
        uy2_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(i, k, mx, my)
            uy1 += d1y_l[k] * u1[idx_k]
            uy2 += d1y_l[k] * u2[idx_k]
            uy2_wide += d1y_l_wide[k] * u2[idx_k]
        end
        
        idx_bdry = idx(i, 1, mx, my)
        u1x = u1_x[idx_bdry]
        u2x = u2_x[idx_bdry]
        
        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        
        # Compute stress components
        sigma_12 = μ_idx*(uy1 + u2x)
        sigma_22 = schema_is_mixed ? λ_idx*u1x + λ_idx*uy2_wide + 2.0f0*μ_idx*uy2 : λ_idx*u1x + (λ_idx + 2.0f0*μ_idx)*uy2
        
        # Update velocity time derivatives
        iD_val = iD[idx_bdry]
        v1_t[idx_bdry] += tau_S1 * quad_term_y * iD_val * (-sigma_12)
        v2_t[idx_bdry] += tau_S2 * quad_term_y * iD_val * (-sigma_22)
    end
    
    return nothing
end

function apply_traction_free_north_kernel!(v1_t, v2_t, u1, u2, u1_x, u2_x, 
                                          d1y_r, d1y_r_wide, d_radius, λ, μ, iD, 
                                          quad_term_y, tau_N1, tau_N2, mx, my, schema_is_mixed)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inline idx(i, j, mx, my) = (j - 1) * mx + i
    
    if i <= mx
        # Compute derivatives
        uy1 = 0.0f0
        uy2 = 0.0f0
        uy2_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(i, my - d_radius + k, mx, my)
            uy1 += d1y_r[k] * u1[idx_k]
            uy2 += d1y_r[k] * u2[idx_k]
            uy2_wide += d1y_r_wide[k] * u2[idx_k]
        end
        
        idx_bdry = idx(i, my, mx, my)
        u1x = u1_x[idx_bdry]
        u2x = u2_x[idx_bdry]
        
        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        
        # Compute stress components
        sigma_12 = μ_idx*(uy1 + u2x)
        sigma_22 = schema_is_mixed ? λ_idx*u1x + λ_idx*uy2_wide + 2.0f0*μ_idx*uy2 : λ_idx*u1x + (λ_idx + 2.0f0*μ_idx)*uy2
        
        # Update velocity time derivatives
        iD_val = iD[idx_bdry]
        v1_t[idx_bdry] += tau_N1 * quad_term_y * iD_val * (sigma_12)
        v2_t[idx_bdry] += tau_N2 * quad_term_y * iD_val * (sigma_22)
    end
    
    return nothing
end


###########################
#### NON-REFLECTING BC ####
###########################

function apply_non_reflecting_west_kernel!(v1_t, v2_t, v1, v2, u1, u2, u1_y, u2_y, 
                                         d1x_l, d1x_l_wide, d_radius, λ, μ, D, iD, 
                                         quad_term_x, tau_W1, tau_W2, mx, my, schema_is_mixed)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    if j <= my
        # Compute stress components at boundary using displacement derivatives
        ux1 = 0.0f0
        ux2 = 0.0f0
        ux1_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(k, j, mx, my)
            ux1 += d1x_l[k] * u1[idx_k]
            ux2 += d1x_l[k] * u2[idx_k]
            ux1_wide += d1x_l_wide[k] * u1[idx_k]
        end
        
        idx_bdry = idx(1, j, mx, my)
        u1y = u1_y[idx_bdry]
        u2y = u2_y[idx_bdry]

        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        rho_idx = D[idx_bdry]
        c_p = sqrt((λ_idx + 2.0f0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        # sigma_11 = (λ+2μ)*u1_x + λ*u2_y
        # sigma_12 = μ*(u1_y + u2_x)
        # Absorbing BC: v1 = sigma_11/(ρ*c_p), v2 = sigma_12/(ρ*c_s)
        sigma_11 = schema_is_mixed ? λ_idx*ux1_wide + 2.0f0*μ_idx*ux1 + λ_idx*u2y : (λ_idx + 2.0f0*μ_idx)*ux1 + λ_idx*u2y
        sigma_12 = μ_idx*(u1y + ux2)

        v1_t[idx_bdry] += tau_W1 * quad_term_x * iD[idx_bdry] * (v1[idx_bdry]*rho_idx*c_p - sigma_11)
        v2_t[idx_bdry] += tau_W2 * quad_term_x * iD[idx_bdry] * (v2[idx_bdry]*rho_idx*c_s - sigma_12)
    end
    
    return nothing
end

function apply_non_reflecting_east_kernel!(v1_t, v2_t, v1, v2, u1, u2, u1_y, u2_y, 
                                        d1x_r, d1x_r_wide, d_radius, λ, μ, D, iD, 
                                        quad_term_x, tau_E1, tau_E2, mx, my, schema_is_mixed)
    j = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    if j <= my
        ux1 = 0.0f0
        ux2 = 0.0f0
        ux1_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(mx - d_radius + k, j, mx, my)
            ux1 += d1x_r[k] * u1[idx_k]
            ux2 += d1x_r[k] * u2[idx_k]
            ux1_wide += d1x_r_wide[k] * u1[idx_k]
        end

        idx_bdry = idx(mx, j, mx, my)
        u1y = u1_y[idx_bdry]
        u2y = u2_y[idx_bdry]

        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        rho_idx = D[idx_bdry]
        c_p = sqrt((λ_idx + 2.0f0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        sigma_11 = schema_is_mixed ? λ_idx*ux1_wide + 2.0f0*μ_idx*ux1 + λ_idx*u2y : (λ_idx + 2.0f0*μ_idx)*ux1 + λ_idx*u2y
        sigma_12 = μ_idx*(u1y + ux2)

        v1_t[idx_bdry] += tau_E1 * quad_term_x * iD[idx_bdry] * (v1[idx_bdry]*rho_idx*c_p + sigma_11)
        v2_t[idx_bdry] += tau_E2 * quad_term_x * iD[idx_bdry] * (v2[idx_bdry]*rho_idx*c_s + sigma_12)
    end

    return nothing
end

function apply_non_reflecting_south_kernel!(v1_t, v2_t, v1, v2, u1, u2, u1_x, u2_x, 
                                         d1y_l, d1y_l_wide, d_radius, λ, μ, D, iD, 
                                         quad_term_y, tau_S1, tau_S2, mx, my, schema_is_mixed)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    if i <= mx
        # Compute derivatives
        uy1 = 0.0f0
        uy2 = 0.0f0
        uy2_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(i, k, mx, my)
            uy1 += d1y_l[k] * u1[idx_k]
            uy2 += d1y_l[k] * u2[idx_k]
            uy2_wide += d1y_l_wide[k] * u2[idx_k]
        end
        
        idx_bdry = idx(i, 1, mx, my)
        u1x = u1_x[idx_bdry]
        u2x = u2_x[idx_bdry]

        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        rho_idx = D[idx_bdry]
        c_p = sqrt((λ_idx + 2.0f0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        # sigma_12 = μ*(u1_y + u2_x)
        # sigma_22 = λ*u1_x + (λ+2μ)*u2_y
        sigma_12 = μ_idx*(uy1 + u2x)
        sigma_22 = schema_is_mixed ? λ_idx*u1x + λ_idx*uy2_wide + 2.0f0*μ_idx*uy2 : λ_idx*u1x + (λ_idx + 2.0f0*μ_idx)*uy2

        v1_t[idx_bdry] += tau_S1 * quad_term_y * iD[idx_bdry] * (v1[idx_bdry]*rho_idx*c_s - sigma_12)
        v2_t[idx_bdry] += tau_S2 * quad_term_y * iD[idx_bdry] * (v2[idx_bdry]*rho_idx*c_p - sigma_22)
    end

    return nothing
end

function apply_non_reflecting_north_kernel!(v1_t, v2_t, v1, v2, u1, u2, u1_x, u2_x, 
                                         d1y_r, d1y_r_wide, d_radius, λ, μ, D, iD, 
                                         quad_term_y, tau_N1, tau_N2, mx, my, schema_is_mixed)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inline idx(i, j, mx, my) = (j - 1) * mx + i

    if i <= mx
        # Compute derivatives
        uy1 = 0.0f0
        uy2 = 0.0f0
        uy2_wide = 0.0f0
        for k in 1:d_radius
            idx_k = idx(i, my - d_radius + k, mx, my)
            uy1 += d1y_r[k] * u1[idx_k]
            uy2 += d1y_r[k] * u2[idx_k]
            uy2_wide += d1y_r_wide[k] * u2[idx_k]
        end
        
        idx_bdry = idx(i, my, mx, my)
        u1x = u1_x[idx_bdry]
        u2x = u2_x[idx_bdry]

        λ_idx = λ[idx_bdry]
        μ_idx = μ[idx_bdry]
        rho_idx = D[idx_bdry]
        c_p = sqrt((λ_idx + 2.0f0 * μ_idx) / rho_idx)
        c_s = sqrt(μ_idx / rho_idx)

        # sigma_12 = μ*(u1_y + u2_x)
        # sigma_22 = λ*u1_x + (λ+2μ)*u2_y
        sigma_12 = μ_idx*(uy1 + u2x)
        sigma_22 = schema_is_mixed ? λ_idx*u1x + λ_idx*uy2_wide + 2.0f0*μ_idx*uy2 : λ_idx*u1x + (λ_idx + 2.0f0*μ_idx)*uy2

        v1_t[idx_bdry] += tau_N1 * quad_term_y * iD[idx_bdry] * (v1[idx_bdry]*rho_idx*c_s + sigma_12)
        v2_t[idx_bdry] += tau_N2 * quad_term_y * iD[idx_bdry] * (v2[idx_bdry]*rho_idx*c_p + sigma_22)
    end

    return nothing
end
