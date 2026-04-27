function RK4(v1, v2, u1, u2, t, dt, grid::Grid, stencils::Stencils, 
    material::Material, bcs::BoundaryConditions, derivatives::Derivatives, 
    device_data::Union{DeviceData, Nothing}, D2_stencils::Union{D2Stencils, Nothing}, cache::Union{RK4_Cache, Nothing})
    
    if cache !== nothing
        k1v1, k1v2, k1u1, k1u2 = cache.k1v1, cache.k1v2, cache.k1u1, cache.k1u2
        k2v1, k2v2, k2u1, k2u2 = cache.k2v1, cache.k2v2, cache.k2u1, cache.k2u2
        k3v1, k3v2, k3u1, k3u2 = cache.k3v1, cache.k3v2, cache.k3u1, cache.k3u2
        k4v1, k4v2, k4u1, k4u2 = cache.k4v1, cache.k4v2, cache.k4u1, cache.k4u2
        rhs_cache = cache.rhs_cache
    else
        rhs_cache = nothing
    end
    
    k1v1_res, k1v2_res, k1u1_res, k1u2_res = RHS(v1, v2, u1, u2, t, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rhs_cache)
    if cache !== nothing
        k1v1 .= k1v1_res; k1v2 .= k1v2_res; k1u1 .= k1u1_res; k1u2 .= k1u2_res
    else
        k1v1, k1v2, k1u1, k1u2 = k1v1_res, k1v2_res, k1u1_res, k1u2_res
    end
    
    k2v1_res, k2v2_res, k2u1_res, k2u2_res = RHS(v1 .+ 0.5f0*dt .* k1v1, v2 .+ 0.5f0*dt .* k1v2, u1 .+ 0.5f0*dt .* k1u1, u2 .+ 0.5f0*dt .* k1u2, t .+ 0.5f0*dt, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rhs_cache)
    if cache !== nothing
        k2v1 .= k2v1_res; k2v2 .= k2v2_res; k2u1 .= k2u1_res; k2u2 .= k2u2_res
    else
        k2v1, k2v2, k2u1, k2u2 = k2v1_res, k2v2_res, k2u1_res, k2u2_res
    end
    
    k3v1_res, k3v2_res, k3u1_res, k3u2_res = RHS(v1 .+ 0.5f0*dt .* k2v1, v2 .+ 0.5f0*dt .* k2v2, u1 .+ 0.5f0*dt .* k2u1, u2 .+ 0.5f0*dt .* k2u2, t .+ 0.5f0*dt, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rhs_cache)
    if cache !== nothing
        k3v1 .= k3v1_res; k3v2 .= k3v2_res; k3u1 .= k3u1_res; k3u2 .= k3u2_res
    else
        k3v1, k3v2, k3u1, k3u2 = k3v1_res, k3v2_res, k3u1_res, k3u2_res
    end
    
    k4v1_res, k4v2_res, k4u1_res, k4u2_res = RHS(v1 .+ dt .* k3v1, v2 .+ dt .* k3v2, u1 .+ dt .* k3u1, u2 .+ dt .* k3u2, t .+ dt, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rhs_cache)
    if cache !== nothing
        k4v1 .= k4v1_res; k4v2 .= k4v2_res; k4u1 .= k4u1_res; k4u2 .= k4u2_res
    else
        k4v1, k4v2, k4u1, k4u2 = k4v1_res, k4v2_res, k4u1_res, k4u2_res
    end
    
    v1_new = v1 .+ dt/6.0f0 .* (k1v1 .+ 2.0f0 .* k2v1 .+ 2.0f0 .* k3v1 .+ k4v1)
    v2_new = v2 .+ dt/6.0f0 .* (k1v2 .+ 2.0f0 .* k2v2 .+ 2.0f0 .* k3v2 .+ k4v2)
    u1_new = u1 .+ dt/6.0f0 .* (k1u1 .+ 2.0f0 .* k2u1 .+ 2.0f0 .* k3u1 .+ k4u1)
    u2_new = u2 .+ dt/6.0f0 .* (k1u2 .+ 2.0f0 .* k2u2 .+ 2.0f0 .* k3u2 .+ k4u2)

    return v1_new, v2_new, u1_new, u2_new
end