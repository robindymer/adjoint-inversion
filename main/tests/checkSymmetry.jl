using LinearAlgebra
using Printf
# using MAT
# NOTE: Make sure const D_tag = "constant" is set in params for this to work

# Check if x^T HD y = y^T HD x holds (for verification)
# NOTE: Not sure about this approach, therefore we have something
# slightly more rigorous down below
########
# function checkSymmetry(grid::Grid, stencils::Stencils, material::Material, derivatives::Derivatives, device_data::Union{DeviceData, Nothing})
#     mx = grid.mx
#     my = grid.my
#     iD = 1.0 ./ material.D

#     x1 = randn(mx*my)
#     x2 = randn(mx*my)
#     y1 = randn(mx*my)
#     y2 = randn(mx*my)
#     lhs = 0.0
#     rhs = 0.0

#     u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy, u1_xλy, u1_yμx = compute_derivatives(x1, x2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data)
#     Dx_1 = iD .* (u1_xλx .+ 2.0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
#     Dx_2 = iD .* (u2_xμx .+ u2_yλy .+ 2.0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
#     u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy, u1_xλy, u1_yμx = compute_derivatives(y1, y2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data)
#     Dy_1 = iD .* (u1_xλx .+ 2.0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
#     Dy_2 = iD .* (u2_xμx .+ u2_yλy .+ 2.0 .* u2_yμy .+ u1_xλy .+ u1_yμx)

#     for j in 1:grid.my
#         for i in 1:grid.mx
#             lhs +=  x1[idx(i, j, mx, my)] * stencils.Hx[i] * stencils.Hy[j] * Dy_1[idx(i, j, mx, my)]
#             lhs +=  x2[idx(i, j, mx, my)] * stencils.Hx[i] * stencils.Hy[j] * Dy_2[idx(i, j, mx, my)]

#             rhs +=  y1[idx(i, j, mx, my)] * stencils.Hx[i] * stencils.Hy[j] * Dx_1[idx(i, j, mx, my)]
#             rhs +=  y2[idx(i, j, mx, my)] * stencils.Hx[i] * stencils.Hy[j] * Dx_2[idx(i, j, mx, my)]
#         end
#     end
#     # Should it instead be lhs + rhs ≈ 0 because of skew-symmetry?
#     println("Symmetry violation (should be close to 0): ", lhs - rhs)
#     println("Skew symmetry violation (should be close to 0): ", lhs + rhs)
# end
########

# Test symmetry, matrix approach
function buildOperatorMatrix(grid::Grid, stencils::Stencils, material::Material, 
                             derivatives::Derivatives, bcs::BoundaryConditions, device_data::Union{DeviceData, Nothing},
                             D2_stencils::Union{D2Stencils, Nothing}, cache::Union{RK4_Cache, Nothing})
    n = grid.mx * grid.my
    D = zeros(2*n, 2*n)  # 2n x 2n for the coupled (u1, u2) system
    iD = 1.0 ./ material.D

    mx = grid.mx
    my = grid.my
    λ = material.λ
    μ = material.μ
    quad_term_x = stencils.quad_term_x
    quad_term_y = stencils.quad_term_y
    
    println("Building operator matrix of size $(2*n) x $(2*n)")
    
    for col in 1:(2*n)
        if col % 100 == 0
            println("  Processing column $col / $(2*n)")
        end
        
        # Create col-th canonical basis vector
        e = zeros(2*n)
        e[col] = 1.0
        
        # Split into u1 and u2 components
        u1 = e[1:n]
        u2 = e[n+1:end]
        
        v1 = zeros(n)
        v2 = zeros(n)
        t = 0.0
        rhs_cache = cache === nothing ? nothing : cache.rhs_cache
        Du1, Du2, _, _ = RHS(v1, v2, u1, u2, t, grid, stencils, material, bcs, derivatives, device_data, D2_stencils, rhs_cache)

        # Store as col-th column
        D[1:n, col] = Du1
        D[n+1:end, col] = Du2
    end

    # file = matopen("operator_matrix.mat", "w")
    # write(file, "M", M)
    # close(file)   
    
    # -- TEST OPERATOR (broken right now) --
    # gaussian(x, y) = exp(-5.0 * (x^2 + y^2))
    # gaussian_vec = zeros(2*n)
    # for j in 1:grid.my
    #     for i in 1:grid.mx
    #         k = idx(i, j, grid.mx, grid.my)
    #         x = grid.x_values[i]
    #         y = grid.y_values[j]
    #         gaussian_vec[k] = gaussian(x, y)
    #         gaussian_vec[n + k] = gaussian(x, y)
    #     end
    # end
    # D_vec = D * gaussian_vec
    # u1 = gaussian_vec[1:n]
    # u2 = gaussian_vec[n+1:end]
    # u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy, u1_xλy, u1_yμx = 
    #     compute_derivatives(u1, u2, Params.SCHEMA, Params.PERIODIC, grid, stencils, material, derivatives, device_data)
    # Du1 = iD .* (u1_xλx .+ 2.0 .* u1_xμx .+ u1_yμy .+ u2_yλx .+ u2_xμy)
    # Du2 = iD .* (u2_xμx .+ u2_yλy .+ 2.0 .* u2_yμy .+ u1_xλy .+ u1_yμx)
    # D_vec_test = vcat(Du1, Du2)
    # diff = maximum(abs.(D_vec - D_vec_test))
    # println("Max difference when applying operator matrix vs direct computation: $diff")

    # -- PLOTTING --
    # heatmap(reshape(D_vec_test[1:n], grid.mx, grid.my)', title="D applied to Gaussian (u1 component)", xlabel="x", ylabel="y")
    # savefig("Du1_test_gaussian.png")
    # heatmap(reshape(D_vec_test[n+1:end], grid.mx, grid.my)', title="D applied to Gaussian (u2 component)", xlabel="x", ylabel="y")
    # savefig("Du2_test_gaussian.png")

    println("Operator matrix built")
    return D
end

function buildQuadratureMatrix(grid::Grid, stencils::Stencils)
    n = grid.mx * grid.my
    
    # Build 2D quadrature weights as diagonal matrix
    H_2D = zeros(n, n)
    for j in 1:grid.my
        for i in 1:grid.mx
            k = idx(i, j, grid.mx, grid.my)
            H_2D[k, k] = stencils.Hx[i] * stencils.Hy[j]
        end
    end
    
    # For coupled system (u1, u2): H = kron(I(2), H_2D)
    # This creates a block diagonal matrix:
    # H = [H_2D   0   ]
    #     [0     H_2D]
    H = zeros(2*n, 2*n)
    H[1:n, 1:n] = H_2D
    H[n+1:end, n+1:end] = H_2D
    
    return H
end

function testSymmetry(grid::Grid, stencils::Stencils, material::Material, 
                     derivatives::Derivatives, bcs::BoundaryConditions, device_data::Union{DeviceData, Nothing}, 
                     D2_stencils::Union{D2Stencils, Nothing}, cache::Union{RK4_Cache, Nothing})
    println("\n=== Testing SBP Symmetry Properties ===")
    
    # Build matrices
    D = buildOperatorMatrix(grid, stencils, material, derivatives, bcs, device_data, D2_stencils, cache)
    H = buildQuadratureMatrix(grid, stencils)
    
    # Check H*D symmetry
    HD = H * D

    # -- CHECK THAT H = H^T >= 0 --
    # println("\nChecking quadrature matrix H:")
    # symmetry_H = maximum(abs.(H - H'))
    # positive_definite = minimum(eigvals(Symmetric(H))) >= -1e-12
    # println("  Max |H - H'|:              $symmetry_H")
    # println("  Positive definite:         $positive_definite")
    
    # Symmetry: HD - (HD)' should be ~0
    symmetry_violation = maximum(abs.(HD - (HD)'))
    relative_symmetry = symmetry_violation / maximum(abs.(HD))

    symmetry_vialation_counter = 0
    symmetry_counter = 0
    for i in 1:size(HD, 1)
        for j in 1:size(HD, 2)
            if abs(HD[i, j] - HD[j, i]) > 1e-10
                # println("Symmetry violation at ($i, $j): HD[i,j] = $(HD[i, j]), HD[j,i] = $(HD[j, i])")
                symmetry_vialation_counter += 1
            else
                symmetry_counter += 1
            end
        end
    end

    println("  Number of symmetry violations > 1e-10: $symmetry_vialation_counter")
    println("  Number of symmetric entries <= 1e-10: $symmetry_counter")
    
    # Skew-symmetry: HD + (HD)' should be large (not skew-symmetric)
    skew_violation = maximum(abs.(HD + (HD)'))
    
    println("\nResults:")
    println("  Max |HD - HD'|:            $symmetry_violation")
    println("  Relative symmetry error:   $relative_symmetry")
    println("  Max |HD + HD'|:            $skew_violation (should be large)")
    println("  Max |HD|:                  $(maximum(abs.(HD)))")

    # file = matopen("M_diff.mat", "w")
    # write(file, "M_diff", HD - (HD)')
    # close(file)  
    
    # # For small matrices, print them
    # if size(HD, 1) <= 32
    #     println("\nOperator HD (first 32x32 block):")
    #     for i in 1:min(32, size(HD, 1))
    #         for j in 1:min(32, size(HD, 2))
    #             print(@sprintf("%7.3f ", HD[i, j]))
    #         end
    #         println()
    #     end
    # end
    
    # Check eigenvalues for energy stability
    println("\nComputing eigenvalues (this may take a while)...")
    evals = eigvals(Symmetric(HD))
    max_eval = maximum(evals)
    min_eval = minimum(evals)
    
    println("  Max eigenvalue:            $max_eval  (should be ≤ 0)")
    println("  Min eigenvalue:            $min_eval")
    println("  Energy stable:             $(max_eval <= 1e-10)")
    
    # Print summary
    println("\n=== Summary ===")
    if relative_symmetry < 1e-10
        println("✓ Operator is symmetric (excellent)")
    elseif relative_symmetry < 1e-6
        println("✓ Operator is nearly symmetric (acceptable)")
    else
        println("✗ WARNING: Large symmetry violation!")
    end
    
    if max_eval <= 1e-10
        println("✓ Operator is energy stable")
    else
        println("✗ WARNING: Operator has positive eigenvalues!")
    end
    
    return D, H, HD, evals
end