include("../types.jl")
include("../stencils.jl")
include("boundary_derivatives.jl")
include("derivative_helpers.jl")
include("derivative_kernels.jl")
include("inner_derivatives.jl")
include("periodic_derivatives.jl")

using CUDA

# If non-periodic, divide domain into 9 blocks (in 2D) and compute derivatives point-wise,
# where we check which block point belong to
function compute_derivatives(u1, u2, SCHEMA, PERIODIC, grid::Grid, stencils::Stencils, material::Material,
    derivatives::Derivatives, device_data::Union{DeviceData,Nothing}, D2_stencils::Union{D2Stencils, Nothing}; plus_or_minus=nothing, skip_narrow=false)
    mx = grid.mx
    my = grid.my

    # Pre-allocate rows for material properties IF material parameters are time dependent
    # Also needed for the periodic case where we still run with this implementation
    if PERIODIC
        rowλ_x = Vector{Float64}(undef, mx)
        rowμ_x = Vector{Float64}(undef, mx)
        rowλ_y = Vector{Float64}(undef, my)
        rowμ_y = Vector{Float64}(undef, my)
    end

    # Process boundary regions (serial - small fraction of work)
    BW_D1 = stencils.BW_D1
    BW_D2 = stencils.BW_D2

    for field in fieldnames(Derivatives)
        fill!(getfield(derivatives, field), 0.0)
    end

    u1_x = copy(derivatives.u1_x)
    u1_y = copy(derivatives.u1_y)
    u2_x = copy(derivatives.u2_x)
    u2_y = copy(derivatives.u2_y)

    u1_xλx = copy(derivatives.u1_xλx)
    u1_xμx = copy(derivatives.u1_xμx)
    u1_yμy = copy(derivatives.u1_yμy)
    u2_yλy = copy(derivatives.u2_yλy)
    u2_yμy = copy(derivatives.u2_yμy)
    u2_xμx = copy(derivatives.u2_xμx)

    u2_yλx = copy(derivatives.u2_yλx)
    u2_xμy = copy(derivatives.u2_xμy)
    u1_xλy = copy(derivatives.u1_xλy)
    u1_yμx = copy(derivatives.u1_yμx)

    if Params.USE_CUDA && CUDA.functional()
        device_data.u1_x .= 0.0f0
        device_data.u1_y .= 0.0f0
        device_data.u2_x .= 0.0f0
        device_data.u2_y .= 0.0f0

        device_data.u1_xλx .= 0.0f0
        device_data.u1_xμx .= 0.0f0
        device_data.u1_yμy .= 0.0f0
        device_data.u2_yλy .= 0.0f0
        device_data.u2_yμy .= 0.0f0
        device_data.u2_xμx .= 0.0f0

        device_data.u2_yλx .= 0.0f0
        device_data.u2_xμy .= 0.0f0
        device_data.u1_xλy .= 0.0f0
        device_data.u1_yμx .= 0.0f0

        u1_x = device_data.u1_x
        u1_y = device_data.u1_y
        u2_x = device_data.u2_x
        u2_y = device_data.u2_y

        u1_xλx = device_data.u1_xλx
        u1_xμx = device_data.u1_xμx
        u1_yμy = device_data.u1_yμy
        u2_yλy = device_data.u2_yλy
        u2_yμy = device_data.u2_yμy
        u2_xμx = device_data.u2_xμx

        u2_yλx = device_data.u2_yλx
        u2_xμy = device_data.u2_xμy
        u1_xλy = device_data.u1_xλy
        u1_yμx = device_data.u1_yμx

        λ = device_data.λ
        μ = device_data.μ
        D1x = device_data.D1x
        D1y = device_data.D1y
        Dpx = device_data.Dpx
        Dpy = device_data.Dpy
        Dmx = device_data.Dmx
        Dmy = device_data.Dmy

        # Launch kernel for interior points with 2D grid
        threads_per_block = (16, 16)  # 256 threads total per block
        blocks_x = ceil(Int, mx / threads_per_block[1])
        blocks_y = ceil(Int, my / threads_per_block[2])

        hx = Params.PRECISION(grid.hx)
        hy = Params.PRECISION(grid.hy)
        if PERIODIC
            # These first three have identical arguments
            if Params.SCHEMA == :narrow
                # CUDA.@device_code_warntype
                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_narrow_kernel!(
                    u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                    u1, u2, mx, my, hx, hy, stencils.Ω_radius, D1x, D1y, λ, μ
                )
            elseif Params.SCHEMA == :wide
                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_wide_kernel!(
                    u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                    u1, u2, mx, my, hx, hy, stencils.Ω_radius, D1x, D1y, λ, μ
                )
            elseif Params.SCHEMA == :mixed
                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_mixed_kernel!(
                    u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                    u1, u2, mx, my, hx, hy, stencils.Ω_radius, D1x, D1y, λ, μ
                )
            elseif Params.SCHEMA == :upwind
                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_upwind_kernel!(
                    u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                    u1, u2, mx, my, hx, hy, stencils.Ω_radius,
                    Dpx, Dpy, Dmx, Dmy, λ, μ
                )
            elseif Params.SCHEMA == :mixed_upwind_μ
                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_mixed_mu_kernel!(
                    u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                    u1, u2, mx, my, hx, hy, stencils.Ω_radius, stencils.narrow_Ω_radius,
                    Dpx, Dpy, Dmx, Dmy, λ, μ
                )
            elseif Params.SCHEMA == :mixed_upwind_λ
                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_mixed_lambda_kernel!(
                    u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                    u1, u2, mx, my, hx, hy, stencils.Ω_radius, stencils.narrow_Ω_radius,
                    Dpx, Dpy, Dmx, Dmy, λ, μ
                )
            else
                error("CUDA kernel not implemented for schema $(Params.SCHEMA)")
            end
        else
            # Flatten all 2D boundary stencils since CUDA requires flat arrays
            if SCHEMA in (:narrow, :wide, :mixed)
                D1x_∂Ω_W = device_data.D1x_∂Ω_W
                D1x_∂Ω_E = device_data.D1x_∂Ω_E
                D1y_∂Ω_N = device_data.D1y_∂Ω_N
                D1y_∂Ω_S = device_data.D1y_∂Ω_S

                D2x_λ_∂Ω_W_vec = device_data.D2x_λ_∂Ω_W_vec
                D2x_μ_∂Ω_W_vec = device_data.D2x_μ_∂Ω_W_vec
                D2x_λ_∂Ω_E_vec = device_data.D2x_λ_∂Ω_E_vec
                D2x_μ_∂Ω_E_vec = device_data.D2x_μ_∂Ω_E_vec
                D2x_λ_vec = device_data.D2x_λ_vec
                D2x_μ_vec = device_data.D2x_μ_vec

                D2y_λ_∂Ω_N_vec = device_data.D2y_λ_∂Ω_N_vec
                D2y_μ_∂Ω_N_vec = device_data.D2y_μ_∂Ω_N_vec
                D2y_λ_∂Ω_S_vec = device_data.D2y_λ_∂Ω_S_vec
                D2y_μ_∂Ω_S_vec = device_data.D2y_μ_∂Ω_S_vec
                D2y_μ_vec = device_data.D2y_μ_vec
                D2y_λ_vec = device_data.D2y_λ_vec
            elseif SCHEMA in (:upwind, :mixed_upwind_μ, :mixed_upwind_λ)
                if plus_or_minus == :plus
                    D1x_∂Ω_W = device_data.Dpx_∂Ω_W
                    D1x_∂Ω_E = device_data.Dpx_∂Ω_E
                    D1y_∂Ω_N = device_data.Dpy_∂Ω_N
                    D1y_∂Ω_S = device_data.Dpy_∂Ω_S
                    D1x = device_data.Dpx
                    D1y = device_data.Dpy
                elseif plus_or_minus == :minus
                    D1x_∂Ω_W = device_data.Dmx_∂Ω_W
                    D1x_∂Ω_E = device_data.Dmx_∂Ω_E
                    D1y_∂Ω_N = device_data.Dmy_∂Ω_N
                    D1y_∂Ω_S = device_data.Dmy_∂Ω_S
                    D1x = device_data.Dmx
                    D1y = device_data.Dmy
                else
                    error("plus_or_minus must be :plus or :minus")
                end
            else
                error("Unknown SCHEMA: $SCHEMA")
            end
            # println("Non-periodic CUDA not yet implemented")

            @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_D1_boundary_kernel!(
                u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                u1, u2, mx, my,
                λ, μ, stencils.∂Ω_radius, stencils.Ω_radius, BW_D1,
                D1x_∂Ω_W, D1x_∂Ω_E,
                D1y_∂Ω_N, D1y_∂Ω_S,
                D1x, D1y
            )
            # CUDA.synchronize()

            @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_D1_inner_kernel!(
                u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                u1, u2, mx, my, stencils.Ω_radius, λ, μ, D1x, D1y, BW_D1
            )
            # CUDA.synchronize()


            if SCHEMA in (:narrow, :mixed, :mixed_upwind_μ, :mixed_upwind_λ)
                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_D2_narrow_boundary_kernel!(
                    u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u1, u2, mx, my, λ, μ, stencils.∂Ω_radius, stencils.Ω_radius, BW_D2, hx, hy,
                    D2x_λ_∂Ω_W_vec, D2x_μ_∂Ω_W_vec, D2x_λ_∂Ω_E_vec, D2x_μ_∂Ω_E_vec,
                    D2y_λ_∂Ω_N_vec, D2y_μ_∂Ω_N_vec, D2y_λ_∂Ω_S_vec, D2y_μ_∂Ω_S_vec,
                    D2x_λ_vec, D2x_μ_vec, D2y_λ_vec, D2y_μ_vec
                )
                # CUDA.synchronize()

                # Probably need something like this
                if Params.SCHEMA in (:mixed_upwind_μ, :mixed_upwind_λ)
                    Ω_radius = stencils.narrow_Ω_radius
                else
                    Ω_radius = stencils.Ω_radius
                end

                @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_D2_narrow_inner_kernel!(
                    u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                    u1, u2, mx, my, λ, μ, Ω_radius, BW_D2, hx, hy,
                    D2x_λ_vec, D2x_μ_vec, D2y_λ_vec, D2y_μ_vec
                )
                # -- Fused kernel alternative, for narrow --
                # @cuda threads = threads_per_block blocks = (blocks_x, blocks_y) apply_narrow_inner_kernel!(
                #     u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx,
                #     u1, u2, mx, my, Ω_radius, λ, μ, D1x, D1y, BW_D1,
                #     u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
                #     BW_D2, hx, hy, D2x_λ_vec, D2x_μ_vec, D2y_λ_vec, D2y_μ_vec
                # )
                # CUDA.synchronize()
            end
        end
        # error("Non-periodic CUDA not yet fully implemented")
        # # CUDA.synchronize()
    else # NOT GPU, USING CPU
        if PERIODIC
            for i in 1:mx
                for j in 1:my
                    apply_D_periodic!(SCHEMA, u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy,
                        u1_xλy, u1_yμx, u1, u2, i, j, mx, my, grid, stencils, material, rowλ_x, rowμ_x, rowλ_y, rowμ_y)
                end
            end
        else # NOT PERIODIC
            if SCHEMA in (:narrow, :wide, :mixed, :mixed_upwind_μ, :mixed_upwind_λ)
                D1x_∂Ω_W = stencils.D1x_∂Ω_W
                D1x_∂Ω_E = stencils.D1x_∂Ω_E
                D1y_∂Ω_N = stencils.D1y_∂Ω_N
                D1y_∂Ω_S = stencils.D1y_∂Ω_S
                D1x = stencils.D1x
                D1y = stencils.D1y
            elseif SCHEMA == :upwind
                if plus_or_minus == :plus
                    D1x_∂Ω_W = stencils.Dpx_∂Ω_W
                    D1x_∂Ω_E = stencils.Dpx_∂Ω_E
                    D1y_∂Ω_N = stencils.Dpy_∂Ω_N
                    D1y_∂Ω_S = stencils.Dpy_∂Ω_S
                    D1x = stencils.Dpx
                    D1y = stencils.Dpy
                elseif plus_or_minus == :minus
                    D1x_∂Ω_W = stencils.Dmx_∂Ω_W
                    D1x_∂Ω_E = stencils.Dmx_∂Ω_E
                    D1y_∂Ω_N = stencils.Dmy_∂Ω_N
                    D1y_∂Ω_S = stencils.Dmy_∂Ω_S
                    D1x = stencils.Dmx
                    D1y = stencils.Dmy
                else
                    error("plus_or_minus must be :plus or :minus")
                end
            # elseif SCHEMA in (:mixed_upwind_μ, :mixed_upwind_λ)
            #     D1x_∂Ω_W = stencils.Dpx_∂Ω_W
            #     D1x_∂Ω_E = stencils.Dpx_∂Ω_E
            #     D1y_∂Ω_N = stencils.Dpy_∂Ω_N
            #     D1y_∂Ω_S = stencils.Dpy_∂Ω_S
            #     D1x = stencils.D1x
            #     D1y = stencils.D1y
            else
                error("Unknown SCHEMA: $SCHEMA")
            end

            # D1 boundaries: West and East boundaries (full height)
            for i in 1:BW_D1
                for j in 1:my
                    apply_D1_boundary!(u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx, u1, u2, i, j, mx, my, stencils, material, plus_or_minus,
                        D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S, D1x, D1y)
                end
            end
            for i in mx-BW_D1+1:mx
                for j in 1:my
                    apply_D1_boundary!(u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx, u1, u2, i, j, mx, my, stencils, material, plus_or_minus,
                        D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S, D1x, D1y)
                end
            end
            # North and South boundaries (only interior width, corners already done)
            for j in 1:BW_D1
                for i in BW_D1+1:mx-BW_D1
                    apply_D1_boundary!(u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx, u1, u2, i, j, mx, my, stencils, material, plus_or_minus,
                        D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S, D1x, D1y)
                end
            end
            for j in my-BW_D1+1:my
                for i in BW_D1+1:mx-BW_D1
                    apply_D1_boundary!(u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx, u1, u2, i, j, mx, my, stencils, material, plus_or_minus,
                        D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S, D1x, D1y)
                end
            end

            # Interior points - PARALLEL (bulk of computation, vectorizable)
            Threads.@threads for i in BW_D1+1:mx-BW_D1
                for j in BW_D1+1:my-BW_D1
                    apply_D1_inner!(u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx, u1, u2, i, j, mx, my, stencils, material, D1x, D1y)
                end
            end

            Ω_radius = SCHEMA in (:mixed_upwind_μ, :mixed_upwind_λ) ? stencils.narrow_Ω_radius : stencils.Ω_radius

            if SCHEMA in (:narrow, :mixed, :mixed_upwind_μ, :mixed_upwind_λ) && !skip_narrow
                # D2 boundaries (serial)
                for i in 1:BW_D2
                    for j in 1:my
                        apply_D2_narrow_boundary!(u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u1, u2, i, j, mx, my, Ω_radius, stencils, D2_stencils)
                    end
                end
                for i in mx-BW_D2+1:mx
                    for j in 1:my
                        apply_D2_narrow_boundary!(u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u1, u2, i, j, mx, my, Ω_radius, stencils, D2_stencils)
                    end
                end
                for j in 1:BW_D2
                    for i in BW_D2+1:mx-BW_D2
                        apply_D2_narrow_boundary!(u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u1, u2, i, j, mx, my, Ω_radius, stencils, D2_stencils)
                    end
                end
                for j in my-BW_D2+1:my
                    for i in BW_D2+1:mx-BW_D2
                        apply_D2_narrow_boundary!(u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u1, u2, i, j, mx, my, Ω_radius, stencils, D2_stencils)
                    end
                end

                Threads.@threads for i in BW_D2+1:mx-BW_D2
                    # Each thread needs its own buffers to avoid race conditions
                    # Pre-allocate rows for material properties IF material parameters are time dependent
                    # rowλ_x_local = Vector{Float64}(undef, mx)
                    # rowμ_x_local = Vector{Float64}(undef, mx)
                    # rowλ_y_local = Vector{Float64}(undef, my)
                    # rowμ_y_local = Vector{Float64}(undef, my)

                    for j in BW_D2+1:my-BW_D2
                        apply_D2_narrow_inner!(u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u1, u2, i, j, mx, my, Ω_radius, D2_stencils)
                    end
                end
            end
        end
    end

    return u1_x, u1_y, u2_x, u2_y, u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx, u2_yλx, u2_xμy, u1_xλy, u1_yμx
end
