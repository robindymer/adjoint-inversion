#=
Isotropic elastic wave equation in 2D
Using SBP-SAT with traction-free BC
We're solving (in form of constant material parameters, D=1):
u1_tt = (λ+2μ)u1_xx + μu1_yy + (λ+μ)u2_xy
u2_tt = μu2_xx + (λ+2μ)u2_yy + (λ+μ)u1_xy

In first time second order space order form:
v1_t = (λ+2μ)u1_xx + μu1_yy + (λ+μ)u2_xy
v2_t = μu2_xx + (λ+2μ)u2_yy + (λ+μ)u1_xy
u1_t = v1
u2_t = v2

Solution vector: w = [v1, v2, u1, u2]^T

With variable material parameters:
u1_tt = 1/D * [∂x(λ∂x u1) + 2∂x(μ∂x u1) + ∂y(μ∂y u1) + ∂x(λ∂y u2) + ∂y(μ∂x u2)]
u2_tt = 1/D * [∂x(μ∂x u2) + ∂y(λ∂y u2) + 2∂y(μ∂y u2) + ∂y(λ∂x u1) + ∂x(μ∂y u1)]
=#

config_file = length(ARGS) > 0 ? "configs/params_$(ARGS[1]).jl" : "configs/params_default.jl"
include(config_file)

include("types.jl")
include("derivatives/derivatives.jl")
include("helpers.jl")
include("stencils.jl")
include("diracDisc.jl")
include("boundaryConditions.jl")
include("boundaryConditionKernels.jl")
include("material_functions.jl")
include("MMS_2D.jl")

include("core/indexing.jl")
include("core/initial_conditions.jl")
include("core/rhs.jl")
include("core/simulation.jl")
include("core/time_integration.jl")
include("tests/checkSymmetry.jl")

using .MMS_2D
using .MaterialFunctions
using .Params
using Plots
using WriteVTK
using CUDA
using BenchmarkTools
# using Profile
# using PProf

mms = MMS_2D.make_MMS(Params.λ_tag, Params.μ_tag, Params.D_tag)

function MMS_measure_convergence()
    SBP_ORDER = Params.SBP_ORDER
    x_l = Params.x_l
    x_r = Params.x_r
    y_l = Params.y_l
    y_r = Params.y_r

    ω = Params.ω
    T = Params.T
    CFL = Params.CFL

    m_values = Params.m_values
    L2_errors = []
    for m in m_values
        mx = m
        my = m
        λ, μ, D = get_populated_material(Params.λ_tag, Params.μ_tag, Params.D_tag, mx, my)
        println("Running simulation with mx=$mx, my=$my")
        L2_err = simulate(x_l, x_r, y_l, y_r, mx, my, D, λ, μ, ω, T, CFL, SBP_ORDER)
        push!(L2_errors, L2_err)
        println("L2 error at T=$T: $L2_err")
    end

    println("Convergence results:")
    for i in 2:length(m_values)
        rate = log(L2_errors[i-1] / L2_errors[i]) / log(m_values[i] / m_values[i-1])
        println("From mx=$(m_values[i-1]) to mx=$(m_values[i]): L2 error=$(L2_errors[i]), rate=$rate")
    end
end

function dry_run()
    SBP_ORDER = Params.SBP_ORDER
    x_l = Params.x_l
    x_r = Params.x_r
    y_l = Params.y_l
    y_r = Params.y_r

    mx = Params.mx
    my = Params.my

    if Params.USE_MARMOUSI
        λ, μ, D = get_populated_marmousi_material(Params.DATA_DIR, mx, my, Params.DOWNSAMPLE_FACTOR)
        
        # Plot lambda to verify
        λ_2d = reshape(λ, mx, my)
        heatmap(λ_2d', 
            title="Marmousi Lambda (λ) - $mx x $my grid",
            xlabel="X", ylabel="Y",
            c=:viridis, yflip=true)
        savefig("$(Params.ASSETS_PATH)/marmousi_lambda_verification.png")
        println("Saved verification plot to $(Params.ASSETS_PATH)/marmousi_lambda_verification.png")
    else
        λ, μ, D = get_populated_material(Params.λ_tag, Params.μ_tag, Params.D_tag, mx, my)
    end
    ω = Params.ω
    T = Params.T
    CFL = Params.CFL

    simulate(x_l, x_r, y_l, y_r, mx, my, D, λ, μ, ω, T, CFL, SBP_ORDER)
end

if Params.USE_CUDA
    try
        CUDA.reclaim()  # Reset CUDA context if corrupted
    catch
        # First time, nothing to reclaim
    end

    if !CUDA.functional()
        @warn "CUDA not functional! Falling back to CPU or restart Julia"
    end
end

if Params.CONVERGENCE_STUDY
    @time MMS_measure_convergence()
    # println("\n=== Second run (without compilation) ===")
    # @time MMS_measure_convergence()
else
    # dry_run()
    # println("\n=== Second run (without compilation) ===")
    # @btime dry_run()
    @time dry_run()
end
