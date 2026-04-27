module MMS_2D

include("material_functions.jl")
using .MaterialFunctions
using Symbolics

@variables x y t ω

function make_MMS(λ_tag::String, μ_tag::String, D_tag::String)
    Dt = Differential(t)
    Dx = Differential(x)
    Dy = Differential(y)

    # Pick material symbolic expressions from MaterialFunctions
    λ_expr = MaterialFunctions.λ_expr_dict[λ_tag]
    μ_expr = MaterialFunctions.μ_expr_dict[μ_tag]
    D_expr = MaterialFunctions.D_expr_dict[D_tag]
    b_expr = 1 / D_expr  # Buoyancy

    # --- Define manufactured fields ---
    u1_expr = sin(π*x)^2 * sin(π*y)^2 * cos(ω*t)
    u2_expr = sin(π*x)^2 * sin(π*y)^2 * cos(ω*t)
    u1_t_expr = Dt(u1_expr)
    u2_t_expr = Dt(u2_expr)

    sigma_11 = λ_expr*(Dx(u1_expr) + Dy(u2_expr)) + 2*μ_expr*Dx(u1_expr)
    sigma_22 = λ_expr*(Dx(u1_expr) + Dy(u2_expr)) + 2*μ_expr*Dy(u2_expr)
    sigma_12 = μ_expr*(Dy(u1_expr) + Dx(u2_expr))

    u1_tt_expr = b_expr*(Dx(sigma_11) + Dy(sigma_12))
    u2_tt_expr = b_expr*(Dx(sigma_12) + Dy(sigma_22))

    fx_expr = Dt(Dt(u1_expr)) - u1_tt_expr
    fy_expr = Dt(Dt(u2_expr)) - u2_tt_expr

    # Build fast functions
    u1 = build_function(u1_expr, x, y, t, ω; expression=Val{false}) |> eval
    u2 = build_function(u2_expr, x, y, t, ω; expression=Val{false}) |> eval
    u1_t = build_function(expand_derivatives(u1_t_expr), x, y, t, ω; expression=Val{false}) |> eval
    u2_t = build_function(expand_derivatives(u2_t_expr), x, y, t, ω; expression=Val{false}) |> eval
    fx = build_function(expand_derivatives(fx_expr), x, y, t, ω; expression=Val{false}) |> eval
    fy = build_function(expand_derivatives(fy_expr), x, y, t, ω; expression=Val{false}) |> eval

    return (
        u1=u1, u2=u2, u1_t=u1_t, u2_t=u2_t,
        fx=fx, fy=fy
    )
end

end
