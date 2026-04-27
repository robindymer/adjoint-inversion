module MaterialFunctions

using Symbolics

@variables x y

# --- Symbolic definitions for known materials ---
λ_expr_dict = Dict(
    "trigonometric" => 1.0 + 0.5*sin(pi*x)*sin(pi*y),
    "constant"      => 1.0,
)

μ_expr_dict = Dict(
    "trigonometric" => 0.5 + 0.25*cos(pi*x)*sin(pi*y),
    "heaviside" => 0.5 - 0.45*0.5*(sign(y) + 1.0),
    "heaviside_los" => 0.5*0.5*(sign(-y) + 1.0), # liquid over solid interface
    "heaviside_corner" => 0.05 + 0.45*0.5*(sign(-x) + 1.0)*0.5*(sign(-y) + 1.0), # corner solid
    "constant_05"      => 0.5,
    "constant_005"      => 0.05,
    "constant_0"      => 0.0,
)

D_expr_dict = Dict(
    "constant"      => 1.0,
    "trigonometric" => 1.0 + 1.0 + 0.25*sin(pi*x)*sin(pi*y),
)

# --- Numerical functions (for main solver) ---
function λ_func(tag::String, x::Float64, y::Float64)
    if tag == "trigonometric"
        return 1.0 + 0.5*sin(pi*x)*sin(pi*y)
    elseif tag == "constant"
        return 1.0
    else
        error("Unknown λ tag: $tag")
    end
end

function μ_func(tag::String, x::Float64, y::Float64)
    if tag == "trigonometric"
        return 0.5 + 0.25*cos(pi*x)*sin(pi*y)
    elseif tag == "heaviside"
        return 0.5 - 0.45*0.5*(sign(y) + 1.0)
    elseif tag == "heaviside_los"
        return 0.5*0.5*(sign(-y) + 1.0)
    elseif tag == "heaviside_corner"
        return 0.05 + 0.45*0.5*(sign(-x) + 1.0)*0.5*(sign(-y) + 1.0)
    elseif tag == "constant_05"
        return 0.5
    elseif tag == "constant_005"
        return 0.05
    elseif tag == "constant_0"
        return 0.0
    else
        error("Unknown μ tag: $tag")
    end
end

function D_func(tag::String, x::Float64, y::Float64)
    if tag == "constant"
        return 1.0
    elseif tag == "trigonometric"
        return 1.0 + 1.0 + 0.25*sin(pi*x)*sin(pi*y)
    else
        error("Unknown D tag: $tag")
    end
end

end
