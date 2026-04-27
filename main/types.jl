using CUDA
using StaticArrays

const FloatType = Params.PRECISION

struct Grid
    x_l::Float64
    x_r::Float64
    y_l::Float64
    y_r::Float64
    mx::Int
    my::Int
    hx::Float64
    hy::Float64
    x_values::Vector{Float64}
    y_values::Vector{Float64}

    # # Inner constructor with defaults
    # function Grid(x_l, x_r, y_l, y_r, mx, my, hx, hy, x_values, y_values, 
    #               x_values_dual=nothing, y_values_dual=nothing)
    #     # Determine dual grid dimensions
    #     mx_d = x_values_dual === nothing ? mx : length(x_values_dual)
    #     my_d = y_values_dual === nothing ? my : length(y_values_dual)
    #     new(x_l, x_r, y_l, y_r, mx, my, mx_d, my_d, hx, hy, x_values, y_values, 
    #         x_values_dual, y_values_dual)
    # end
end

# Abstract type for all stencil types
abstract type AbstractStencils end

# Centered stencils (SBP order 2, 4)
struct StencilsCentered <: AbstractStencils
    SBP_ORDER::Int
    ∂Ω_radius::Int
    Ω_radius::Int
    BW_D1::Int
    BW_D2::Int
    D1x::Vector{Float64}
    D1y::Vector{Float64}
    D1x_∂Ω_W::Union{Vector{Float64}, Vector{Vector{Float64}}}
    D1x_∂Ω_E::Union{Vector{Float64}, Vector{Vector{Float64}}}
    D1y_∂Ω_N::Union{Vector{Float64}, Vector{Vector{Float64}}}
    D1y_∂Ω_S::Union{Vector{Float64}, Vector{Vector{Float64}}}
    d1x_l::Vector{Float64}
    d1x_r::Vector{Float64}
    d1y_l::Vector{Float64}
    d1y_r::Vector{Float64}
    d1x_l_wide::Vector{Float64}
    d1x_r_wide::Vector{Float64}
    d1y_l_wide::Vector{Float64}
    d1y_r_wide::Vector{Float64}
    quad_term_x::Float64
    quad_term_y::Float64
    Hx::Vector{Float64}
    Hy::Vector{Float64}
end

# Upwind stencils (SBP order 5, 7, 9)
struct StencilsUpwind <: AbstractStencils
    SBP_ORDER::Int
    ∂Ω_radius::Int
    Ω_radius::Int
    narrow_Ω_radius::Int
    BW_D1::Int
    BW_D2::Int
    Dpx::Vector{Float64}  # D-plus x
    Dpy::Vector{Float64}  # D-plus y
    Dmx::Vector{Float64}  # D-minus x
    Dmy::Vector{Float64}  # D-minus y
    Dpx_∂Ω_W::Vector{Vector{Float64}}
    Dpx_∂Ω_E::Vector{Vector{Float64}}
    Dpy_∂Ω_N::Vector{Vector{Float64}}
    Dpy_∂Ω_S::Vector{Vector{Float64}}
    Dmx_∂Ω_W::Vector{Vector{Float64}}
    Dmx_∂Ω_E::Vector{Vector{Float64}}
    Dmy_∂Ω_N::Vector{Vector{Float64}}
    Dmy_∂Ω_S::Vector{Vector{Float64}}
    d1x_l::Vector{Float64}
    d1x_r::Vector{Float64}
    d1y_l::Vector{Float64}
    d1y_r::Vector{Float64}
    d1x_l_central::Vector{Float64}
    d1x_r_central::Vector{Float64}
    d1y_l_central::Vector{Float64}
    d1y_r_central::Vector{Float64}
    quad_term_x::Float64
    quad_term_y::Float64
    Hx::Vector{Float64}
    Hy::Vector{Float64}
end

# return StencilsMixedUpwind(SBP_ORDER, ∂Ω_radius, Ω_radius, narrow_Ω_radius, BW_D1, BW_D2, D1x, D1y, Dpx, Dpy, Dmx, Dmy,
#                     D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S, Dpx_∂Ω_W, Dpx_∂Ω_E, Dpy_∂Ω_N, Dpy_∂Ω_S,
#                     Dmx_∂Ω_W, Dmx_∂Ω_E, Dmy_∂Ω_N, Dmy_∂Ω_S, d1x_l_upwind, d1x_r_upwind, 
#                     d1y_l_upwind, d1y_r_upwind, d1x_l_wide, d1x_r_wide, d1y_l_wide, 
#                     d1y_r_wide, d1x_l_narrow, d1x_r_narrow, d1y_l_narrow, d1y_r_narrow,
#                     quad_term_x, quad_term_y, Hx, Hy)
struct StencilsMixedUpwind <: AbstractStencils
    SBP_ORDER::Int
    ∂Ω_radius::Int
    Ω_radius::Int
    narrow_Ω_radius::Int
    BW_D1::Int
    BW_D2::Int
    D1x::Vector{Float64}
    D1y::Vector{Float64}
    Dpx::Vector{Float64}  # D-plus x
    Dpy::Vector{Float64}  # D-plus y
    Dmx::Vector{Float64}  # D-minus x
    Dmy::Vector{Float64}  # D-minus y
    D1x_∂Ω_W::Vector{Vector{Float64}}
    D1x_∂Ω_E::Vector{Vector{Float64}}
    D1y_∂Ω_N::Vector{Vector{Float64}}
    D1y_∂Ω_S::Vector{Vector{Float64}}
    Dpx_∂Ω_W::Vector{Vector{Float64}}
    Dpx_∂Ω_E::Vector{Vector{Float64}}
    Dpy_∂Ω_N::Vector{Vector{Float64}}
    Dpy_∂Ω_S::Vector{Vector{Float64}}
    Dmx_∂Ω_W::Vector{Vector{Float64}}
    Dmx_∂Ω_E::Vector{Vector{Float64}}
    Dmy_∂Ω_N::Vector{Vector{Float64}}
    Dmy_∂Ω_S::Vector{Vector{Float64}}
    d1x_l_upwind::Vector{Float64}
    d1x_r_upwind::Vector{Float64}
    d1y_l_upwind::Vector{Float64}
    d1y_r_upwind::Vector{Float64}
    d1x_l_wide::Vector{Float64}
    d1x_r_wide::Vector{Float64}
    d1y_l_wide::Vector{Float64}
    d1y_r_wide::Vector{Float64}
    d1x_l_narrow::Vector{Float64}
    d1x_r_narrow::Vector{Float64}
    d1y_l_narrow::Vector{Float64}
    d1y_r_narrow::Vector{Float64}
    quad_term_x::Float64
    quad_term_y::Float64
    Hx::Vector{Float64}
    Hy::Vector{Float64}
end

# Alias for backward compatibility
if Params.SCHEMA == :upwind
    const Stencils = StencilsUpwind
elseif Params.SCHEMA in (:mixed_upwind_μ, :mixed_upwind_λ)
    const Stencils = StencilsMixedUpwind
else
    const Stencils = StencilsCentered
end

struct Material
    D::Vector{Float64} # Density
    λ::Vector{Float64} # Lame parameter
    μ::Vector{Float64} # Lame parameter
    ω::Float64 # Angular frequency
end

struct BoundaryConditions
    west::Symbol
    east::Symbol
    south::Symbol
    north::Symbol
end

struct Derivatives
    u1_x::Vector{Float64}
    u1_y::Vector{Float64}
    u2_x::Vector{Float64}
    u2_y::Vector{Float64}
    u2_yλx::Vector{Float64}
    u2_xμy::Vector{Float64}
    u1_xλy::Vector{Float64}
    u1_yμx::Vector{Float64}
    u1_xλx::Vector{Float64}
    u1_xμx::Vector{Float64}
    u1_yμy::Vector{Float64}
    u2_yλy::Vector{Float64}
    u2_yμy::Vector{Float64}
    u2_xμx::Vector{Float64}
end

struct D2Stencils{M <: SMatrix, V <: SVector}
    D2x_λ_∂Ω_W_vec::Vector{M}
    D2x_λ_∂Ω_E_vec::Vector{M}
    D2x_μ_∂Ω_W_vec::Vector{M}
    D2x_μ_∂Ω_E_vec::Vector{M}

    D2x_λ_vec::Matrix{V}
    D2x_μ_vec::Matrix{V}

    D2y_λ_∂Ω_N_vec::Vector{M}
    D2y_λ_∂Ω_S_vec::Vector{M}
    D2y_μ_∂Ω_N_vec::Vector{M}
    D2y_μ_∂Ω_S_vec::Vector{M}

    D2y_μ_vec::Matrix{V}
    D2y_λ_vec::Matrix{V}
end

# TODO: Remove Union{..., Nothing} and see if perfomance increases (type stability)
struct DeviceData{M <: SMatrix, V <: SVector}
    u1_x::CuArray{FloatType}
    u1_y::CuArray{FloatType}
    u2_x::CuArray{FloatType}
    u2_y::CuArray{FloatType}
    u2_yλx::CuArray{FloatType}
    u2_xμy::CuArray{FloatType}
    u1_xλy::CuArray{FloatType}
    u1_yμx::CuArray{FloatType}
    u1_xλx::CuArray{FloatType}
    u1_xμx::CuArray{FloatType}
    u1_yμy::CuArray{FloatType}
    u2_yλy::CuArray{FloatType}
    u2_yμy::CuArray{FloatType}
    u2_xμx::CuArray{FloatType}

    v1::CuArray{FloatType}
    v2::CuArray{FloatType}
    u1::CuArray{FloatType}
    u2::CuArray{FloatType}

    λ::CuArray{FloatType}
    μ::CuArray{FloatType}
    D::CuArray{FloatType}
    iD::CuArray{FloatType}

    D1x::Union{CuArray{FloatType}, Nothing}
    D1y::Union{CuArray{FloatType}, Nothing}

    Dpx::Union{CuArray{FloatType}, Nothing}
    Dpy::Union{CuArray{FloatType}, Nothing}
    Dmx::Union{CuArray{FloatType}, Nothing}
    Dmy::Union{CuArray{FloatType}, Nothing}

    delta_2D::Union{CuArray{FloatType}, Nothing}

    D1x_∂Ω_W::Union{CuMatrix{FloatType}, Nothing}
    D1x_∂Ω_E::Union{CuMatrix{FloatType}, Nothing}
    D1y_∂Ω_N::Union{CuMatrix{FloatType}, Nothing}
    D1y_∂Ω_S::Union{CuMatrix{FloatType}, Nothing}

    Dpx_∂Ω_W::Union{CuMatrix{FloatType}, Nothing}
    Dpx_∂Ω_E::Union{CuMatrix{FloatType}, Nothing}
    Dpy_∂Ω_N::Union{CuMatrix{FloatType}, Nothing}
    Dpy_∂Ω_S::Union{CuMatrix{FloatType}, Nothing}

    Dmx_∂Ω_W::Union{CuMatrix{FloatType}, Nothing}
    Dmx_∂Ω_E::Union{CuMatrix{FloatType}, Nothing}
    Dmy_∂Ω_N::Union{CuMatrix{FloatType}, Nothing}
    Dmy_∂Ω_S::Union{CuMatrix{FloatType}, Nothing}

    d1x_l::Union{CuArray{FloatType}, Nothing}
    d1x_r::Union{CuArray{FloatType}, Nothing}
    d1y_l::Union{CuArray{FloatType}, Nothing}
    d1y_r::Union{CuArray{FloatType}, Nothing}

    d1x_l_wide::Union{CuArray{FloatType}, Nothing}
    d1x_r_wide::Union{CuArray{FloatType}, Nothing}
    d1y_l_wide::Union{CuArray{FloatType}, Nothing}
    d1y_r_wide::Union{CuArray{FloatType}, Nothing}

    D2x_λ_∂Ω_W_vec::Union{CuArray{M}, Nothing}
    D2x_λ_∂Ω_E_vec::Union{CuArray{M}, Nothing}
    D2x_μ_∂Ω_W_vec::Union{CuArray{M}, Nothing}
    D2x_μ_∂Ω_E_vec::Union{CuArray{M}, Nothing}
    D2x_λ_vec::Union{CuMatrix{V}, Nothing}
    D2x_μ_vec::Union{CuMatrix{V}, Nothing}

    D2y_λ_∂Ω_N_vec::Union{CuArray{M}, Nothing}
    D2y_λ_∂Ω_S_vec::Union{CuArray{M}, Nothing}
    D2y_μ_∂Ω_N_vec::Union{CuArray{M}, Nothing}
    D2y_μ_∂Ω_S_vec::Union{CuArray{M}, Nothing}
    D2y_μ_vec::Union{CuMatrix{V}, Nothing}
    D2y_λ_vec::Union{CuMatrix{V}, Nothing}
end

# Outer constructor that infers type parameters from the D2 stencil arrays
function DeviceData(
    u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx,
    u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
    v1, v2, u1, u2, λ, μ, D, iD,
    D1x, D1y, Dpx, Dpy, Dmx, Dmy, delta_2D,
    D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S,
    Dpx_∂Ω_W, Dpx_∂Ω_E, Dpy_∂Ω_N, Dpy_∂Ω_S,
    Dmx_∂Ω_W, Dmx_∂Ω_E, Dmy_∂Ω_N, Dmy_∂Ω_S,
    d1x_l, d1x_r, d1y_l, d1y_r,
    d1x_l_wide, d1x_r_wide, d1y_l_wide, d1y_r_wide,
    D2x_λ_∂Ω_W_vec, D2x_λ_∂Ω_E_vec, D2x_μ_∂Ω_W_vec, D2x_μ_∂Ω_E_vec,
    D2x_λ_vec, D2x_μ_vec,
    D2y_λ_∂Ω_N_vec, D2y_λ_∂Ω_S_vec, D2y_μ_∂Ω_N_vec, D2y_μ_∂Ω_S_vec,
    D2y_μ_vec, D2y_λ_vec
)
    # Infer M and V from the actual array element types, handling Nothing case
    M = D2x_λ_∂Ω_W_vec === nothing ? SMatrix{1,1,Float64,1} : eltype(D2x_λ_∂Ω_W_vec)
    V = D2x_λ_vec === nothing ? SVector{1,Float64} : eltype(D2x_λ_vec)
    
    DeviceData{M, V}(
        u1_x, u1_y, u2_x, u2_y, u2_yλx, u2_xμy, u1_xλy, u1_yμx,
        u1_xλx, u1_xμx, u1_yμy, u2_yλy, u2_yμy, u2_xμx,
        v1, v2, u1, u2, λ, μ, D, iD,
        D1x, D1y, Dpx, Dpy, Dmx, Dmy, delta_2D,
        D1x_∂Ω_W, D1x_∂Ω_E, D1y_∂Ω_N, D1y_∂Ω_S,
        Dpx_∂Ω_W, Dpx_∂Ω_E, Dpy_∂Ω_N, Dpy_∂Ω_S,
        Dmx_∂Ω_W, Dmx_∂Ω_E, Dmy_∂Ω_N, Dmy_∂Ω_S,
        d1x_l, d1x_r, d1y_l, d1y_r,
        d1x_l_wide, d1x_r_wide, d1y_l_wide, d1y_r_wide,
        D2x_λ_∂Ω_W_vec, D2x_λ_∂Ω_E_vec, D2x_μ_∂Ω_W_vec, D2x_μ_∂Ω_E_vec,
        D2x_λ_vec, D2x_μ_vec,
        D2y_λ_∂Ω_N_vec, D2y_λ_∂Ω_S_vec, D2y_μ_∂Ω_N_vec, D2y_μ_∂Ω_S_vec,
        D2y_μ_vec, D2y_λ_vec
    )
end

struct RHS_Cache{T}
    v1_t::T
    v2_t::T
    u1_t::T
    u2_t::T
end

struct RK4_Cache{T}
    # Intermediate stages (Pre-allocated!)
    k1v1::T; k1v2::T; k1u1::T; k1u2::T
    k2v1::T; k2v2::T; k2u1::T; k2u2::T
    k3v1::T; k3v2::T; k3u1::T; k3u2::T
    k4v1::T; k4v2::T; k4u1::T; k4u2::T
    
    # The RHS internal cache (what you already have)
    rhs_cache::RHS_Cache
end