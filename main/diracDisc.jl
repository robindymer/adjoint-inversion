"""
    diracDiscr(g::Grid, x_s::Vector{Float64}, m_order::Int, s_order::Int, H::Union{quadrature{Float64}, Vector{quadrature{Float64}}})

n-dimensional delta function discretization

# Arguments
- `g::Grid`: Cartesian grid
- `x_s::Vector{Float64}`: source point coordinate vector, e.g. [x, y] or [x, y, z]
- `m_order::Int`: Number of moment conditions
- `s_order::Int`: Number of smoothness conditions
- `H::Union{quadrature{Float64}, Vector{quadrature{Float64}}}`: Norm matrix (or vector of 1D norm matrices)

# Returns
- `d::Vector{Float64}`: Discretized delta function
"""
function diracDiscr(g::Grid, x_s::Vector{Float64}, m_order::Int, s_order::Int, quadrature::Union{Vector{Float64}, Vector{Vector{Float64}}})
    dim = 2  # TODO: Fix hardcoded dimension
    d_1D = Vector{Vector{Float64}}(undef, dim)

    # Allow for non-vector input in 1D
    quadrature_vec = quadrature isa Vector ? quadrature : [quadrature for _ in 1:dim]
    
    # Create 1D dirac discr for each coordinate direction
    # for i in 1:dim
    d_1D[1] = diracDiscr1D(x_s[1], g.x_values, m_order, s_order, quadrature_vec[1])
    d_1D[2] = diracDiscr1D(x_s[2], g.y_values, m_order, s_order, quadrature_vec[2])
    # end

    # Create 2D delta function using outer product: d_1D[1] ⊗ d_1D[2]
    # This creates a matrix where delta[i,j] corresponds to x[i], y[j]
    # vec() then vectorizes in column-major order (x varies fastest)
    d = d_1D[1] * d_1D[2]'
    d = vec(d)  # Turn into column vector
    
    # # Start with the last dimension
    # d = d_1D[dim]
    
    # # Perform outer products from dim-1 down to 1
    # for i in dim-1:-1:1
    #     # Perform outer product, transpose, and then turn into column vector
    #     d = (d_1D[i] * d')'
    #     d = vec(d)  # Turn into column vector
    # end
    
    return d
end


"""
    diracDiscr1D(x_s::Float64, x::Vector{Float64}, m_order::Int, s_order::Int, H::Matrix{Float64})

Helper function for 1D delta functions

# Arguments
- `x_s::Float64`: Source location
- `x::Vector{Float64}`: Grid points
- `m_order::Int`: Number of moment conditions
- `s_order::Int`: Number of smoothness conditions
- `H::Matrix{Float64}`: Norm matrix

# Returns
- `ret::Vector{Float64}`: 1D discretized delta function
"""
function diracDiscr1D(x_s::Float64, x::Vector{Float64}, m_order::Int, s_order::Int, quadrature::Vector{Float64})
    # Return zeros if x_s is outside grid
    if x_s < x[1] || x_s > x[end]
        return zeros(length(x))
    end
    
    tot_order = m_order + s_order  # This is equiv. to the number of equations solved for

    # Get interior grid spacing
    middle = floor(Int, length(x) / 2)
    h = x[middle + 1] - x[middle]  # Use middle point to allow for staggered grids

    index = sourceIndices(x_s, x, tot_order, h)

    polynomial = (x[index] .- x[index[1]]) / (x[index[end]] - x[index[1]])
    x_0 = (x_s - x[index[1]]) / (x[index[end]] - x[index[1]])

    quadrature_weights = quadrature[index] / h

    h_polynomial = polynomial[2] - polynomial[1]
    
    # Initialize moment equations right-hand side
    b = zeros(tot_order)
    for i in 1:m_order
        b[i] = x_0^(i - 1)
    end

    # Build moment matrix M
    M = zeros(m_order, tot_order)
    for i in 1:tot_order
        for j in 1:m_order
            M[j, i] = polynomial[i]^(j - 1) * h_polynomial * quadrature_weights[i]
        end
    end

    # Build smoothness matrix S
    S = zeros(s_order, tot_order)
    for i in 1:tot_order
        for j in 1:s_order
            S[j, i] = (-1)^(i - 1) * polynomial[i]^(j - 1)
        end
    end

    # Combine into single system
    A = vcat(M, S)

    # Solve for delta function coefficients
    d = A \ b
    
    # Create return vector
    ret = zeros(length(x))
    ret[index] = d / h * h_polynomial
    
    return ret
end


"""
    sourceIndices(x_s::Float64, x::Vector{Float64}, tot_order::Int, h::Float64)

Find the indices that are within range of the point source location

# Arguments
- `x_s::Float64`: Source location
- `x::Vector{Float64}`: Grid points
- `tot_order::Int`: Total order (moment + smoothness)
- `h::Float64`: Grid spacing

# Returns
- `I::Vector{Int}`: Indices of grid points to use for delta function
"""
function sourceIndices(x_s::Float64, x::Vector{Float64}, tot_order::Int, h::Float64)
    # Find the indices that are within range of the point source location
    I = findall(abs.(x .- x_s) .<= tot_order * h / 2)

    if length(I) > tot_order
        if length(I) == tot_order + 2
            I = I[2:end-1]
        elseif length(I) == tot_order + 1
            I = I[1:end-1]
        end
    elseif length(I) < tot_order
        if x_s < x[1] + ceil(tot_order / 2) * h
            I = collect(1:tot_order)
        elseif x_s > x[end] - ceil(tot_order / 2) * h
            I = collect(length(x) - tot_order + 1:length(x))
        else
            if I[end] < length(x)
                push!(I, I[end] + 1)
            else
                pushfirst!(I, I[1] - 1)
            end
        end
    end
    
    return I
end