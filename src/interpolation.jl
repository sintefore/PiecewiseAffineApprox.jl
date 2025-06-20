# A general 1D piecewise linear function defined by points 
# in sorted order and their associated function value 
struct PWAFunction1D
    x::Vector{Float64}
    z::Vector{Float64}

    function PWAFunction1D(x::Vector{Float64}, z::Vector{Float64})
        @assert length(x) == length(z)
        @assert issorted(x)
        return new(x, z)
    end
end

function PWAFunction1D(x::Vector, z::Vector)
    return PWAFunction1D(
        convert(Vector{Float64}, x),
        convert(Vector{Float64}, z),
    )
end

function _asconvex(pwa::PWAFunction1D)
    N = length(pwa.x)
    x = pwa.x
    z = pwa.z

    c = [(z[i+1] - z[i]) / (x[i+1] - x[i]) for i ∈ 1:(N-1)]
    d = [z[i] - c[i] * x[i] for i ∈ 1:(N-1)]

    return PWAFunc{Convex,1}(Plane.(c, d))
end

# Create a convex approximation to a general pwa function
# by perturbing function values as little as possible (l1-deviation)
function _convexify1D(pwa::PWAFunction1D, options)
    N = length(pwa.x)
    𝒩 = 1:N
    x = pwa.x
    z = pwa.z

    m = Model()

    δ⁺ = @variable(m, [𝒩], lower_bound = 0)
    δ⁻ = @variable(m, [𝒩], lower_bound = 0)

    @objective(m, Min, sum(δ⁺[i] + δ⁻[i] for i ∈ 𝒩))
    for i ∈ 2:(N-1)
        @constraint(
            m,
            (z[i] + δ⁺[i] - δ⁻[i] - z[i-1] - δ⁺[i-1] + δ⁻[i-1]) /
            (x[i] - x[i-1]) <=
            (z[i+1] + δ⁺[i+1] - δ⁻[i+1] - z[i] - δ⁺[i] + δ⁻[i]) /
            (x[i+1] - x[i])
        )
    end

    set_optimizer(m, options.optimizer)
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    dp = value.(δ⁺)
    dn = value.(δ⁻)

    return PWAFunction1D(x, [z[i] + dp[i] - dn[i] for i ∈ 𝒩])
end

_convexify1D(x, z, options) = _convexify1D(PWAFunction1D(x, z), options)

# Create a pwa interpolant to the given points with a maximum number of 
# segments 
function _interpolatepw(x, z, options)
    @assert(length(x) == length(z))

    N = length(x)
    𝒩 = 1:N

    # Find slopes
    c = [(z[j] - z[i]) / (x[j] - x[i]) for i ∈ 𝒩, j ∈ 𝒩]

    # Calculate penalties
    if options.metric == :l1
        p = [
            (
                i < j ?
                sum(abs(c[i, j] * (x[k] - x[i]) + z[i] - z[k]) for k ∈ i:j) : 0
            ) for i ∈ 𝒩, j ∈ 𝒩
        ]
    elseif options.metric == :l2
        p = [
            (
                i < j ?
                sum((c[i, j] * (x[k] - x[i]) + z[i] - z[k])^2 for k ∈ i:j) :
                0
            ) for i ∈ 𝒩, j ∈ 𝒩
        ]
    elseif options.metric == :max
        p = [
            (
                i < j ?
                maximum(
                    abs(c[i, j] * (x[k] - x[i]) + z[i] - z[k]) for k ∈ i:j
                ) : 0
            ) for i ∈ 𝒩, j ∈ 𝒩
        ]
    else
        error("Unrecognized/unsupported metric $(options.metric)")
    end

    m = Model()
    @variable(m, 𝑢[𝒩, 𝒩], Bin)

    # Minimize total metric
    @objective(m, Min, sum(p[i, j] * 𝑢[i, j] for i ∈ 𝒩, j ∈ 𝒩))

    # Number of line segments in interpolant
    @constraint(m, sum(𝑢[i, j] for i ∈ 𝒩, j ∈ 𝒩) == options.planes)

    # Only forward segments allowed
    for i ∈ 𝒩, j ∈ 𝒩
        if i >= j
            @constraint(m, 𝑢[i, j] == 0)
        end
    end

    # Path modelling
    for j ∈ 𝒩
        @constraint(
            m,
            (j > 1 ? sum(𝑢[i, j] for i ∈ 𝒩) : 1) ==
            (j < N ? sum(𝑢[j, i] for i ∈ 𝒩) : 1)
        )
    end

    set_optimizer(m, options.optimizer)
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    𝑢ᴼᵖᵗ = value.(m[:𝑢])
    xᴼᵖᵗ = collect(x[i] for i ∈ 𝒩, j ∈ 𝒩 if 𝑢ᴼᵖᵗ[i, j] == 1)
    push!(xᴼᵖᵗ, x[N])
    zᴼᵖᵗ = collect(z[i] for i ∈ 𝒩, j ∈ 𝒩 if 𝑢ᴼᵖᵗ[i, j] == 1)
    push!(zᴼᵖᵗ, z[N])

    return PWAFunction1D(xᴼᵖᵗ, zᴼᵖᵗ)
end

function _convex_linearization_ipol(x, z, options)
    return _asconvex(_convexify1D(_interpolatepw(x, z, options), options))
end
