
# Approximation error using the value of the pwa function
# in the given points. Support for multiple metrics (l1, l2, max).
function _approx_error(X::Matrix, z::Vector, pwa::PWAFunc, metric = :l1)
    err = 0.0
    for (i, x̄) ∈ enumerate(eachcol(X))
        v = evaluate(pwa, x̄)

        if metric == :l1
            err += abs(v - z[i])
        elseif metric == :l2 || metric == :rms
            err += (v - z[i])^2
        elseif metric == :max
            err = max(err, abs(v - z[i]))
        end
    end

    if metric == :rms
        err = err / size(X, 2)
    end

    if metric == :l2 || metric == :rms
        err = sqrt(err)
    end

    return err
end

function _convex_linearization_mb_single(
    X::Matrix,
    z::Vector,
    K,
    lᵐᵃˣ,
    metric,
    optimizer,
    strict,
    rng,
)
    𝒫 = _random_partition(rng, X, K)
    pwa = _refine_partition(X, z, 𝒫, lᵐᵃˣ, metric, optimizer, strict)
    e = _approx_error(X, z, pwa, metric)
    return e => pwa
end

# Finds a pwa convex approximation for the provided data
# X is a matrix with a column for each data point and z is a vector with the 
# corresponding function values
function _convex_linearization_mb(X::Matrix, z::Vector, options::Cluster)
    @assert(size(X, 2) == length(z))

    Nᵗʳ = options.trials    # Number of trials
    lᵐᵃˣ = options.itlim    # Iteration limit
    K = options.planes
    metric = options.metric
    strict = options.strict
    optimizer = options.optimizer

    if isnothing(optimizer) && !(strict == :none && metric ∈ (:l2, :rms))
        error(
            "Cluster with metric = :$metric and strict = :$strict requires an optimizer. " *
            "Only metric = :l2 with strict = :none can be used without one.",
        )
    end

    # Work on normalized data to make the heuristic invariant to affine
    # transformations of the input (Magnani & Boyd, Sect. 1.3)
    μ, T = _normalization(X)
    X̃ = T * (X .- μ)

    # Draw one seed per trial up front so that the result does not depend on
    # how the trials are scheduled across threads
    seeds = rand(options.rng, UInt64, Nᵗʳ)

    @debug "Starting heuristic search "
    approxes = collect(
        fetch.(
            @spawn _convex_linearization_mb_single(
                X̃,
                z,
                K,
                lᵐᵃˣ,
                metric,
                optimizer,
                strict,
                StableRNG(seeds[i]),
            ) for i ∈ 1:Nᵗʳ
        ),
    )
    min_error, pwa_best = argmin(first, approxes)

    @debug "Terminating search - best approximation error = $(min_error) ($metric)"
    return _denormalize(pwa_best, μ, T)
end

# Affine transformation x̃ = T(x - μ) giving data with zero mean and identity
# covariance. T = Σ^(-1/2) (ZCA whitening) is used since it is unique and leaves
# data that is only scaled along the axes in identical normalized coordinates.
# Directions with (near) zero variance are mapped to zero.
function _normalization(X)
    μ = vec(mean(X, dims = 2))
    λ, V = eigen(Symmetric(cov(X, dims = 2)))
    keep = λ .> 1e-10 * max(maximum(λ), 0)
    if !any(keep)
        # All points coincide, nothing to normalize
        return μ, Matrix{Float64}(I, length(μ), length(μ))
    end
    Vᵣ = V[:, keep]
    return μ, Vᵣ * Diagonal(1 ./ sqrt.(λ[keep])) * Vᵣ'
end

# Transform a pwa function fitted to normalized data x̃ = T(x - μ) back to the
# original coordinates: ãᵀx̃ + b̃ = (Tᵀã)ᵀx + (b̃ - ãᵀTμ)
function _denormalize(pwa::PWAFunc{Convex,D}, μ, T) where {D}
    planes = map(pwa.planes) do p
        ã = collect(p.α)
        return Plane(T' * ã, p.β - dot(ã, T * μ))
    end
    return PWAFunc{Convex,D}(planes)
end

# Euclidean distance between two points
_dist(x, y) = sqrt(sum((x - y) .^ 2))

# Draw K samples from a normal distribution with the same mean and covariance
# as the columns of X. The covariance is factored with an eigendecomposition
# rather than Cholesky so that singular covariances (e.g. collinear points) work.
function _sample_normal(rng::AbstractRNG, X, K)
    μ = vec(mean(X, dims = 2))
    λ, V = eigen(Symmetric(cov(X, dims = 2)))
    L = V * Diagonal(sqrt.(max.(λ, 0)))
    return μ .+ L * randn(rng, length(μ), K)
end

# Create a random partition of the points into K sets
# by generating K random points and group each data point to the
# closest of these random points.
function _random_partition(rng::AbstractRNG, X, K)
    P = _sample_normal(rng, X, K)

    𝒫 = [Int[] for j ∈ 1:K]
    for (i, x) ∈ enumerate(eachcol(X))
        # Find the nearest point amongst the p's
        jmin = argmin(_dist(x, P[:, j]) for j ∈ 1:K)
        push!(𝒫[jmin], i)
    end
    return 𝒫
end

# Improve a partition by fitting a hyperplane for all points 
# belonging to each subset of the partition and then 
# updating the partition such that each point is associated with the
# plane being active at the point. 
#
# The process terminates after a given number of iterations returning
# the pwa function corresponding to the final partition.
function _refine_partition(X, z, 𝒫, lᵐᵃˣ, metric, optimizer, strict)
    D = size(X, 1)
    pwa = nothing
    for it ∈ 1:lᵐᵃˣ
        pwa = PWAFunc{Convex,D}()
        for p ∈ 𝒫
            if length(p) > 0
                x̄ = X[:, p]
                z̄ = z[p]
                a, b = _local_fit(x̄, z̄, metric, optimizer, strict)
                _addplane!(pwa, a, b)
            end
        end
        𝒫ⁿᵉʷ = _update_partition(X, pwa)
        if 𝒫ⁿᵉʷ == 𝒫
            break
        end
        𝒫 = 𝒫ⁿᵉʷ
    end

    return pwa
end

# Finds the hypeplane best fitting the given subset of points 
# using the given metric and possible restrictions on whether
# the plane should be strictly above or below the points
function _local_fit(X̄, z̄, metric, optimizer, strict)
    if strict == :none && (metric == :l2 || metric == :rms)
        return _least_squares_fit(X̄, z̄)
    end

    M, N = size(X̄)

    # Create an optimization model to find the best a and b such that  ax + b ≈ z
    m = Model()

    if Threads.nthreads() > 1
        # Use threading for individual problems, avoid oversubscribing by limiting threads for each problem
        MOI.set(m, MOI.NumberOfThreads(), 1)
    end

    @variable(m, a[1:M])
    @variable(m, b)
    @variable(m, ẑ[1:N])

    for i ∈ 1:N
        @constraint(m, ẑ[i] == sum(a[j] * X̄[j, i] for j ∈ 1:M) + b)
    end

    if strict == :inner
        @constraint(m, ẑ .≥ z̄)
    elseif strict == :outer
        @constraint(m, ẑ .≤ z̄)
    end

    obj = AffExpr()
    if metric == :l2 || metric == :rms
        obj = sum((z̄[i] - ẑ[i])^2 for i ∈ 1:N)
    elseif metric == :max
        @variable(m, t)
        obj = t
        for i ∈ 1:N
            @constraint(m, t ≥ (z̄[i] - ẑ[i]))
            @constraint(m, t ≥ (ẑ[i] - z̄[i]))
        end
    elseif metric == :l1
        @variable(m, t[1:N])
        obj = sum(t[i] for i ∈ 1:N)
        for i ∈ 1:N
            @constraint(m, t[i] ≥ (z̄[i] - ẑ[i]))
            @constraint(m, t[i] ≥ (ẑ[i] - z̄[i]))
        end
    else
        error("Unrecognized/unsupported metric type $(metric)")
    end

    # TODO: consider adding regularization term

    @objective(m, Min, obj)

    set_optimizer(m, optimizer)
    set_silent(m)
    optimize!(m)

    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    ā = JuMP.value.(a)
    b̄ = JuMP.value.(b)

    return ā, b̄
end

# Closed form least squares fit of the plane aᵀx + b to the points (Magnani & Boyd, eq. 9).
# Solved via the centered normal equations, with a small ridge term on a to handle
# subsets that do not span all directions (e.g. fewer than D+1 points). This gives
# zero slope along directions with no variation in the data.
function _least_squares_fit(X̄, z̄)
    x̄ = vec(mean(X̄, dims = 2))
    Xc = X̄ .- x̄
    zc = z̄ .- mean(z̄)
    C = Xc * Xc'
    λ = 1e-8 * tr(C) / size(C, 1)
    a = λ > 0 ? (C + λ * I) \ (Xc * zc) : zeros(size(C, 1))
    b = mean(z̄) - dot(a, x̄)
    return a, b
end

# The hyperplane being active at the point x for a convex pwa function
function _active(pwa::PWAFunc{Convex,D}, x) where {D}
    return argmax(collect(evaluate(p, x) for p ∈ pwa.planes))
end

# Creating an updated partition by associating each data point to
# the hyperplane being active for the data point
function _update_partition(X, pwa)
    𝒫 = [Int[] for _ ∈ 1:_planes(pwa)]
    for (i, x) ∈ enumerate(eachcol(X))
        push!(𝒫[_active(pwa, x)], i)
    end
    return 𝒫
end
