
function _convex_linearization_mb_single(
    X::Matrix,
    z::Vector,
    K,
    lᵐᵃˣ,
    metric,
    optimizer,
    strict,
    rng;
    polish = 5,
)
    𝒫 = _random_partition(rng, X, K)
    if strict == :none && metric ∈ (:l2, :rms)
        pwa = _refine_partition(X, z, 𝒫, lᵐᵃˣ, metric, optimizer, strict)
    else
        # Two phases: the partition is first refined with closed form least
        # squares, which is fast and moves the planes smoothly, and then for a
        # few iterations with the requested fit. This is both faster and better
        # than the requested fit alone, whose solutions at vertices of the
        # feasible region tend to get stuck in limit cycles early.
        pwa₁ = _refine_partition(X, z, 𝒫, lᵐᵃˣ, :l2, nothing, :none)
        𝒫 = _update_partition(X, pwa₁)
        l = min(lᵐᵃˣ, polish)
        pwa = _refine_partition(X, z, 𝒫, l, metric, optimizer, strict)
    end
    e = approximation_error(pwa, X, z, metric)
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

    _check_option(options, "metric", metric, (:l1, :l2, :rms, :max))
    _check_option(options, "strict", strict, (:none, :inner, :outer))
    K ≥ 1 || throw(ArgumentError("Cluster requires at least one plane"))
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
    tasks = [
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
    ]
    approxes = map(tasks) do task
        try
            fetch(task)
        catch e
            # Report the error of the trial rather than the failure of the task
            e isa TaskFailedException ? throw(e.task.result) : rethrow()
        end
    end
    min_error, pwa_best = argmin(first, approxes)

    @debug "Terminating search - best approximation error = $(min_error) ($metric)"
    return _denormalize(pwa_best, μ, T)
end

# Affine transformation x̃ = T(x - μ) giving data with zero mean and identity
# covariance. Each coordinate is first scaled by its range, so that coordinates
# with a small but genuine range are not treated as degenerate, and the result
# is then whitened with T = Σ^(-1/2) (ZCA whitening), which is unique and leaves
# data that is only scaled along the axes in identical normalized coordinates.
# Directions with (near) zero variance are mapped to zero.
function _normalization(X)
    D, N = size(X)
    μ = vec(mean(X, dims = 2))
    N < 2 && return μ, Matrix{Float64}(I, D, D)
    r = vec(maximum(X, dims = 2) .- minimum(X, dims = 2))
    S = Diagonal([rᵢ > 0 ? 1 / rᵢ : 1.0 for rᵢ ∈ r])
    λ, V = eigen(Symmetric(cov(S * (X .- μ), dims = 2)))
    keep = λ .> 1e-10 * max(maximum(λ), 0)
    if !any(keep)
        # All points coincide, nothing to normalize
        return μ, Matrix{Float64}(I, D, D)
    end
    Vᵣ = V[:, keep]
    return μ, Vᵣ * Diagonal(1 ./ sqrt.(λ[keep])) * Vᵣ' * S
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
    size(X, 2) < 2 && return repeat(μ, 1, K)
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
# Subsets that do not spread out in all directions spanned by the data (e.g. fewer
# than D+1 points, or points lying in a lower dimensional subspace) do not determine
# a plane. Fitting one anyway can give a plane that dominates all others (e.g. a
# flat plane through a single high point). Such subsets are dropped and their
# points reassigned to the active plane among the remaining ones. To keep the number of planes,
# the partition is then refilled by splitting the subsets with the largest error
# (deviates from Magnani & Boyd, who drop empty subsets and continue with
# fewer planes).
#
# The process terminates when the partition no longer changes, or after a given
# number of iterations, returning the best pwa function found (as measured by
# metric). Only fits to complete partitions are considered, to retain the
# properties given by strict.
function _refine_partition(
    X,
    z,
    𝒫,
    lᵐᵃˣ,
    metric,
    optimizer,
    strict;
    refill = true,
    minvar = 1e-4,
)
    D, N = size(X)
    K = length(𝒫)
    determined = _determined_test(X, minvar)
    𝒫 = filter(determined, 𝒫)
    if isempty(𝒫)
        𝒫 = [collect(1:N)]
    end
    refill && _refill!(𝒫, X, ones(N), K, determined)
    fit(q) = _fit_planes(X, z, q, metric, optimizer, strict)

    pwaᵇᵉˢᵗ = fit(𝒫)
    eᵇᵉˢᵗ = approximation_error(pwaᵇᵉˢᵗ, X, z, metric)
    pwa = pwaᵇᵉˢᵗ
    for it ∈ 1:lᵐᵃˣ
        if it > 1
            pwa = fit(𝒫)
            e = approximation_error(pwa, X, z, metric)
            if e < eᵇᵉˢᵗ
                pwaᵇᵉˢᵗ, eᵇᵉˢᵗ = pwa, e
            end
        end

        # Update the partition with the active planes, reassigning points
        # in undetermined subsets to the remaining planes
        𝒫ⁿᵉʷ = _update_partition(X, pwa)
        active = pwa
        if !all(determined, 𝒫ⁿᵉʷ)
            keep = determined.(𝒫ⁿᵉʷ)
            if any(keep)
                active = PWAFunc{Convex,D}(pwa.planes[keep])
                𝒫ⁿᵉʷ = _update_partition(X, active)
            end
            filter!(p -> !isempty(p), 𝒫ⁿᵉʷ)
        end
        if refill && length(𝒫ⁿᵉʷ) < K
            r = [evaluate(active, x) - z[i] for (i, x) ∈ enumerate(eachcol(X))]
            w = metric ∈ (:l2, :rms) ? r .^ 2 : abs.(r)
            _refill!(𝒫ⁿᵉʷ, X, w, K, determined)
        end
        if 𝒫ⁿᵉʷ == 𝒫
            # Converged, pwa is the fit to this partition
            return pwaᵇᵉˢᵗ
        end
        𝒫 = 𝒫ⁿᵉʷ
    end

    # Iteration limit reached, also consider the fit to the final partition
    pwa = fit(𝒫)
    return approximation_error(pwa, X, z, metric) < eᵇᵉˢᵗ ? pwa : pwaᵇᵉˢᵗ
end

# Fit a plane to each subset of the partition 𝒫. For strict = :outer, all planes
# are required to lie below all data points, so that the pwa function is an outer
# approximation also when the partition has not converged.
function _fit_planes(X, z, 𝒫, metric, optimizer, strict)
    pwa = PWAFunc{Convex,size(X, 1)}()
    for p ∈ 𝒫
        a, b = _local_fit(X[:, p], z[p], metric, optimizer, strict, X, z)
        _addplane!(pwa, a, b)
    end
    return pwa
end

# Returns a function testing whether a subset of the points (columns of X)
# determines a plane: its sample variance must be at least minvar in all directions
# where the full data set has variation. X is assumed to be normalized, i.e. to
# have unit variance in these directions.
function _determined_test(X, minvar)
    r = if size(X, 2) < 2
        0
    else
        count(>(0.5), eigvals(Symmetric(cov(X, dims = 2))))
    end
    return function (p)
        length(p) > r || return false
        r == 0 && return true
        λ = eigvals(Symmetric(cov(X[:, p], dims = 2, corrected = false)))
        return count(≥(minvar), λ) ≥ r
    end
end

# Increase the number of subsets in the partition 𝒫 to K by repeatedly splitting
# the subset with the largest total weight w (e.g. approximation error), such that
# all subsets still determine a plane.
function _refill!(𝒫, X, w, K, determined)
    splittable = trues(length(𝒫))
    while length(𝒫) < K
        candidates = [j for j ∈ eachindex(𝒫) if splittable[j]]
        isempty(candidates) && break
        j = argmax(j -> sum(w[𝒫[j]]), candidates)
        p₁, p₂ = _split_principal(X, 𝒫[j])
        if !determined(p₁) || !determined(p₂)
            splittable[j] = false
        else
            𝒫[j] = p₁
            push!(𝒫, p₂)
            push!(splittable, true)
        end
    end
    return 𝒫
end

# Split the points p in two halves at the median along their principal direction
function _split_principal(X, p)
    Xc = X[:, p] .- mean(X[:, p], dims = 2)
    d = eigen(Symmetric(Xc * Xc')).vectors[:, end]
    s = vec(d' * Xc)
    m = median(s)
    return p[s .< m], p[s .>= m]
end

# Finds the hypeplane best fitting the given subset of points 
# using the given metric and possible restrictions on whether
# the plane should be strictly above the points (strict = :inner) or
# below the points Xᵒ with values zᵒ (strict = :outer)
function _local_fit(X̄, z̄, metric, optimizer, strict, Xᵒ = X̄, zᵒ = z̄)
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
        for i ∈ axes(Xᵒ, 2)
            @constraint(m, sum(a[j] * Xᵒ[j, i] for j ∈ 1:M) + b ≤ zᵒ[i])
        end
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

    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) ($(raw_status(m)))",
        )
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

function _fit(a::Cluster, input::FunctionEvaluations{D}) where {D}
    x = [p[i] for i ∈ 1:D, p ∈ input.points]
    return _convex_linearization_mb(x, input.values, a)
end
