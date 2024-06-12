
# Approximation error using the value of the pwa function
# in the given points. Support for multiple metrics (l1, l2, max).
function _approx_error(X::Matrix, z::Vector, pwa::PWAFunc, penalty = :l1)
    err = 0.0
    for (i, x̄) ∈ enumerate(eachcol(X))
        v = evaluate(pwa, x̄)

        if penalty == :l1
            err += abs(v - z[i])
        elseif penalty == :l2 || penalty == :rms
            err += (v - z[i])^2
        elseif penalty == :max
            err = max(err, abs(v - z[i]))
        end
    end

    if penalty == :rms
        err = err / size(X, 2)
    end

    if penalty == :l2 || penalty == :rms
        err = sqrt(err)
    end

    return err
end

function _convex_linearization_mb_single(
    X::Matrix,
    z::Vector,
    K,
    lᵐᵃˣ,
    penalty,
    optimizer,
    strict,
)
    𝒫 = _random_partition(X, K)
    pwa = _refine_partition(X, z, 𝒫, lᵐᵃˣ, penalty, optimizer, strict)
    e = _approx_error(X, z, pwa, penalty)
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
    penalty = options.pen
    strict = options.strict
    optimizer = options.optimizer

    @info "Starting heuristic search "
    approxes = collect(
        fetch.(
            @spawn _convex_linearization_mb_single(
                X,
                z,
                K,
                lᵐᵃˣ,
                penalty,
                optimizer,
                strict,
            ) for i ∈ 1:Nᵗʳ
        ),
    )
    min_error, pwa_best = argmin(first, approxes)

    @info "Terminating search - best approximation error = $(min_error) ($penalty)"
    return pwa_best
end

# Euclidean distance between two points
_dist(x, y) = sqrt(sum((x - y) .^ 2))

# Create a random partition of the points into K sets
# by generating K random points and group each data point to the
# closest of these random points.
function _random_partition(X, K)
    μ = vec(mean(X, dims = 2))
    σ² = cov(X, dims = 2)
    n = MvNormal(μ, σ²)
    P = rand(n, K)

    𝒫 = [[] for j ∈ 1:K]
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
function _refine_partition(X, z, 𝒫, lᵐᵃˣ, penalty, optimizer, strict)
    D = size(X, 1)
    pwa = nothing
    for it ∈ 1:lᵐᵃˣ
        pwa = PWAFunc{Convex,D}()
        for p ∈ 𝒫
            if length(p) > 0
                x̄ = X[:, p]
                z̄ = z[p]
                a, b = _local_fit(x̄, z̄, penalty, optimizer, strict)
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
# using the given penalty and possible restrictions on whether
# the plane should be strictly above or below the points
function _local_fit(X̄, z̄, penalty, optimizer, strict)
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

    if strict == :outer
        @constraint(m, ẑ .≥ z̄)
    elseif strict == :inner
        @constraint(m, ẑ .≤ z̄)
    end

    obj = AffExpr()
    if penalty == :l2 || penalty == :rms
        obj = sum((z̄[i] - ẑ[i])^2 for i ∈ 1:N)
    elseif penalty == :max
        @variable(m, t)
        obj = t
        for i ∈ 1:N
            @constraint(m, t ≥ (z̄[i] - ẑ[i]))
            @constraint(m, t ≥ (ẑ[i] - z̄[i]))
        end
    elseif penalty == :l1
        @variable(m, t[1:N])
        obj = sum(t[i] for i ∈ 1:N)
        for i ∈ 1:N
            @constraint(m, t[i] ≥ (z̄[i] - ẑ[i]))
            @constraint(m, t[i] ≥ (ẑ[i] - z̄[i]))
        end
    else
        error("Unrecognized/unsupported penalty type $(penalty)")
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

# The hyperplane being active at the point x for a convex pwa function 
function _active(pwa::PWAFunc{Convex,D}, x) where {D}
    return argmax(collect(evaluate(p, x) for p ∈ pwa.planes))
end

# Creating an updated partition by associating each data point to
# the hyperplane being active for the data point
function _update_partition(X, pwa)
    𝒫 = [[] for _ ∈ 1:_planes(pwa)]
    for (i, x) ∈ enumerate(eachcol(X))
        push!(𝒫[_active(pwa, x)], i)
    end
    return 𝒫
end
