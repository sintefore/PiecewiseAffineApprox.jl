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
    if !is_solved_and_feasible(m)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end

    dp = value.(δ⁺)
    dn = value.(δ⁻)

    return PWAFunction1D(x, [z[i] + dp[i] - dn[i] for i ∈ 𝒩])
end

_convexify1D(x, z, options) = _convexify1D(PWAFunction1D(x, z), options)

# Create a pwa interpolant to the given points (sorted by x) using exactly
# options.planes segments with end points at data points. The segments are
# chosen to minimize the total deviation in options.metric: the sum over the
# segments for :l1 and :l2, and the maximum over the segments for :max.
function _interpolatepw(x, z, options)
    @assert length(x) == length(z)
    p = _segment_penalties(x, z, options.metric)
    minimax = options.metric == :max
    knots, _ = _shortest_path_segments(p, options.planes; minimax)
    return PWAFunction1D(x[knots], z[knots])
end

# Penalty p[i, j] (i < j) for approximating the data points (xₖ, zᵈₖ),
# k = i, ..., j-1, by the line segment between the knots (xᵢ, zᵏᵢ) and
# (xⱼ, zᵏⱼ):
# the sum of absolute deviations (:l1), the sum of squared deviations (:l2) or
# the maximum absolute deviation (:max). Each data point is thereby assigned to
# exactly one segment, except the last point that is always a knot.
function _segment_penalties(x, zᵏ, zᵈ, metric)
    if metric == :l2
        return _segment_penalties_l2(x, zᵏ, zᵈ)
    elseif metric ∉ (:l1, :max)
        error("Unrecognized/unsupported metric $(metric)")
    end
    N = length(x)
    p = zeros(N, N)
    for i ∈ 1:N-1, j ∈ i+1:N
        c = (zᵏ[j] - zᵏ[i]) / (x[j] - x[i])
        acc = 0.0
        for k ∈ i:j-1
            r = abs(zᵏ[i] + c * (x[k] - x[i]) - zᵈ[k])
            acc = metric == :l1 ? acc + r : max(acc, r)
        end
        p[i, j] = acc
    end
    return p
end
_segment_penalties(x, z, metric) = _segment_penalties(x, z, z, metric)

# Sum of squared deviations for all segments in O(N²) using prefix sums.
# The data is centered first to reduce cancellation errors.
function _segment_penalties_l2(x, zᵏ, zᵈ)
    N = length(x)
    x = x .- mean(x)
    zᵏ = zᵏ .- mean(zᵈ)
    zᵈ = zᵈ .- mean(zᵈ)
    Sx = [0.0; cumsum(x)]
    Sz = [0.0; cumsum(zᵈ)]
    Sxx = [0.0; cumsum(x .^ 2)]
    Szz = [0.0; cumsum(zᵈ .^ 2)]
    Sxz = [0.0; cumsum(x .* zᵈ)]
    p = zeros(N, N)
    for i ∈ 1:N-1, j ∈ i+1:N
        # Deviations at points k = i, ..., j-1 are α + c xₖ - zᵈₖ
        c = (zᵏ[j] - zᵏ[i]) / (x[j] - x[i])
        α = zᵏ[i] - c * x[i]
        n = j - i
        sx, sz = Sx[j] - Sx[i], Sz[j] - Sz[i]
        sxx, szz, sxz = Sxx[j] - Sxx[i], Szz[j] - Szz[i], Sxz[j] - Sxz[i]
        r² =
            n * α^2 + c^2 * sxx + szz + 2α * c * sx - 2α * sz - 2c * sxz
        p[i, j] = max(r², 0.0)
    end
    return p
end

# Find the path 1 → N through increasing indices with exactly K edges that
# minimizes the sum of edge costs p[i, j], or with minimax = true the maximum
# edge cost. Returns the indices on the path and the minimal cost of reaching
# N using 1, ..., K edges.
#
# For minimax, ties are broken by the sum of edge costs: a second pass
# minimizes the sum using only edges with cost at most the minimax value.
function _shortest_path_segments(p, K; minimax = false)
    if !minimax
        return _dp_segments(p, K, +)
    end
    _, costs = _dp_segments(p, K, max)
    q = [v ≤ costs[K] ? v : Inf for v ∈ p]
    path, _ = _dp_segments(q, K, +)
    return path, costs
end

# Dynamic program in O(KN²) for the path 1 → N with exactly K edges that
# minimizes the edge costs p[i, j] aggregated by op (+ or max). Ties are broken
# towards smaller indices.
function _dp_segments(p, K, op)
    cost, pred = _dp_tables(p, K, op)
    N = size(p, 1)
    return _dp_backtrack(pred, K, N), cost[2:end, N]
end

# The optimal paths with exactly s = 1, ..., K edges
function _dp_all_paths(p, K, op)
    _, pred = _dp_tables(p, K, op)
    return [_dp_backtrack(pred, s, size(p, 1)) for s ∈ 1:K]
end

function _dp_tables(p, K, op)
    N = size(p, 1)
    @assert 1 ≤ K ≤ N - 1
    # cost[s+1, j]: minimal cost of reaching j using s edges
    cost = fill(Inf, K + 1, N)
    pred = zeros(Int, K + 1, N)
    cost[1, 1] = 0.0
    for s ∈ 1:K, j ∈ s+1:N
        for i ∈ s:j-1
            v = op(cost[s, i], p[i, j])
            if v < cost[s+1, j]
                cost[s+1, j] = v
                pred[s+1, j] = i
            end
        end
    end
    return cost, pred
end

function _dp_backtrack(pred, s, N)
    path = [N]
    for t ∈ s:-1:1
        pushfirst!(path, pred[t+1, path[1]])
    end
    return path
end

# Find the path 1 → N through increasing indices with at most K edges that
# minimizes the edge costs p[i, j] (aggregated by + or, with minimax = true,
# by max), subject to the slopes σ[i, j] being nondecreasing along the path.
# With slopes of the line segments between data points, this gives the best
# convex interpolant through a subset of the data points.
#
# The state is the last edge (i, j). For each j, the incoming edges and the
# outgoing edges are sorted by slope, so that the best feasible incoming edge
# for all outgoing edges is found in a single sweep. This gives O(KN²) time
# after sorting. Returns an empty path if no convex path exists.
function _convex_shortest_path(σ, p, K; minimax = false)
    if !minimax
        return _convex_dp_segments(σ, p, K, +)
    end
    path = _convex_dp_segments(σ, p, K, max)
    isempty(path) && return path
    M = maximum(p[path[k], path[k+1]] for k ∈ 1:length(path)-1)
    q = [v ≤ M ? v : Inf for v ∈ p]
    return _convex_dp_segments(σ, q, K, +)
end

function _convex_dp_segments(σ, p, K, op)
    best, last, pred = _convex_dp_tables(σ, p, K, op)
    s = argmin(best)
    best[s] == Inf && return Int[]
    return _convex_backtrack(pred, s, last[s], size(p, 1))
end

# The optimal convex paths with at most k = 1, ..., K edges
function _convex_dp_all_paths(σ, p, K, op)
    best, last, pred = _convex_dp_tables(σ, p, K, op)
    N = size(p, 1)
    return map(1:min(K, N - 1)) do k
        s = argmin(best[1:k])
        return _convex_backtrack(pred, s, last[s], N)
    end
end

# Returns the minimal cost of reaching N with exactly s edges, the first node
# of the corresponding last edge, and the predecessor table
function _convex_dp_tables(σ, p, K, op)
    N = size(p, 1)
    K = min(K, N - 1)
    # cost[i, j]: minimal cost of a path to j with s edges, the last (i, j)
    cost = fill(Inf, N, N)
    for j ∈ 2:N
        cost[1, j] = p[1, j]
    end
    pred = zeros(N ≤ typemax(Int16) ? Int16 : Int32, K, N, N)
    incoming = [sortperm(σ[1:j-1, j]) for j ∈ 1:N]
    outgoing = [j .+ sortperm(σ[j, j+1:N]) for j ∈ 1:N]
    best = fill(Inf, K)
    last = ones(Int, K)
    best[1] = cost[1, N]
    for s ∈ 2:K
        next = fill(Inf, N, N)
        for j ∈ s:N-1
            in, out = incoming[j], outgoing[j]
            ptr, c, cᵢ = 1, Inf, 0
            for l ∈ out
                # Include all incoming edges with slope at most the outgoing
                while ptr ≤ length(in) && _leq(σ[in[ptr], j], σ[j, l])
                    i = in[ptr]
                    if cost[i, j] < c
                        c, cᵢ = cost[i, j], i
                    end
                    ptr += 1
                end
                if c < Inf
                    next[j, l] = op(c, p[j, l])
                    pred[s, j, l] = cᵢ
                end
            end
        end
        cost = next
        for i ∈ 1:N-1
            if cost[i, N] < best[s]
                best[s], last[s] = cost[i, N], i
            end
        end
    end
    return best, last, pred
end

function _convex_backtrack(pred, s, i, N)
    path = [i, N]
    for t ∈ s:-1:2
        pushfirst!(path, pred[t, path[1], path[2]])
    end
    return path
end

# Slope comparison with a small relative tolerance, so that collinear points
# are not rejected due to rounding
_leq(a, b) = a ≤ b + 1e-12 * max(abs(a), abs(b), 1.0)

# Function values closest to z (in metric) such that the points (x, z) are
# convex, i.e. with nondecreasing slopes between consecutive points
function _convexify_values(x, z, metric, optimizer)
    N = length(x)
    m = Model(optimizer)
    set_silent(m)
    @variable(m, v[1:N])
    for i ∈ 2:N-1
        @constraint(
            m,
            (v[i] - v[i-1]) / (x[i] - x[i-1]) ≤
            (v[i+1] - v[i]) / (x[i+1] - x[i])
        )
    end
    if metric == :l2
        @objective(m, Min, sum((v[i] - z[i])^2 for i ∈ 1:N))
    elseif metric == :l1
        @variable(m, t[1:N])
        @constraint(m, t .≥ v .- z)
        @constraint(m, t .≥ z .- v)
        @objective(m, Min, sum(t))
    elseif metric == :max
        @variable(m, t)
        @constraint(m, t .≥ v .- z)
        @constraint(m, t .≥ z .- v)
        @objective(m, Min, t)
    else
        error("Unrecognized/unsupported metric $(metric)")
    end
    optimize!(m)
    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end
    return value.(v)
end

# Deviation of a convex pwa function from the data in the given metric
function _ipol_error(pwa, x, z, metric)
    r = [abs(evaluate(pwa, (x[k],)) - z[k]) for k ∈ eachindex(x)]
    return metric == :max ? maximum(r) : metric == :l2 ? sum(abs2, r) : sum(r)
end

# A: Interpolate the data, then adjust the knot values to make it convex
function _ipol_knots(x, z, options)
    return _asconvex(_convexify1D(_interpolatepw(x, z, options), options))
end

# B: Make the data convex, then interpolate. Any interpolant through convex
# points is convex. Deviations are measured against the original data.
function _ipol_data(x, z, options)
    v = _convexify_values(x, z, options.metric, options.optimizer)
    p = _segment_penalties(x, v, z, options.metric)
    minimax = options.metric == :max
    knots, _ = _shortest_path_segments(p, options.planes; minimax)
    return _asconvex(PWAFunction1D(x[knots], v[knots]))
end

# C: The best convex interpolant through a subset of the data points, with at
# most options.planes segments. Returns nothing if no convex interpolant exists.
function _ipol_dp(x, z, options)
    N = length(x)
    σ = [i < j ? (z[j] - z[i]) / (x[j] - x[i]) : 0.0 for i ∈ 1:N, j ∈ 1:N]
    p = _segment_penalties(x, z, options.metric)
    minimax = options.metric == :max
    knots = _convex_shortest_path(σ, p, options.planes; minimax)
    isempty(knots) && return nothing
    return _asconvex(PWAFunction1D(x[knots], z[knots]))
end

# Fits with k = 1, ..., K planes for the given convexify method, computed from
# a single dynamic program (without the minimax tie-breaking used for the final
# fit). Used for selecting the number of planes.
function _ipol_candidates(x, z, options, K)
    method, metric = options.convexify, options.metric
    op = metric == :max ? max : +
    if method == :dp || (method == :auto && isnothing(options.optimizer))
        N = length(x)
        σ = [i < j ? (z[j] - z[i]) / (x[j] - x[i]) : 0.0 for i ∈ 1:N, j ∈ 1:N]
        p = _segment_penalties(x, z, metric)
        return [
            _asconvex(PWAFunction1D(x[q], z[q])) for
            q ∈ _convex_dp_all_paths(σ, p, K, op)
        ]
    elseif method == :data
        v = _convexify_values(x, z, metric, options.optimizer)
        p = _segment_penalties(x, v, z, metric)
        return [
            _asconvex(PWAFunction1D(x[q], v[q])) for q ∈ _dp_all_paths(p, K, op)
        ]
    elseif method == :knots
        p = _segment_penalties(x, z, metric)
        return [
            _asconvex(_convexify1D(x[q], z[q], options)) for
            q ∈ _dp_all_paths(p, K, op)
        ]
    elseif method == :auto
        dp = _ipol_candidates(x, z, _with(options; convexify = :dp), K)
        _isconvex(x, z) && return dp
        data = _ipol_candidates(x, z, _with(options; convexify = :data), K)
        return map(dp, data) do a, b
            ea, eb = (_ipol_error(f, x, z, metric) for f ∈ (a, b))
            return eb < ea ? b : a
        end
    end
    throw(ArgumentError("Unknown convexify option :$method"))
end

# Copy of the Interpol options with some fields changed
function _with(options::Interpol; kwargs...)
    fields = (f => getfield(options, f) for f ∈ fieldnames(Interpol))
    return Interpol(; fields..., kwargs...)
end

# Select the number of planes (at most K) according to options.selection.
# Returns the selected number and the score for each number of planes.
function _select_planes(x, z, options, K)
    N, metric, selection = length(x), options.metric, options.selection
    if selection == :cv
        # One standard error rule: the fewest planes with a cross validation
        # error within one standard error (over the folds) of the minimum
        E = _cv_errors(x, z, options, K)
        F = size(E, 1)
        scores = vec(mean(E, dims = 1))
        se = vec(std(E, dims = 1)) ./ sqrt(F)
        kᵐⁱⁿ = argmin(scores)
        return findfirst(≤(scores[kᵐⁱⁿ] + se[kᵐⁱⁿ]), scores), scores
    end
    fits = _ipol_candidates(x, z, options, K)
    e = [_ipol_error(f, x, z, metric) for f ∈ fits]
    if selection == :tolerance
        options.tolerance > 0 || throw(
            ArgumentError("selection = :tolerance requires tolerance > 0"),
        )
        # Error per data point, comparable to the tolerance
        scores = metric == :l1 ? e ./ N : metric == :l2 ? sqrt.(e ./ N) : e
        k = findfirst(≤(options.tolerance), scores)
        if isnothing(k)
            @warn "Tolerance $(options.tolerance) not reached with $K planes " *
                  "(error $(minimum(scores)))"
            k = K
        end
        return k, scores
    elseif selection ∈ (:aic, :bic)
        metric == :max && throw(
            ArgumentError("selection = :$selection requires metric :l1 or :l2"),
        )
        # -2 log-likelihood for Gaussian (:l2) or Laplace (:l1) errors, with
        # 2 parameters per plane (knot values and breakpoints)
        m2ll = if metric == :l2
            N .* log.(2π .* e ./ N) .+ N
        else
            2N .* log.(2 .* e ./ N) .+ 2N
        end
        nparams = [2length(f.planes) for f ∈ fits]
        penalty = selection == :aic ? 2.0 : log(N)
        scores = m2ll .+ penalty .* nparams
        return argmin(scores), scores
    end
    throw(ArgumentError("Unknown selection option :$selection"))
end

# Cross validation errors E[f, k] for fold f and k = 1, ..., K planes: the mean
# absolute deviation (:l1), mean squared deviation (:l2) or maximum absolute
# deviation (:max) on the held-out points. The points (sorted by x) are assigned
# to the folds in an interleaved manner, so that all folds cover the full range
# of the data.
function _cv_errors(x, z, options, K)
    N, metric = length(x), options.metric
    F = min(options.folds, N)
    F ≥ 2 || throw(ArgumentError("selection = :cv requires folds ≥ 2"))
    fold = [mod(i - 1, F) + 1 for i ∈ 1:N]
    E = fill(Inf, F, K)
    for f ∈ 1:F
        train, test = findall(≠(f), fold), findall(==(f), fold)
        length(train) ≥ 2 ||
            throw(ArgumentError("Too few points for $F-fold cross validation"))
        Kᶠ = min(K, length(train) - 1)
        fits = _ipol_candidates(x[train], z[train], options, Kᶠ)
        for k ∈ 1:Kᶠ
            r = [abs(evaluate(fits[k], (x[t],)) - z[t]) for t ∈ test]
            E[f, k] =
                metric == :max ? maximum(r) :
                metric == :l2 ? mean(abs2, r) : mean(r)
        end
    end
    return E
end

function _convex_linearization_ipol(x, z, options)
    @assert length(x) == length(z)
    x, z = _sort_unique(x, z)
    N = length(x)
    method = options.convexify
    if method ∈ (:data, :knots) && isnothing(options.optimizer)
        error("Interpol with convexify = :$method requires an optimizer")
    end
    if options.selection != :none
        # planes is the maximum number of planes to consider
        if options.planes < 1 || N < 2
            throw(ArgumentError("Interpol requires planes ≥ 1 and 2 points"))
        end
        K = min(options.planes, N - 1)
        k, scores = _select_planes(x, z, options, K)
        @debug "Interpol selected $k planes" options.selection scores
        options = _with(options; planes = k, selection = :none)
    elseif !(1 ≤ options.planes ≤ N - 1)
        throw(
            ArgumentError(
                "Interpol with planes = $(options.planes) requires at least " *
                "planes + 1 distinct points, got $N",
            ),
        )
    end
    return _ipol_fit(x, z, options)
end

function _ipol_fit(x, z, options)
    method = options.convexify
    if method == :auto
        return _ipol_auto(x, z, options)
    elseif method == :dp
        return _ipol_dp(x, z, options)
    elseif method == :data
        return _ipol_data(x, z, options)
    elseif method == :knots
        return _ipol_knots(x, z, options)
    end
    throw(ArgumentError("Unknown convexify option :$method"))
end

# Use the convex dynamic program if the data is convex (where all methods give
# the same result) or no optimizer is available. Otherwise also make the data
# convex before interpolating, and return the best of the two.
function _ipol_auto(x, z, options)
    pwa = _ipol_dp(x, z, options)
    if _isconvex(x, z) || isnothing(options.optimizer)
        return pwa
    end
    alt = _ipol_data(x, z, options)
    e, eₐ = (_ipol_error(f, x, z, options.metric) for f ∈ (pwa, alt))
    return eₐ < e ? alt : pwa
end

# Whether the points (sorted by x) have nondecreasing slopes
function _isconvex(x, z)
    s = [(z[i+1] - z[i]) / (x[i+1] - x[i]) for i ∈ 1:length(x)-1]
    return all(_leq(s[i], s[i+1]) for i ∈ 1:length(s)-1)
end

# Sort the points by x and merge points with identical x, using the average
# of their function values
function _sort_unique(x, z)
    xs = Float64[]
    zs = Float64[]
    n = Int[]
    for i ∈ sortperm(x)
        if !isempty(xs) && x[i] == xs[end]
            n[end] += 1
            zs[end] += (z[i] - zs[end]) / n[end]
        else
            push!(xs, x[i])
            push!(zs, z[i])
            push!(n, 1)
        end
    end
    return xs, zs
end
