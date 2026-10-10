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

# The three methods below return the knots of a convex interpolant.

# A: Interpolate the data, then adjust the knot values to make it convex
function _ipol_knots(x, z, options)
    return _convexify1D(_interpolatepw(x, z, options), options)
end

# B: Make the data convex, then interpolate. Any interpolant through convex
# points is convex. Deviations are measured against the original data.
function _ipol_data(x, z, options)
    v = _convexify_values(x, z, options.metric, options.optimizer)
    p = _segment_penalties(x, v, z, options.metric)
    minimax = options.metric == :max
    knots, _ = _shortest_path_segments(p, options.planes; minimax)
    return PWAFunction1D(x[knots], v[knots])
end

# C: The best convex interpolant through a subset of the data points, with at
# most options.planes segments. A convex interpolant always exists, since the
# line segment between the end points is convex.
function _ipol_dp(x, z, options)
    N = length(x)
    σ = [i < j ? (z[j] - z[i]) / (x[j] - x[i]) : 0.0 for i ∈ 1:N, j ∈ 1:N]
    p = _segment_penalties(x, z, options.metric)
    minimax = options.metric == :max
    knots = _convex_shortest_path(σ, p, options.planes; minimax)
    return PWAFunction1D(x[knots], z[knots])
end

# Knots of the fits with k = 1, ..., K planes for the given convexify method,
# computed from a single dynamic program (without the minimax tie-breaking used
# for the final fit). Used for selecting the number of planes.
function _ipol_candidates(x, z, options, K)
    method, metric = options.convexify, options.metric
    op = metric == :max ? max : +
    if method == :dp || (method == :auto && isnothing(options.optimizer))
        N = length(x)
        σ = [i < j ? (z[j] - z[i]) / (x[j] - x[i]) : 0.0 for i ∈ 1:N, j ∈ 1:N]
        p = _segment_penalties(x, z, metric)
        paths = _convex_dp_all_paths(σ, p, K, op)
        return [PWAFunction1D(x[q], z[q]) for q ∈ paths]
    elseif method == :data
        v = _convexify_values(x, z, metric, options.optimizer)
        p = _segment_penalties(x, v, z, metric)
        return [PWAFunction1D(x[q], v[q]) for q ∈ _dp_all_paths(p, K, op)]
    elseif method == :knots
        p = _segment_penalties(x, z, metric)
        paths = _dp_all_paths(p, K, op)
        return [_convexify1D(x[q], z[q], options) for q ∈ paths]
    elseif method == :auto
        dp = _ipol_candidates(x, z, _with(options; convexify = :dp), K)
        _isconvex(x, z) && return dp
        data = _ipol_candidates(x, z, _with(options; convexify = :data), K)
        return map(dp, data) do a, b
            ea, eb = (_knots_error(k, x, z, metric) for k ∈ (a, b))
            return eb < ea ? b : a
        end
    end
    throw(ArgumentError("Unknown convexify option :$method"))
end

# Candidate fits used for selection. If the knots are to be adjusted, the
# function values of the candidates are refitted, as a cheaper proxy also for
# moving the breakpoints.
function _selection_fits(x, z, options, K)
    fits = _ipol_candidates(x, z, options, K)
    options.adjust == :none && return fits
    return [_adjust_values(x, z, k, options) for k ∈ fits]
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
    fits = _selection_fits(x, z, options, K)
    e = [_knots_error(k, x, z, metric) for k ∈ fits]
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
        # Treat errors at the level of rounding errors as exact fits, so that
        # these are compared equal
        scale = max(maximum(abs, z), floatmin())
        e = max.(e, metric == :l2 ? N * (1e-12 * scale)^2 : N * 1e-12 * scale)
        # -2 log-likelihood for Gaussian (:l2) or Laplace (:l1) errors, with
        # 2 parameters per plane (knot values and breakpoints)
        m2ll = if metric == :l2
            N .* log.(2π .* e ./ N) .+ N
        else
            2N .* log.(2 .* e ./ N) .+ 2N
        end
        nparams = [2(length(k.x) - 1) for k ∈ fits]
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
        fits = _selection_fits(x[train], z[train], options, Kᶠ)
        for k ∈ 1:Kᶠ
            r = abs.(_interpolate(fits[k].x, fits[k].z, x[test]) .- z[test])
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
    if options.adjust != :none && options.metric != :l2 &&
       isnothing(options.optimizer)
        error(
            "Interpol with adjust = :$(options.adjust) and metric = " *
            ":$(options.metric) requires an optimizer",
        )
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
    knots = if method == :auto
        _ipol_auto(x, z, options)
    elseif method == :dp
        _ipol_dp(x, z, options)
    elseif method == :data
        _ipol_data(x, z, options)
    elseif method == :knots
        _ipol_knots(x, z, options)
    else
        throw(ArgumentError("Unknown convexify option :$method"))
    end
    if options.adjust == :values
        knots = _adjust_values(x, z, knots, options)
    elseif options.adjust == :breakpoints
        knots = _adjust_breakpoints(x, z, knots, options)
    elseif options.adjust != :none
        throw(ArgumentError("Unknown adjust option :$(options.adjust)"))
    end
    return _asconvex(knots)
end

# Use the convex dynamic program if the data is convex (where all methods give
# the same result) or no optimizer is available. Otherwise also make the data
# convex before interpolating, and return the best of the two.
function _ipol_auto(x, z, options)
    knots = _ipol_dp(x, z, options)
    if _isconvex(x, z) || isnothing(options.optimizer)
        return knots
    end
    alt = _ipol_data(x, z, options)
    e, eₐ = (_knots_error(k, x, z, options.metric) for k ∈ (knots, alt))
    return eₐ < e ? alt : knots
end

# Whether the points (sorted by x) have nondecreasing slopes
function _isconvex(x, z)
    s = [(z[i+1] - z[i]) / (x[i+1] - x[i]) for i ∈ 1:length(x)-1]
    return all(_leq(s[i], s[i+1]) for i ∈ 1:length(s)-1)
end

# Deviation of the piecewise linear function through the knots from the data
function _knots_error(knots::PWAFunction1D, x, z, metric)
    r = abs.(_interpolate(knots.x, knots.z, x) .- z)
    return metric == :max ? maximum(r) : metric == :l2 ? sum(abs2, r) : sum(r)
end

# Linear interpolation through (t, v) at the points x within [t[1], t[end]]
function _interpolate(t, v, x)
    return map(x) do xₖ
        j = clamp(searchsortedlast(t, xₖ), 1, length(t) - 1)
        λ = (xₖ - t[j]) / (t[j+1] - t[j])
        return (1 - λ) * v[j] + λ * v[j+1]
    end
end

# Keep the breakpoints, and refit the function values at the knots to minimize
# the deviation from the data, subject to convexity. As the current values are
# feasible, the deviation does not increase.
function _adjust_values(x, z, knots, options)
    v = _fit_values(x, z, knots.x, options.metric, options.optimizer)
    adjusted = PWAFunction1D(knots.x, v)
    e, eₐ = (_knots_error(k, x, z, options.metric) for k ∈ (knots, adjusted))
    return eₐ ≤ e ? adjusted : knots
end

# Also move the interior breakpoints, by repeated sweeps where each breakpoint
# is moved by a golden section search between its neighbours, refitting the
# function values for each candidate position. A move is only accepted if it
# reduces the deviation. Finds a local optimum.
function _adjust_breakpoints(x, z, knots, options; sweeps = 10, steps = 25)
    metric, optimizer = options.metric, options.optimizer
    t = copy(knots.x)
    fit(t) = PWAFunction1D(t, _fit_values(x, z, t, metric, optimizer))
    err(t) = _knots_error(fit(t), x, z, metric)
    # Start from the best of the given knots and refitted values
    refit = fit(t)
    e, eᵣ = (_knots_error(k, x, z, metric) for k ∈ (knots, refit))
    best, e = eᵣ < e ? (refit, eᵣ) : (knots, e)
    φ = (sqrt(5) - 1) / 2
    for _ ∈ 1:sweeps
        e₀ = e
        for i ∈ 2:length(t)-1
            # Keep a small margin to the neighbouring breakpoints
            δ = 1e-3 * (t[i+1] - t[i-1])
            a, b = t[i-1] + δ, t[i+1] - δ
            errᵢ(τ) = err([t[1:i-1]; τ; t[i+1:end]])
            c, d = b - φ * (b - a), a + φ * (b - a)
            ec, ed = errᵢ(c), errᵢ(d)
            for _ ∈ 1:steps
                if ec < ed
                    b, d, ed = d, c, ec
                    c = b - φ * (b - a)
                    ec = errᵢ(c)
                else
                    a, c, ec = c, d, ed
                    d = a + φ * (b - a)
                    ed = errᵢ(d)
                end
            end
            τ, eτ = ec < ed ? (c, ec) : (d, ed)
            if eτ < e
                t[i], e = τ, eτ
                best = fit(t)
            end
        end
        e ≥ e₀ - 1e-9 * max(e₀, 1e-12) && break
    end
    return best
end

# Function values v at the breakpoints t (with t[1] = x[1] and t[end] = x[end])
# of the convex piecewise linear function with minimal deviation from the data
function _fit_values(x, z, t, metric, optimizer)
    if metric == :l2
        return _fit_values_l2(x, z, t)
    end
    isnothing(optimizer) && error(
        "Adjusting the interpolant with metric :$metric requires an optimizer",
    )
    n, N = length(t), length(x)
    m = Model(optimizer)
    set_silent(m)
    @variable(m, v[1:n])
    for i ∈ 2:n-1
        @constraint(
            m,
            (v[i] - v[i-1]) / (t[i] - t[i-1]) ≤
            (v[i+1] - v[i]) / (t[i+1] - t[i])
        )
    end
    # The value at x[k] is a convex combination of the values at the
    # surrounding breakpoints
    seg = [clamp(searchsortedlast(t, x[k]), 1, n - 1) for k ∈ 1:N]
    λ = [(x[k] - t[seg[k]]) / (t[seg[k]+1] - t[seg[k]]) for k ∈ 1:N]
    f = [(1 - λ[k]) * v[seg[k]] + λ[k] * v[seg[k]+1] for k ∈ 1:N]
    if metric == :l1
        @variable(m, r[1:N])
        @constraint(m, r .≥ f .- z)
        @constraint(m, r .≥ z .- f)
        @objective(m, Min, sum(r))
    elseif metric == :max
        @variable(m, r)
        @constraint(m, r .≥ f .- z)
        @constraint(m, r .≥ z .- f)
        @objective(m, Min, r)
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

# Least squares fit of a convex piecewise linear function with breakpoints t.
# Written as a + b x + Σᵢ cᵢ max(0, x - tᵢ) over the interior breakpoints, the
# function is convex if and only if c ≥ 0. With (a, b) projected out, this is a
# nonnegative least squares problem.
function _fit_values_l2(x, z, t)
    A = [ones(length(x)) x]
    H = [max(0.0, xₖ - tᵢ) for xₖ ∈ x, tᵢ ∈ t[2:end-1]]
    # Project out the affine part
    Hₚ = H .- A * (A \ H)
    zₚ = z .- A * (A \ z)
    c = isempty(t[2:end-1]) ? Float64[] : _nnls(Hₚ, zₚ)
    ab = A \ (z .- H * c)
    v = ab[1] .+ ab[2] .* t
    for i ∈ eachindex(c)
        v .+= c[i] .* max.(0.0, t .- t[i+1])
    end
    return v
end

# Nonnegative least squares min ‖Ax - b‖ s.t. x ≥ 0 (Lawson & Hanson, 1974)
function _nnls(A, b; maxiter = 3size(A, 2) + 10)
    n = size(A, 2)
    tol = 10 * eps() * norm(A, 1) * max(size(A)...)
    x = zeros(n)
    P = falses(n)
    w = A' * (b - A * x)
    for _ ∈ 1:maxiter
        candidates = findall(j -> !P[j] && w[j] > tol, 1:n)
        isempty(candidates) && break
        P[candidates[argmax(w[candidates])]] = true
        while true
            s = zeros(n)
            s[P] = A[:, P] \ b
            if all(s[P] .> tol)
                x = s
                break
            end
            # Step towards s until a variable hits zero, and remove it
            J = findall(j -> P[j] && s[j] ≤ tol, 1:n)
            α = minimum(x[j] / (x[j] - s[j]) for j ∈ J)
            x .+= α .* (s .- x)
            P .&= x .> tol
            x[.!P] .= 0.0
        end
        w = A' * (b - A * x)
    end
    return x
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
