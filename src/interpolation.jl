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

# One-dimensional data with distinct x in increasing order. Repeated
# measurements at the same x are kept as a group: z is the mean of the group,
# n the number of measurements and vals the measurements themselves, so that
# deviations can be measured against every measurement.
struct _Data1D
    x::Vector{Float64}
    z::Vector{Float64}
    n::Vector{Int}
    vals::Vector{Vector{Float64}}
end
_Data1D(x::AbstractVector, z::AbstractVector) = _sort_unique(x, z)
Base.length(d::_Data1D) = length(d.x)

# The total number of measurements
_total(d::_Data1D) = sum(d.n)

# The data at the points with the given indices
_subset(d::_Data1D, idx) = _Data1D(d.x[idx], d.z[idx], d.n[idx], d.vals[idx])

# Sort the points by x and group points with the same x. Points closer than
# a small fraction of the range of x are also grouped, as the slopes between
# them are not meaningful, and the group is placed at their mean.
function _sort_unique(x, z; rtol = 1e-6)
    xs = Float64[]
    vals = Vector{Float64}[]
    sums = Float64[]
    tol = rtol * (maximum(x) - minimum(x))
    for i ∈ sortperm(x)
        if !isempty(xs) && x[i] - xs[end] ≤ tol
            push!(vals[end], z[i])
            sums[end] += x[i]
        else
            push!(xs, x[i])
            push!(vals, [z[i]])
            push!(sums, x[i])
        end
    end
    return _Data1D(sums ./ length.(vals), mean.(vals), length.(vals), vals)
end

function _asconvex(pwa::PWAFunction1D)
    N = length(pwa.x)
    x = pwa.x
    z = pwa.z

    c = [(z[i+1] - z[i]) / (x[i+1] - x[i]) for i ∈ 1:(N-1)]
    d = [z[i] - c[i] * x[i] for i ∈ 1:(N-1)]

    return PWAFunc{Convex,1}(Plane.(c, d))
end

# Create a pwa interpolant to the data using exactly options.planes segments
# with end points at data points. The segments are chosen to minimize the
# total deviation in options.metric: the sum over the segments for :l1 and :l2,
# and the maximum over the segments for :max.
function _interpolatepw(d::_Data1D, options)
    p = _segment_penalties(d, options.metric)
    minimax = options.metric == :max
    knots, _ = _shortest_path_segments(p, options.planes; minimax)
    return PWAFunction1D(d.x[knots], d.z[knots])
end
function _interpolatepw(x::AbstractVector, z::AbstractVector, options)
    return _interpolatepw(_Data1D(x, z), options)
end

# Penalty p[i, j] (i < j) for approximating the measurements at the points
# k = i, ..., j-1 by the line segment between the knots (xᵢ, zᵏᵢ) and
# (xⱼ, zᵏⱼ): the sum of absolute deviations (:l1), the sum of squared
# deviations (:l2) or the maximum absolute deviation (:max). Each point is
# thereby assigned to exactly one segment, except the last point that is
# always a knot. Deviations are measured against every measurement, so
# repeated measurements at a point count with their number and spread.
function _segment_penalties(d::_Data1D, zᵏ, metric)
    if metric == :l2
        return _segment_penalties_l2(d, zᵏ)
    elseif metric ∉ (:l1, :max)
        error("Unrecognized/unsupported metric $(metric)")
    end
    x = d.x
    N = length(x)
    p = zeros(N, N)
    for i ∈ 1:N-1, j ∈ i+1:N
        c = (zᵏ[j] - zᵏ[i]) / (x[j] - x[i])
        acc = 0.0
        for k ∈ i:j-1
            f = zᵏ[i] + c * (x[k] - x[i])
            for zₖ ∈ d.vals[k]
                r = abs(f - zₖ)
                acc = metric == :l1 ? acc + r : max(acc, r)
            end
        end
        p[i, j] = acc
    end
    return p
end
_segment_penalties(d::_Data1D, metric) = _segment_penalties(d, d.z, metric)
function _segment_penalties(x::AbstractVector, zᵏ, zᵈ::AbstractVector, metric)
    return _segment_penalties(_Data1D(x, zᵈ), zᵏ, metric)
end
function _segment_penalties(x::AbstractVector, z::AbstractVector, metric)
    return _segment_penalties(_Data1D(x, z), metric)
end

# Sum of squared deviations for all segments in O(N²) using prefix sums over
# the groups of measurements. The data is centered first to reduce
# cancellation errors.
function _segment_penalties_l2(d::_Data1D, zᵏ)
    N = length(d)
    n = Float64.(d.n)
    x̄ = sum(n .* d.x) / sum(n)
    z̄ = sum(sum.(d.vals)) / sum(n)
    x = d.x .- x̄
    zᵏ = zᵏ .- z̄
    # Sums of the centered measurements and their squares per group
    sz = [sum(v .- z̄) for v ∈ d.vals]
    szz = [sum(abs2, v .- z̄) for v ∈ d.vals]
    Sn = [0.0; cumsum(n)]
    Sx = [0.0; cumsum(n .* x)]
    Sxx = [0.0; cumsum(n .* x .^ 2)]
    Sz = [0.0; cumsum(sz)]
    Szz = [0.0; cumsum(szz)]
    Sxz = [0.0; cumsum(x .* sz)]
    p = zeros(N, N)
    for i ∈ 1:N-1, j ∈ i+1:N
        # Deviations at the points k = i, ..., j-1 are α + c xₖ - z
        c = (zᵏ[j] - zᵏ[i]) / (x[j] - x[i])
        α = zᵏ[i] - c * x[i]
        m, sx, sz₁ = Sn[j] - Sn[i], Sx[j] - Sx[i], Sz[j] - Sz[i]
        sxx, szz₁, sxz = Sxx[j] - Sxx[i], Szz[j] - Szz[i], Sxz[j] - Sxz[i]
        r² = m * α^2 + c^2 * sxx + szz₁ + 2α * c * sx - 2α * sz₁ - 2c * sxz
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
    # The costs of the previous and the current number of edges alternate
    # between two buffers
    next = similar(cost)
    for s ∈ 2:K
        fill!(next, Inf)
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
        cost, next = next, cost
        for i ∈ 1:N-1
            if cost[i, N] < best[s]
                best[s], last[s] = cost[i, N], i
            end
        end
    end
    return best, last, pred
end

# The memory needed by the convex dynamic program for N points and K planes
# (penalties, slopes, two cost buffers and the predecessors), which grows as
# N² and limits the size of the data that can be handled
function _check_size(N, K; limit = 4 * 2^30)
    bytes = (4 * 8 + 2K) * N^2
    bytes ≤ limit || throw(
        ArgumentError(
            "Interpol with $N distinct points and up to $K planes needs " *
            "about $(round(bytes / 2^30; digits = 1)) GiB of memory. " *
            "Reduce the number of points (e.g. by sampling) or planes.",
        ),
    )
    return
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

# Function values v at the points of d closest to the measurements (in
# metric) such that the points (x, v) are convex, i.e. with nondecreasing
# slopes between consecutive points
function _convexify_values(d::_Data1D, metric, optimizer)
    x, N = d.x, length(d)
    m = Model(optimizer)
    set_silent(m)
    @variable(m, v[1:N])
    _slope_constraints!(m, x, v)
    _deviation_objective!(m, v, d, metric)
    optimize!(m)
    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end
    return value.(v)
end
function _convexify_values(x::AbstractVector, z::AbstractVector, metric, opt)
    return _convexify_values(_Data1D(x, z), metric, opt)
end

# Constraints that the values v (variables or expressions) at the increasing
# points x are convex, i.e. have nondecreasing slopes between consecutive
# points. The slope condition is multiplied by the lengths of the two intervals
# and scaled by the longer one, so that the coefficients lie in (0, 1] instead
# of growing as the inverse of the interval lengths.
function _slope_constraints!(m, x, v)
    for i ∈ 2:length(x)-1
        h₋, h₊ = x[i] - x[i-1], x[i+1] - x[i]
        h = max(h₋, h₊)
        @constraint(
            m,
            (v[i] - v[i-1]) * (h₊ / h) ≤ (v[i+1] - v[i]) * (h₋ / h)
        )
    end
    return
end

# Sets the objective of the model to the deviation of the function values f
# (expressions or variables, one per point of d) from all the measurements, in
# the given metric
function _deviation_objective!(m, f, d::_Data1D, metric)
    N = length(d)
    if metric == :l2
        # Σ (f - z)² over the measurements, up to a constant
        @objective(
            m,
            Min,
            sum(d.n[g] * f[g]^2 - 2 * sum(d.vals[g]) * f[g] for g ∈ 1:N)
        )
    elseif metric == :l1
        @variable(m, t[g = 1:N, k = 1:d.n[g]])
        @constraint(m, [g = 1:N, k = 1:d.n[g]], t[g, k] ≥ f[g] - d.vals[g][k])
        @constraint(m, [g = 1:N, k = 1:d.n[g]], t[g, k] ≥ d.vals[g][k] - f[g])
        @objective(m, Min, sum(t))
    elseif metric == :max
        @variable(m, t)
        @constraint(m, [g = 1:N], t ≥ f[g] - minimum(d.vals[g]))
        @constraint(m, [g = 1:N], t ≥ maximum(d.vals[g]) - f[g])
        @objective(m, Min, t)
    else
        error("Unrecognized/unsupported metric $(metric)")
    end
    return
end

# Absolute deviations of the function values f at the points of d from all
# the measurements, and their aggregation in a metric
function _deviations(d::_Data1D, f)
    return [abs(f[g] - zₖ) for g ∈ eachindex(f) for zₖ ∈ d.vals[g]]
end
function _aggregate(r, metric)
    return metric == :max ? maximum(r) : metric == :l2 ? sum(abs2, r) : sum(r)
end

# Deviation of a convex pwa function from the data in the given metric
function _ipol_error(pwa, d::_Data1D, metric)
    f = [evaluate(pwa, (xg,)) for xg ∈ d.x]
    return _aggregate(_deviations(d, f), metric)
end
function _ipol_error(pwa, x::AbstractVector, z::AbstractVector, metric)
    return _ipol_error(pwa, _Data1D(x, z), metric)
end

# The two methods below return the knots of a convex interpolant.

# B: Make the data convex, then interpolate. Any interpolant through convex
# points is convex. Deviations are measured against the original data.
function _ipol_data(d::_Data1D, options)
    v = _convexify_values(d, options.metric, options.optimizer)
    p = _segment_penalties(d, v, options.metric)
    minimax = options.metric == :max
    knots, _ = _shortest_path_segments(p, options.planes; minimax)
    return PWAFunction1D(d.x[knots], v[knots])
end

# C: The best convex interpolant through a subset of the data points, with at
# most options.planes segments. A convex interpolant always exists, since the
# line segment between the end points is convex.
function _ipol_dp(d::_Data1D, options)
    σ = _slopes(d)
    p = _segment_penalties(d, options.metric)
    minimax = options.metric == :max
    knots = _convex_shortest_path(σ, p, options.planes; minimax)
    return PWAFunction1D(d.x[knots], d.z[knots])
end

# Slopes of the line segments between the points of d
function _slopes(d::_Data1D)
    x, z, N = d.x, d.z, length(d)
    return [i < j ? (z[j] - z[i]) / (x[j] - x[i]) : 0.0 for i ∈ 1:N, j ∈ 1:N]
end

# Knots of the fits with k = 1, ..., K planes for the given convexify method,
# computed from a single dynamic program (without the minimax tie-breaking used
# for the final fit). Used for selecting the number of planes.
function _ipol_candidates(d::_Data1D, options, K)
    method, metric = options.convexify, options.metric
    op = metric == :max ? max : +
    if method == :dp || (method == :auto && isnothing(options.optimizer))
        p = _segment_penalties(d, metric)
        paths = _convex_dp_all_paths(_slopes(d), p, K, op)
        return [PWAFunction1D(d.x[q], d.z[q]) for q ∈ paths]
    elseif method == :data
        v = _convexify_values(d, metric, options.optimizer)
        p = _segment_penalties(d, v, metric)
        return [PWAFunction1D(d.x[q], v[q]) for q ∈ _dp_all_paths(p, K, op)]
    elseif method == :auto
        dp = _ipol_candidates(d, _with(options; convexify = :dp), K)
        _isconvex(d.x, d.z) && return dp
        data = _ipol_candidates(d, _with(options; convexify = :data), K)
        return map(dp, data) do a, b
            ea, eb = (_knots_error(k, d, metric) for k ∈ (a, b))
            return eb < ea ? b : a
        end
    end
    throw(ArgumentError("Unknown convexify option :$method"))
end

# Candidate fits used for selection. If the knots are to be adjusted, the
# function values of the candidates are refitted, as a cheaper proxy also for
# moving the breakpoints.
function _selection_fits(d::_Data1D, options, K)
    fits = _ipol_candidates(d, options, K)
    options.adjust == :none && return fits
    return [_adjust_values(d, k, options) for k ∈ fits]
end

# Copy of the Interpol options with some fields changed
function _with(options::Interpol; kwargs...)
    fields = (f => getfield(options, f) for f ∈ fieldnames(Interpol))
    return Interpol(; fields..., kwargs...)
end

# Select the number of planes (at most K) according to options.selection.
# Returns the selected number, the score for each number of planes, and the
# candidate fits (nothing for cross validation, where the fits are to parts of
# the data).
function _select_planes(d::_Data1D, options, K; scale = 1.0)
    N, metric, selection = _total(d), options.metric, options.selection
    fits = nothing
    k, scores = if selection == :cv
        # One standard error rule: the fewest planes with a cross validation
        # error within one standard error (over the folds) of the minimum
        E = _cv_errors(d, options, K)
        F = size(E, 1)
        scores = vec(mean(E, dims = 1))
        se = vec(std(E, dims = 1)) ./ sqrt(F)
        kᵐⁱⁿ = argmin(scores)
        findfirst(≤(scores[kᵐⁱⁿ] + se[kᵐⁱⁿ]), scores), scores
    elseif selection == :tolerance
        options.tolerance > 0 || throw(
            ArgumentError("selection = :tolerance requires tolerance > 0"),
        )
        fits = _selection_fits(d, options, K)
        e = [_knots_error(k, d, metric) for k ∈ fits]
        # Error per measurement, comparable to the tolerance
        scores = metric == :l1 ? e ./ N : metric == :l2 ? sqrt.(e ./ N) : e
        k = findfirst(≤(options.tolerance), scores)
        if isnothing(k)
            @warn "Tolerance $(options.tolerance * scale) not reached with " *
                  "$K planes (error $(minimum(scores) * scale))"
            k = K
        end
        k, scores
    elseif selection ∈ (:aic, :bic)
        metric == :max && throw(
            ArgumentError("selection = :$selection requires metric :l1 or :l2"),
        )
        fits = _selection_fits(d, options, K)
        e = [_knots_error(k, d, metric) for k ∈ fits]
        # Treat errors at the level of rounding errors as exact fits, so that
        # these are compared equal
        scale = max(maximum(abs, d.z), floatmin())
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
        argmin(scores), scores
    else
        throw(ArgumentError("Unknown selection option :$selection"))
    end
    if selection != :tolerance && k == K && K == options.planes
        @warn "The maximum number of planes ($K) was selected by " *
              ":$selection, consider increasing `planes`"
    end
    return k, scores, fits
end

# Cross validation errors E[f, k] for fold f and k = 1, ..., K planes: the mean
# absolute deviation (:l1), mean squared deviation (:l2) or maximum absolute
# deviation (:max) on the held-out measurements. The points (sorted by x) are
# assigned to the folds in an interleaved manner, so that all folds cover the
# full range of the data, and all measurements at a point are held out
# together.
function _cv_errors(d::_Data1D, options, K)
    N, metric = length(d), options.metric
    F = min(options.folds, N)
    F ≥ 2 || throw(ArgumentError("selection = :cv requires folds ≥ 2"))
    fold = [mod(i - 1, F) + 1 for i ∈ 1:N]
    E = fill(Inf, F, K)
    for f ∈ 1:F
        train, test = findall(≠(f), fold), findall(==(f), fold)
        length(train) ≥ 2 ||
            throw(ArgumentError("Too few points for $F-fold cross validation"))
        Kᶠ = min(K, length(train) - 1)
        fits = _selection_fits(_subset(d, train), options, Kᶠ)
        dtest = _subset(d, test)
        for k ∈ 1:Kᶠ
            r = _deviations(dtest, _interpolate(fits[k].x, fits[k].z, dtest.x))
            E[f, k] =
                metric == :max ? maximum(r) :
                metric == :l2 ? mean(abs2, r) : mean(r)
        end
    end
    return E
end

function _convex_linearization_ipol(x, z, options)
    @assert length(x) == length(z)
    options = _validate(options)
    # Work with points and values normalized to the order of one, so that the
    # tolerances of the algorithms are independent of the scale of the data
    x̄, sx = _shift_scale(x)
    z̄, sz = _shift_scale(z)
    d = _sort_unique((x .- x̄) ./ sx, (z .- z̄) ./ sz)
    options = _with(options; tolerance = options.tolerance / sz)
    knots = _fit_knots(d, options, sz)
    return _denormalize_1d(_asconvex(knots), x̄, sx, z̄, sz)
end

# Shift and scale mapping the values to [-1, 1]
function _shift_scale(v)
    lo, hi = extrema(v)
    return (lo + hi) / 2, hi > lo ? (hi - lo) / 2 : 1.0
end

# The planes fitted to normalized points x̃ = (x - x̄) / sx and values
# z̃ = (z - z̄) / sz, in the original coordinates
function _denormalize_1d(pwa::PWAFunc{Convex,1}, x̄, sx, z̄, sz)
    planes = map(pwa.planes) do q
        a = q.α[1] * sz / sx
        return Plane((a,), sz * q.β - a * x̄ + z̄)
    end
    return PWAFunc{Convex,1}(planes)
end

# The knots of the fit to the normalized data d, with the number of planes
# selected if requested. The scale of the values is used to report the
# tolerance in the original units.
function _fit_knots(d::_Data1D, options, sz)
    N = length(d)
    method = options.convexify
    if method == :data && isnothing(options.optimizer)
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
        N ≥ 2 || throw(ArgumentError("Interpol requires at least 2 points"))
        K = min(options.planes, N - 1)
        _check_size(N, K)
        k, scores, fits = _select_planes(d, options, K; scale = sz)
        @debug "Interpol selected $k planes" options.selection scores
        options = _with(options; planes = k, selection = :none)
        # The candidate fit is the final fit, except for the minimax tie
        # breaking of :max, and the search over the breakpoints
        if !isnothing(fits) && options.metric != :max
            knots = fits[k]
            if options.adjust == :breakpoints
                knots = _adjust_breakpoints(d, knots, options)
            end
            return knots
        end
    elseif !(1 ≤ options.planes ≤ N - 1)
        throw(
            ArgumentError(
                "Interpol with planes = $(options.planes) requires at least " *
                "planes + 1 distinct points, got $N",
            ),
        )
    else
        _check_size(N, options.planes)
    end
    return _ipol_fit(d, options)
end

# Check the options, and return them with :rms replaced by :l2
function _validate(options::Interpol)
    valid(name, value, values) =
        value ∈ values || throw(
            ArgumentError(
                "Unsupported $name = :$value for Interpol, use one of " *
                join(":" .* string.(values), ", "),
            ),
        )
    valid("metric", options.metric, (:l1, :l2, :rms, :max))
    valid("convexify", options.convexify, (:auto, :dp, :data))
    valid("selection", options.selection, (:none, :tolerance, :aic, :bic, :cv))
    valid("adjust", options.adjust, (:none, :values, :breakpoints))
    options.planes ≥ 1 ||
        throw(ArgumentError("Interpol requires planes ≥ 1"))
    options.selection == :cv && options.folds < 2 &&
        throw(ArgumentError("selection = :cv requires folds ≥ 2"))
    options.selection == :tolerance && options.tolerance ≤ 0 &&
        throw(ArgumentError("selection = :tolerance requires tolerance > 0"))
    return options.metric == :rms ? _with(options; metric = :l2) : options
end

function _ipol_fit(d::_Data1D, options)
    method = options.convexify
    knots = if method == :auto
        _ipol_auto(d, options)
    elseif method == :dp
        _ipol_dp(d, options)
    elseif method == :data
        _ipol_data(d, options)
    else
        throw(ArgumentError("Unknown convexify option :$method"))
    end
    if options.adjust == :values
        knots = _adjust_values(d, knots, options)
    elseif options.adjust == :breakpoints
        knots = _adjust_breakpoints(d, knots, options)
    elseif options.adjust != :none
        throw(ArgumentError("Unknown adjust option :$(options.adjust)"))
    end
    return knots
end

# Use the convex dynamic program if the data is convex (where all methods give
# the same result) or no optimizer is available. Otherwise also make the data
# convex before interpolating, and return the best of the two.
function _ipol_auto(d::_Data1D, options)
    knots = _ipol_dp(d, options)
    if _isconvex(d.x, d.z) || isnothing(options.optimizer)
        return knots
    end
    alt = _ipol_data(d, options)
    e, eₐ = (_knots_error(k, d, options.metric) for k ∈ (knots, alt))
    return eₐ < e ? alt : knots
end

# Whether the points (sorted by x) have nondecreasing slopes
function _isconvex(x, z)
    s = [(z[i+1] - z[i]) / (x[i+1] - x[i]) for i ∈ 1:length(x)-1]
    return all(_leq(s[i], s[i+1]) for i ∈ 1:length(s)-1)
end

# Deviation of the piecewise linear function through the knots from the data
function _knots_error(knots::PWAFunction1D, d::_Data1D, metric)
    f = _interpolate(knots.x, knots.z, d.x)
    return _aggregate(_deviations(d, f), metric)
end
function _knots_error(knots::PWAFunction1D, x::AbstractVector, z, metric)
    return _knots_error(knots, _Data1D(x, z), metric)
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
function _adjust_values(d::_Data1D, knots, options)
    v = _fit_values(d, knots.x, options.metric, options.optimizer)
    adjusted = PWAFunction1D(knots.x, v)
    e, eₐ = (_knots_error(k, d, options.metric) for k ∈ (knots, adjusted))
    return eₐ ≤ e ? adjusted : knots
end

# Also move the interior breakpoints, by repeated sweeps where each breakpoint
# is moved by a golden section search between its neighbours, refitting the
# function values for each candidate position. A move is only accepted if it
# reduces the deviation. Finds a local optimum.
function _adjust_breakpoints(
    d::_Data1D,
    knots,
    options;
    sweeps = 10,
    steps = 25,
)
    metric, optimizer = options.metric, options.optimizer
    t = copy(knots.x)
    fit(t) = PWAFunction1D(t, _fit_values(d, t, metric, optimizer))
    err(t) = _knots_error(fit(t), d, metric)
    # Start from the best of the given knots and refitted values
    refit = fit(t)
    e, eᵣ = (_knots_error(k, d, metric) for k ∈ (knots, refit))
    best, e = eᵣ < e ? (refit, eᵣ) : (knots, e)
    φ = (sqrt(5) - 1) / 2
    for _ ∈ 1:sweeps
        e₀ = e
        for i ∈ 2:length(t)-1
            # Keep a small margin to the neighbouring breakpoints
            δ = 1e-3 * (t[i+1] - t[i-1])
            a, b = t[i-1] + δ, t[i+1] - δ
            errᵢ(τ) = err([t[1:i-1]; τ; t[i+1:end]])
            c, d₁ = b - φ * (b - a), a + φ * (b - a)
            ec, ed = errᵢ(c), errᵢ(d₁)
            for _ ∈ 1:steps
                if ec < ed
                    b, d₁, ed = d₁, c, ec
                    c = b - φ * (b - a)
                    ec = errᵢ(c)
                else
                    a, c, ec = c, d₁, ed
                    d₁ = a + φ * (b - a)
                    ed = errᵢ(d₁)
                end
            end
            τ, eτ = ec < ed ? (c, ec) : (d₁, ed)
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
# of the convex piecewise linear function with minimal deviation from the
# measurements
function _fit_values(d::_Data1D, t, metric, optimizer)
    if metric == :l2
        return _fit_values_l2(_expand(d)..., t)
    end
    isnothing(optimizer) && error(
        "Adjusting the interpolant with metric :$metric requires an optimizer",
    )
    n, N = length(t), length(d)
    m = Model(optimizer)
    set_silent(m)
    @variable(m, v[1:n])
    _slope_constraints!(m, t, v)
    # The value at a point is a convex combination of the values at the
    # surrounding breakpoints
    seg = [clamp(searchsortedlast(t, d.x[g]), 1, n - 1) for g ∈ 1:N]
    λ = [(d.x[g] - t[seg[g]]) / (t[seg[g]+1] - t[seg[g]]) for g ∈ 1:N]
    f = [(1 - λ[g]) * v[seg[g]] + λ[g] * v[seg[g]+1] for g ∈ 1:N]
    _deviation_objective!(m, f, d, metric)
    optimize!(m)
    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end
    return value.(v)
end
function _fit_values(x::AbstractVector, z::AbstractVector, t, metric, opt)
    return _fit_values(_Data1D(x, z), t, metric, opt)
end

# All measurements as vectors of points and values
function _expand(d::_Data1D)
    x = reduce(vcat, fill(d.x[g], d.n[g]) for g ∈ eachindex(d.x))
    return x, reduce(vcat, d.vals)
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
