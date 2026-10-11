#=
This file holds an implementation of a progressive fitting heuristic
inspired by the algorithm in the paper "A linear programming approach to
difference-of-convex piecewise linear approximation" by Kazda and Li (2024).

The current implementation requires that the function evaluations are
samples of a convex function, i.e. it is not robust with regard to noisy
and non-convex data.

The main algorithm consists of the following steps:
1. While error < tolerance
    a. Introduce an extra plane in the pwa approximation
    b. Find local optimum for allocating points to plane (using linear programming)
=#

"""
    enforce_curvature(f::FunctionEvaluations, curvature::Curvature, optimizer, metric = :l1)

Create a slightly perturbed version of the function values
to ensure that the data points can be interpolated by a convex/concave
piecewise affine function.

Enforcing the curvature is performed using a linear optimization problem
that adjust the function values and tries to minimize the total deviation.
The total deviation can be measured in different metrics specified
by the `metric` parameter: `:l1` (sum of absolute changes), `:l2` (sum of
squared changes, solved as a quadratic program) or `:max` (largest change).
This function can be useful as a pre-step for running the full-order and progressive
fitting heuristics that require data points that are convex/concave in this sense.
"""
function enforce_curvature(
    f::FunctionEvaluations,
    c::Convex,
    optimizer,
    metric = :l1,
)
    return convexify(f, optimizer, metric)
end

function enforce_curvature(
    f::FunctionEvaluations,
    c::Concave,
    optimizer,
    metric = :l1,
)
    fneg = FunctionEvaluations(f.points, -f.values)
    fconvex = convexify(fneg, optimizer, metric, "concave")
    return FunctionEvaluations(f.points, -fconvex.values)
end

function convexify(
    f::FunctionEvaluations{D},
    optimizer,
    metric = :l1,
    curv_string = "convex",
) where {D}
    if metric ∉ (:l1, :l2, :rms, :max)
        throw(
            ArgumentError(
                "Unsupported metric :$metric for enforcing curvature, " *
                "use :l1, :l2 or :max",
            ),
        )
    end
    if D == 1 && allunique(first.(f.points))
        return _convexify_1d(f, optimizer, metric, curv_string)
    end
    K = length(f)
    X = [pt[i] for i ∈ 1:D, pt ∈ f.points]

    m = Model(optimizer)
    set_silent(m)
    @variable(m, s⁺[1:K] ≥ 0)
    @variable(m, s⁻[1:K] ≥ 0)
    v = f.values .+ s⁺ .- s⁻
    lc = _LazyConvexity(X)
    af, bf = _add_convexity!(m, lc, v)
    obj = AffExpr()
    if metric == :l2 || metric == :rms
        obj = sum(s⁺[l]^2 + s⁻[l]^2 for l ∈ 1:K)
    elseif metric == :max
        @variable(m, smax)
        obj = smax
        for l ∈ 1:K
            @constraint(m, smax ≥ s⁺[l] + s⁻[l])
        end
    else
        obj = sum(s⁺[l] + s⁻[l] for l ∈ 1:K)
    end
    @objective(m, Min, obj)

    _optimize_lazy!(m, lc, v, af, bf)

    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Unable to enforce $(curv_string) curvature: " *
            "$(termination_status(m)) ($(raw_status(m)))",
        )
    end

    return _adjusted_values(f, value.(s⁺) .- value.(s⁻), curv_string)
end

# For one-dimensional data with distinct points, the data is convex if and only
# if the slopes between consecutive points (sorted by x) are nondecreasing
function _convexify_1d(f::FunctionEvaluations, optimizer, metric, curv_string)
    K = length(f)
    x = first.(f.points)
    z = f.values
    o = sortperm(x)
    m = Model(optimizer)
    set_silent(m)
    @variable(m, δ[1:K])
    v = z .+ δ
    for j ∈ 2:K-1
        i, k, l = o[j-1], o[j], o[j+1]
        @constraint(
            m,
            (v[k] - v[i]) / (x[k] - x[i]) ≤ (v[l] - v[k]) / (x[l] - x[k])
        )
    end
    if metric == :l2 || metric == :rms
        @objective(m, Min, sum(δ .^ 2))
    elseif metric == :max
        @variable(m, t)
        @constraint(m, t .≥ δ)
        @constraint(m, t .≥ -δ)
        @objective(m, Min, t)
    else
        @variable(m, t[1:K])
        @constraint(m, t .≥ δ)
        @constraint(m, t .≥ -δ)
        @objective(m, Min, sum(t))
    end
    optimize!(m)
    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Unable to enforce $(curv_string) curvature: " *
            "$(termination_status(m)) ($(raw_status(m)))",
        )
    end
    return _adjusted_values(f, value.(δ), curv_string)
end

# Function evaluations with the values adjusted by δ. Changes at the level of
# the solver tolerance are ignored.
function _adjusted_values(f::FunctionEvaluations, δ, curv_string)
    K = length(f)
    scale = max(maximum(abs, f.values), 1.0)
    adjusted = abs.(δ) .> 1e-7 * scale
    if !any(adjusted)
        @debug "Data points are $(curv_string)"
        return f
    end
    @debug "Data are not $(curv_string), $(count(adjusted)) point(s) " *
           "adjusted, largest change = $(round(maximum(abs, δ); sigdigits = 3))"
    values = [f.values[l] + δ[l] for l ∈ 1:K]
    return FunctionEvaluations(f.points, values)
end

# Create a full order convex approximation for the given data points, i.e.
# using a separate plane for each data point. The plane for point k passes
# through the point and lies below all other points. Among such planes, the one
# with the smallest deviation (in metric) from the other points is chosen. As
# the planes are independent of each other, each plane is found separately: in
# closed form for one-dimensional data, and otherwise by a small optimization
# problem.
function _full_order_pwa(f::FunctionEvaluations{D}, optimizer, metric) where {D}
    if metric ∉ (:l1, :l2, :rms, :max, :none)
        throw(
            ArgumentError("Unsupported metric :$metric for full order fitting"),
        )
    end
    X = [pt[i] for i ∈ 1:D, pt ∈ f.points]
    z = f.values
    planes = if D == 1
        [_full_order_plane_1d(X, z, k, metric) for k ∈ eachindex(z)]
    else
        [_full_order_plane(X, z, k, metric, optimizer) for k ∈ eachindex(z)]
    end
    return PWAFunc{Convex,D}(planes)
end

function _nonconvex_error()
    return error(
        "No full order convex approximation exists, since the data points " *
        "are not convex (for example due to noise). Use " *
        "`enforce_curvature(data, Convex(), optimizer)` to adjust the " *
        "function values to be convex first.",
    )
end

# Tolerance for a plane lying below a data point, as for the feasibility
# tolerance of an optimization solver
_below_tolerance(z) = 1e-7 * max(1.0, maximum(abs, z))

# One-dimensional full order plane through point k. Its slope must lie between
# the largest slope to the points on the left and the smallest slope to the
# points on the right, and the deviation from the other points is a simple
# function of the slope within this interval.
function _full_order_plane_1d(X, z, k, metric; check = true)
    dx = X[1, :] .- X[1, k]
    dz = z .- z[k]
    lo, hi = -Inf, Inf
    for l ∈ eachindex(z)
        if dx[l] < 0
            lo = max(lo, dz[l] / dx[l])
        elseif dx[l] > 0
            hi = min(hi, dz[l] / dx[l])
        end
    end
    a = if lo > hi
        # Not convex, unless within the tolerance (checked below)
        (lo + hi) / 2
    elseif metric == :l1
        # The total deviation decreases linearly with the slope times the sum
        # of the distances to the other points
        w = sum(dx)
        w > 0 ? hi : w < 0 ? lo : isfinite(lo) ? lo : isfinite(hi) ? hi : 0.0
    elseif metric ∈ (:l2, :rms)
        clamp(dot(dx, dz) / max(dot(dx, dx), floatmin()), lo, hi)
    elseif metric == :max
        _minimax_slope(dx, dz, lo, hi)
    else
        isfinite(lo) && isfinite(hi) ? (lo + hi) / 2 :
        isfinite(lo) ? lo : isfinite(hi) ? hi : 0.0
    end
    τ = _below_tolerance(z)
    if check && any(a * dx[l] > dz[l] + τ for l ∈ eachindex(z))
        _nonconvex_error()
    end
    return Plane((a,), z[k] - a * X[1, k])
end

# The slope a in [lo, hi] minimizing the largest gap dz - a dx to the other
# points. The largest gap to the points on the left increases with a, and the
# largest gap to the points on the right decreases with a, so the optimal slope
# balances the two, and is found by bisection.
function _minimax_slope(dx, dz, lo, hi)
    isfinite(lo) || return isfinite(hi) ? hi : 0.0
    isfinite(hi) || return lo
    left(a) = maximum(dz[l] - a * dx[l] for l ∈ eachindex(dx) if dx[l] < 0)
    right(a) = maximum(dz[l] - a * dx[l] for l ∈ eachindex(dx) if dx[l] > 0)
    left(lo) ≥ right(lo) && return lo
    left(hi) ≤ right(hi) && return hi
    for _ ∈ 1:100
        a = (lo + hi) / 2
        left(a) < right(a) ? (lo = a) : (hi = a)
        hi - lo ≤ 1e-14 * max(abs(lo), abs(hi), 1.0) && break
    end
    return (lo + hi) / 2
end

# Full order plane through point k for general dimension, by an optimization
# problem with the slopes a as variables
function _full_order_plane(X, z, k, metric, optimizer; slack = 0.0)
    D, N = size(X)
    dX = X .- X[:, k]
    dz = z .- z[k]
    m = Model(optimizer)
    set_silent(m)
    @variable(m, a[1:D])
    # The gaps between the data points and the plane, which must be nonnegative
    # (up to the given slack)
    r = [dz[l] - sum(a[i] * dX[i, l] for i ∈ 1:D) for l ∈ 1:N]
    @constraint(m, r .≥ -slack)
    if metric == :l1
        @objective(m, Min, sum(r))
    elseif metric ∈ (:l2, :rms)
        @objective(m, Min, sum(rₗ^2 for rₗ ∈ r))
    elseif metric == :max
        @variable(m, t)
        @constraint(m, t .≥ r)
        @objective(m, Min, t)
    end
    optimize!(m)
    if termination_status(m) ∈ (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
        _nonconvex_error()
    elseif !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end
    â = value.(a)
    return Plane(Tuple(â), z[k] - dot(â, X[:, k]))
end

# Fit p planes, where U maps the data points to planes, minimizing the total
# gap between the points and their planes plus μ/p times the total gap between
# all points and all planes. The planes are required to lie below all data
# points, so the gap of point l to plane k is zₗ - (aₖ⋅xₗ + bₖ) ≥ 0, and the
# objective is linear in the plane coefficients and separable over the planes.
# Each plane then maximizes its value at a weighted centroid cₖ of the points,
# subject to lying below all points, i.e. it is a supporting plane of the lower
# convex envelope of the data at cₖ.
#
# With a group fitter (see _group_fitter), each plane is instead fitted to its
# points in the given metric, and the returned objective is the aggregated gap
# of the points to their planes in that metric.
function _reduced_order_pwa(
    f::FunctionEvaluations{D},
    p::Int,
    U,
    optimizer,
    μ = 1e-4;
    support = _envelope_support(f, optimizer),
    fitter = nothing,
    metric = :l1,
) where {D}
    X = [pt[i] for i ∈ 1:D, pt ∈ f.points]
    z = f.values
    N = length(z)
    if !isnothing(fitter)
        groups = [Int[] for _ ∈ 1:p]
        for (l, k) ∈ U
            push!(groups[k], l)
        end
        c̄ = vec(mean(X, dims = 2))
        planes = [isempty(g) ? support(c̄) : fitter(g) for g ∈ groups]
        gaps = [z[l] - evaluate(planes[k], X[:, l]) for (l, k) ∈ U]
        agg = metric == :max ? maximum(gaps) : sum(abs2, gaps)
        # The total gap breaks ties, so that improvements of other planes
        # than the one with the largest gap are also recognized
        return agg + 1e-6 * sum(gaps), PWAFunc{Convex,D}(planes)
    end
    # Weights and weighted centroids of the points for each plane
    n = zeros(p)
    S = zeros(D, p)
    for (l, k) ∈ U
        n[k] += 1
        S[:, k] .+= X[:, l]
    end
    w = n .+ μ * N / p
    C = (S .+ (μ / p) .* sum(X, dims = 2)) ./ w'
    planes = [support(C[:, k]) for k ∈ 1:p]
    gap(l, k) = z[l] - evaluate(planes[k], X[:, l])
    obj =
        sum(gap(l, k) for (l, k) ∈ U) +
        μ / p * sum(gap(l, k) for l ∈ 1:N, k ∈ 1:p)
    return obj, PWAFunc{Convex,D}(planes)
end

# Returns a function giving the supporting plane of the lower convex envelope of
# the data points at a point c, i.e. the plane maximizing its value at c subject
# to lying below all data points. For one-dimensional data, this is the line
# through the segment of the lower convex hull containing c.
function _envelope_support(f::FunctionEvaluations{D}, optimizer) where {D}
    X = [pt[i] for i ∈ 1:D, pt ∈ f.points]
    z = f.values
    if D == 1
        hx, hz = _lower_hull(X[1, :], z)
        return function (c)
            length(hx) == 1 && return Plane((0.0,), hz[1])
            j = clamp(searchsortedlast(hx, c[1]), 1, length(hx) - 1)
            a = (hz[j+1] - hz[j]) / (hx[j+1] - hx[j])
            return Plane((a,), hz[j] - a * hx[j])
        end
    end
    m = Model(optimizer)
    set_silent(m)
    @variable(m, a[1:D])
    @variable(m, b)
    @constraint(
        m,
        [l = 1:length(z)],
        sum(a[i] * X[i, l] for i ∈ 1:D) + b ≤ z[l],
    )
    return function (c)
        yield() # Avoid blocking of main thread causing hangs in some cases, see
        # https://github.com/sintefore/PiecewiseAffineApprox.jl/issues/11
        @objective(m, Max, sum(a[i] * c[i] for i ∈ 1:D) + b)
        optimize!(m)
        if !is_solved_and_feasible(m; allow_local = true)
            error(
                "Optimization failed: $(termination_status(m)) " *
                "($(raw_status(m)))",
            )
        end
        return Plane(Tuple(value.(a)), value(b))
    end
end

# Returns a function fitting a plane to a group of the data points (given by
# their indices), minimizing the gaps between the points of the group and the
# plane in the given metric, subject to the plane lying below all data points.
# A small multiple μ of the mean gap breaks ties. Returns nothing for :l1, where
# the planes of the weighted centroids are used, and for :l2 in two or more
# dimensions (which would require a quadratic program per group).
function _group_fitter(
    f::FunctionEvaluations{D},
    metric,
    optimizer;
    μ = 1e-4,
) where {D}
    if metric == :l1 || (metric ∈ (:l2, :rms) && D > 1)
        return nothing
    end
    X = [pt[i] for i ∈ 1:D, pt ∈ f.points]
    z = f.values
    if D == 1
        x = X[1, :]
        hx, hz = _lower_hull(x, z)
        # Slopes of the hull edges: optimal slopes lie between the extremes
        slopes = diff(hz) ./ diff(hx)
        return function (g)
            isempty(slopes) && return Plane((0.0,), hz[1])
            # The highest line with slope a below all points
            b(a) = minimum(hz[j] - a * hx[j] for j ∈ eachindex(hx))
            function obj(a)
                bₐ = b(a)
                r = (max(z[l] - a * x[l] - bₐ, 0.0) for l ∈ g)
                agg = metric == :max ? maximum(r) : sum(abs2, r)
                return agg + μ * sum(r) / length(g)
            end
            a = _golden_section(obj, minimum(slopes), maximum(slopes))
            return Plane((a,), b(a))
        end
    end
    return function (g)
        m = Model(optimizer)
        set_silent(m)
        @variable(m, a[1:D])
        @variable(m, b)
        @variable(m, t)
        plane(l) = sum(a[i] * X[i, l] for i ∈ 1:D) + b
        @constraint(m, [l ∈ eachindex(z)], plane(l) ≤ z[l])
        @constraint(m, [l ∈ g], t ≥ z[l] - plane(l))
        @objective(m, Min, t + μ * sum(z[l] - plane(l) for l ∈ g) / length(g))
        optimize!(m)
        if !is_solved_and_feasible(m; allow_local = true)
            error(
                "Optimization failed: $(termination_status(m)) " *
                "($(raw_status(m)))",
            )
        end
        return Plane(Tuple(value.(a)), value(b))
    end
end

# Minimizer of a convex function on [lo, hi] by golden section search
function _golden_section(obj, lo, hi; maxiter = 200)
    lo == hi && return lo
    φ = (sqrt(5) - 1) / 2
    c, d = hi - φ * (hi - lo), lo + φ * (hi - lo)
    fc, fd = obj(c), obj(d)
    for _ ∈ 1:maxiter
        hi - lo ≤ 1e-12 * max(abs(lo), abs(hi), 1.0) && break
        if fc ≤ fd
            hi, d, fd = d, c, fc
            c = hi - φ * (hi - lo)
            fc = obj(c)
        else
            lo, c, fc = c, d, fd
            d = lo + φ * (hi - lo)
            fd = obj(d)
        end
    end
    # Also consider the end points, as the minimum may be at the boundary
    candidates = (lo, (lo + hi) / 2, hi)
    return candidates[argmin(obj.(candidates))]
end

# Vertices of the lower convex hull of the points (x, z), sorted by x
function _lower_hull(x, z)
    hx, hz = Float64[], Float64[]
    for i ∈ sortperm(collect(zip(x, z)))
        # Only the lowest point for each x
        !isempty(hx) && x[i] == hx[end] && continue
        while length(hx) ≥ 2 &&
            (hx[end] - hx[end-1]) * (z[i] - hz[end-1]) -
            (hz[end] - hz[end-1]) * (x[i] - hx[end-1]) ≤ 0
            pop!(hx)
            pop!(hz)
        end
        push!(hx, x[i])
        push!(hz, z[i])
    end
    return hx, hz
end

# The hyperplane being active at the point x for a convex pwa function
function _active_plane(pwl::PWAFunc{Convex,D}, x) where {D}
    return argmax(collect(evaluate(p, x) for p ∈ pwl.planes))
end

# Map data points to the active plane for the piecewise approximation
function _update_allocation(f::FunctionEvaluations, pwa_red)
    U = [(l, _active_plane(pwa_red, f.points[l])) for l ∈ 1:length(f)]
    return U
end

# Do a local improvement of point allocation to planes until error
# is not further reduced
function _allocation_improvement(
    f::FunctionEvaluations,
    p,
    Uᴵ,
    optimizer;
    support = _envelope_support(f, optimizer),
    fitter = nothing,
    metric = :l1,
)
    U = Uᴵ
    prev_obj = Inf64
    fit(U, p) = _reduced_order_pwa(f, p, U, optimizer; support, fitter, metric)
    obj, pwa_red = fit(U, p)
    while obj < prev_obj
        prev_obj = obj
        U = _update_allocation(f, pwa_red)
        obj, pwa_red = fit(U, p)
    end
    # Remove planes that are not active at any data point
    active = sort(unique(last.(_update_allocation(f, pwa_red))))
    return PWAFunc(pwa_red.planes[active], Convex())
end

# Add an extra plane to the pwa: the full order plane at the point with the
# largest gap to the approximation. As this plane passes through the point, it
# is not already part of the approximation (the gap would then be zero), and
# the point is allocated to it in the updated allocation.
function _increase_order(f::FunctionEvaluations, pwa_red, pwa)
    s = [z - evaluate(pwa_red, x) for (x, z) ∈ f]
    push!(pwa_red.planes, pwa.planes[argmax(s)])
    return _update_allocation(f, pwa_red)
end

# Approximation error per data point, compared with the tolerance: the mean
# absolute deviation (:l1), the root mean square deviation (:l2, :rms) or the
# maximum absolute deviation (:max)
function _tolerance_error(f::FunctionEvaluations, pwa::PWAFunc, metric)
    r = [abs(evaluate(pwa, x) - z) for (x, z) ∈ f]
    if metric == :l1
        return mean(r)
    elseif metric == :l2 || metric == :rms
        return sqrt(mean(abs2, r))
    elseif metric == :max
        return maximum(r)
    end
    throw(
        ArgumentError(
            "Unsupported metric :$metric for Progressive, use :l1, :l2 or :max",
        ),
    )
end

# Main algorithm
function _progressive_pwa(
    f::FunctionEvaluations,
    optimizer,
    δᵗᵒˡ,
    metric,
    full_order_metric = :l1,
)
    if metric ∉ (:l1, :l2, :rms, :max)
        throw(
            ArgumentError(
                "Unsupported metric :$metric for Progressive, " *
                "use :l1, :l2 or :max",
            ),
        )
    end

    # Find the full order convex approximation, will error if the points are
    # not convex
    pwa = _full_order_pwa(f, optimizer, full_order_metric)

    # Start with an approximation with one plane and allocate all points to that segment
    p = 1
    U = [(l, 1) for l ∈ 1:length(f)]

    # Increase the number of planes until the required tolerance is met. As
    # inactive planes are removed, the number of planes does not necessarily
    # grow with each round, so the number of rounds is limited.
    # The planes are fitted to their points in the metric of the tolerance
    support = _envelope_support(f, optimizer)
    fitter = _group_fitter(f, metric, optimizer)
    improve(U, p) =
        _allocation_improvement(f, p, U, optimizer; support, fitter, metric)
    pwa_red = improve(U, p)
    rounds = 0
    while _tolerance_error(f, pwa_red, metric) > δᵗᵒˡ && rounds < 2length(f)
        rounds += 1
        U = _increase_order(f, pwa_red, pwa)
        pwa_red = improve(U, length(pwa_red.planes))
    end
    _check_tolerance(f, pwa_red, metric, δᵗᵒˡ, rounds)
    return pwa_red
end

# Warn if the fitting stopped without meeting the tolerance
function _check_tolerance(f, pwa, metric, δ, rounds)
    e = _tolerance_error(f, pwa, metric)
    p = length(pwa.planes)
    spread = maximum(f.values) - minimum(f.values)
    if e > δ + 1e-6 * (δ + spread) + 1e-12 * max(1.0, maximum(abs, f.values))
        @warn "Progressive fitting stopped after $rounds rounds without " *
              "meeting the tolerance $δ (error $e with $p planes)"
    else
        @debug "Fitting finished, error = $(round(e; sigdigits = 3)), p = $p"
    end
    return nothing
end

# Map algorithm structure to correct parameters
function _progressive_pwa(f::FunctionEvaluations, options::Progressive)
    if options.strict == :outer
        return _progressive_pwa(
            f,
            options.optimizer,
            options.tolerance,
            options.metric,
        )
    elseif options.strict ∈ (:none, :inner)
        return _progressive_band(
            f,
            options.optimizer,
            options.tolerance,
            options.strict,
        )
    end
    throw(
        ArgumentError(
            "Unsupported strict = :$(options.strict) for Progressive, " *
            "use :outer, :none or :inner",
        ),
    )
end

# Progressive fitting within a band around the data, following Kazda and Li
# (2024) specialized to convex functions. Values v within the band, such that
# the points (x, v) are convex, are found together with p planes below these
# values, minimizing the total gap between the values and the planes allocated
# to the points. When the total gap is zero, the maximum of the planes equals v
# at all data points, so the approximation lies within the band: within ±δ of
# the data (strict = :none) or between the data and δ above it (:inner).
function _progressive_band(
    f::FunctionEvaluations{D},
    optimizer,
    δ,
    strict;
    μ = 1e-4,
) where {D}
    δ ≥ 0 || throw(ArgumentError("The tolerance must be nonnegative"))
    # Work with data normalized to the order of one, so that the thresholds
    # below are independent of the scale and offset of the data
    X, μx, σx = _normalize_points(f)
    z, z̄, s = _normalize_values(f.values, δ)
    δ̂ = δ / s
    N = length(z)
    lo = strict == :inner ? copy(z) : z .- δ̂
    hi = z .+ δ̂
    # The approximation lies within the band when the largest gap between the
    # values and the planes is zero, up to the solver tolerances
    ε = max(1e-6 * δ̂, 1e-9)
    # Convexity constraints in two or more dimensions are generated lazily, and
    # kept for all the linear programs solved below
    lc = D == 1 && allunique(X[1, :]) ? nothing : _LazyConvexity(X)
    # The constraints that the planes lie below the values are also generated
    # lazily, starting with the planes of each point and its neighbours
    nbrs = _nearest_neighbours(X, 2(D + 1))
    lp(U, p) = _band_lp(X, z, lo, hi, U, p, optimizer, μ; lc, nbrs)
    U = ones(Int, N)
    v, planes, gaps = lp(U, 1)
    rounds = 0
    for _ ∈ 1:2N
        rounds += 1
        # Allocation improvement: allocate the points to the active planes
        # (removing planes that are not active at any point) and refit, as
        # long as the total gap decreases
        φ = sum(gaps)
        while maximum(gaps) > ε
            U, planes = _active_allocation(X, planes)
            v₂, planes₂, gaps₂ = lp(U, length(planes))
            sum(gaps₂) < φ * (1 - 1e-9) || break
            v, planes, gaps, φ = v₂, planes₂, gaps₂, sum(gaps₂)
        end
        maximum(gaps) ≤ ε && break
        # Order increase: add the full order planes at the points with the
        # largest gaps, without refitting, until the total gap has decreased
        φ₀ = φ
        while maximum(gaps) > ε && φ ≥ φ₀ * (1 - 1e-6)
            push!(planes, _band_plane(X, v, argmax(gaps), optimizer))
            gaps = [
                max(v[l] - maximum(evaluate(q, X[:, l]) for q ∈ planes), 0.0)
                for l ∈ 1:N
            ]
            φ = sum(gaps)
        end
        maximum(gaps) ≤ ε && break
        U, planes = _active_allocation(X, planes)
        v, planes, gaps = lp(U, length(planes))
    end
    pwa = PWAFunc{Convex,D}([_denormalize(q, μx, σx, z̄, s) for q ∈ planes])
    _check_tolerance(f, pwa, :max, δ, rounds)
    return pwa
end

# The points as columns of a matrix, shifted and scaled to [-1, 1] in each
# coordinate, together with the shift and scale
function _normalize_points(f::FunctionEvaluations{D}) where {D}
    X = [pt[i] for i ∈ 1:D, pt ∈ f.points]
    lo, hi = vec(minimum(X, dims = 2)), vec(maximum(X, dims = 2))
    μx = (lo .+ hi) ./ 2
    σx = [h > l ? (h - l) / 2 : 1.0 for (l, h) ∈ zip(lo, hi)]
    return (X .- μx) ./ σx, μx, σx
end

# The values shifted to mean zero and scaled to a largest absolute value of one
# (or by the tolerance δ for constant values), together with the shift and scale
function _normalize_values(z, δ)
    z̄ = mean(z)
    s = maximum(abs, z .- z̄)
    s = s > 0 ? s : δ > 0 ? δ : 1.0
    return (z .- z̄) ./ s, z̄, s
end

# The plane ẑ = â⋅x̂ + b̂ for normalized points x̂ = (x - μx) ./ σx and values
# ẑ = (z - z̄) / s, in the original coordinates
function _denormalize(q::Plane, μx, σx, z̄, s)
    â = collect(q.α)
    a = s .* â ./ σx
    b = s * (q.β - sum(â .* μx ./ σx)) + z̄
    return Plane(Tuple(a), b)
end

# Allocate each point to the plane that is active at the point, and remove the
# planes that are not active at any point
function _active_allocation(X, planes)
    U = [argmax([evaluate(q, X[:, l]) for q ∈ planes]) for l ∈ axes(X, 2)]
    active = sort(unique(U))
    index = Dict(k => i for (i, k) ∈ enumerate(active))
    return [index[k] for k ∈ U], planes[active]
end

# A plane through the point (xᵢ, vᵢ) and below all other points (xₗ, vₗ), where
# the values v are convex up to the solver tolerance
function _band_plane(X, v, i, optimizer)
    if size(X, 1) == 1
        return _full_order_plane_1d(X, v, i, :l1; check = false)
    end
    τ = _below_tolerance(v)
    return _full_order_plane(X, v, i, :l1, optimizer; slack = τ)
end

# The linear program for given allocation U of the data points to p planes.
# Returns the values v, the planes and the gaps between the values and the
# allocated planes. Among solutions with the same total gap, values close to
# the data are preferred (with weight μ). Without lazy convexity constraints
# lc, the data must be one-dimensional with distinct points. With neighbours
# nbrs of the points, the constraints that the planes lie below the values are
# generated lazily: initially for the planes allocated to each point and its
# neighbours, and then for violated pairs until none remain.
function _band_lp(
    X,
    z,
    lo,
    hi,
    U,
    p,
    optimizer,
    μ;
    lc = nothing,
    nbrs = nothing,
)
    D, N = size(X)
    m = Model(optimizer)
    set_silent(m)
    @variable(m, lo[l] ≤ v[l = 1:N] ≤ hi[l])
    if isnothing(lc)
        _convexity_constraints_1d!(m, X, v)
    else
        af, bf = _add_convexity!(m, lc, v)
    end
    @variable(m, a[1:p, 1:D])
    @variable(m, b[1:p])
    plane(k, l) = sum(a[k, i] * X[i, l] for i ∈ 1:D) + b[k]
    below = if isnothing(nbrs)
        [Set(1:p) for _ ∈ 1:N]
    else
        [Set(U[j] for j ∈ (l, nbrs[l]...)) for l ∈ 1:N]
    end
    for l ∈ 1:N, k ∈ below[l]
        @constraint(m, plane(k, l) ≤ v[l])
    end
    @variable(m, t[1:N])
    @constraint(m, t .≥ v .- z)
    @constraint(m, t .≥ z .- v)
    @objective(m, Min, sum(v[l] - plane(U[l], l) for l ∈ 1:N) + μ * sum(t))
    while true
        optimize!(m)
        is_solved_and_feasible(m; allow_local = true) || break
        # Planes above the values of points without a constraint for them (all
        # solution values are queried before the model is modified)
        v̂ = value.(v)
        P = value.(a) * X .+ value.(b)
        τ = 1e-7 * max(1.0, maximum(abs, v̂))
        new = Tuple{Int,Int}[]
        for l ∈ 1:N
            viol = [
                (P[k, l] - v̂[l], k) for
                k ∈ 1:p if P[k, l] - v̂[l] > τ && k ∉ below[l]
            ]
            append!(new, ((l, k) for (_, k) ∈ first(sort!(viol; rev = true), 3)))
        end
        added = isnothing(lc) ? 0 : _add_violated_convexity!(m, lc, v, af, bf)
        for (l, k) ∈ new
            push!(below[l], k)
            @constraint(m, plane(k, l) ≤ v[l])
        end
        added + length(new) == 0 && break
    end
    if termination_status(m) ∈ (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
        error(
            "No convex function lies within the tolerance of the data. " *
            "Increase the tolerance, or use `enforce_curvature` first.",
        )
    elseif !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end
    v̂ = value.(v)
    planes = [Plane(Tuple(value.(a[k, :])), value(b[k])) for k ∈ 1:p]
    gaps = [max(v̂[l] - evaluate(planes[U[l]], X[:, l]), 0.0) for l ∈ 1:N]
    return v̂, planes, gaps
end

# Constraints ensuring that the one-dimensional points (x, v), with distinct x,
# are convex: nondecreasing slopes between consecutive points
function _convexity_constraints_1d!(m, X, v)
    N = size(X, 2)
    x = X[1, :]
    o = sortperm(x)
    for j ∈ 2:N-1
        i, k, l = o[j-1], o[j], o[j+1]
        @constraint(
            m,
            (v[k] - v[i]) / (x[k] - x[i]) ≤ (v[l] - v[k]) / (x[l] - x[k])
        )
    end
    return
end

# Convexity of the points (x, v) in general: a plane through each point
# (xₖ, vₖ) that lies below all other points. As most of these N² constraints are
# not binding, they are generated lazily, as in the relaxed linear program of
# Kazda and Li (2024): initially for the nearest neighbours of each point, and
# then for violated pairs until none remain. The generated pairs (k, l) are
# kept, and reused for further problems with the same points.
struct _LazyConvexity
    X::Matrix{Float64}
    pairs::Set{Tuple{Int,Int}}
end

function _LazyConvexity(X::AbstractMatrix; neighbours = 2(size(X, 1) + 1))
    nbrs = _nearest_neighbours(X, neighbours)
    pairs = Set((k, l) for k ∈ axes(X, 2) for l ∈ nbrs[k])
    return _LazyConvexity(Matrix{Float64}(X), pairs)
end

# Adds the planes through the points, and the generated constraints that the
# planes lie below the other points, to the model. Returns the slopes and
# intercepts of the planes.
function _add_convexity!(m, lc::_LazyConvexity, v)
    D, N = size(lc.X)
    af = @variable(m, [1:N, 1:D])
    bf = @variable(m, [1:N])
    for k ∈ 1:N
        @constraint(m, sum(af[k, i] * lc.X[i, k] for i ∈ 1:D) + bf[k] == v[k])
    end
    for (k, l) ∈ lc.pairs
        @constraint(m, sum(af[k, i] * lc.X[i, l] for i ∈ 1:D) + bf[k] ≤ v[l])
    end
    return af, bf
end

# Optimizes the model, adding the violated convexity constraints and
# reoptimizing until none are violated
function _optimize_lazy!(m, lc::_LazyConvexity, v, af, bf)
    while true
        optimize!(m)
        is_solved_and_feasible(m; allow_local = true) || return
        _add_violated_convexity!(m, lc, v, af, bf) == 0 && return
    end
end

# Adds the violated convexity constraints for the current solution, at most
# maxadd per plane, and returns the number of constraints added
function _add_violated_convexity!(
    m,
    lc::_LazyConvexity,
    v,
    af,
    bf;
    maxadd = 10,
)
    D, N = size(lc.X)
    v̂ = value.(v)
    # P[k, l]: value of plane k at point l, compared with the value of l
    P = value.(af) * lc.X .+ value.(bf)
    τ = 1e-7 * max(1.0, maximum(abs, v̂))
    added = 0
    for k ∈ 1:N
        viol = [
            (P[k, l] - v̂[l], l) for
            l ∈ 1:N if l != k && P[k, l] - v̂[l] > τ && (k, l) ∉ lc.pairs
        ]
        for (_, l) ∈ first(sort!(viol; rev = true), maxadd)
            push!(lc.pairs, (k, l))
            @constraint(m, sum(af[k, i] * lc.X[i, l] for i ∈ 1:D) + bf[k] ≤ v[l])
            added += 1
        end
    end
    return added
end

# The indices of the n nearest neighbours of each point (columns of X), with
# the coordinates scaled to the same range
function _nearest_neighbours(X, n)
    D, N = size(X)
    r = vec(maximum(X, dims = 2) .- minimum(X, dims = 2))
    Y = X .* [rᵢ > 0 ? 1 / rᵢ : 1.0 for rᵢ ∈ r]
    return map(1:N) do k
        d = [sum(abs2, Y[:, l] .- Y[:, k]) for l ∈ 1:N]
        return [l for l ∈ partialsortperm(d, 1:min(N, n + 1)) if l != k]
    end
end
