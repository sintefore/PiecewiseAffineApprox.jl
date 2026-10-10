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

    m = Model(optimizer)
    set_silent(m)
    @variable(m, a[1:K, 1:D])
    @variable(m, b[1:K])
    @variable(m, s⁺[1:K] ≥ 0)
    @variable(m, s⁻[1:K] ≥ 0)

    for (l, (x, z)) ∈ enumerate(f)
        for k ∈ 1:K
            @constraint(m, z + s⁺[l] - s⁻[l] ≥ dot(a[k, :], x) + b[k])
        end
        @constraint(m, z + s⁺[l] - s⁻[l] ≤ dot(a[l, :], x) + b[l])
    end
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

    optimize!(m)

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
function _full_order_plane_1d(X, z, k, metric)
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
        w > 0 ? hi : w < 0 ? lo : (isfinite(lo) ? lo : hi)
    elseif metric ∈ (:l2, :rms)
        clamp(dot(dx, dz) / max(dot(dx, dx), floatmin()), lo, hi)
    elseif metric == :max
        _minimax_slope(dx, dz, lo, hi)
    else
        isfinite(lo) && isfinite(hi) ? (lo + hi) / 2 :
        isfinite(lo) ? lo : isfinite(hi) ? hi : 0.0
    end
    τ = _below_tolerance(z)
    any(a * dx[l] > dz[l] + τ for l ∈ eachindex(z)) && _nonconvex_error()
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
function _full_order_plane(X, z, k, metric, optimizer)
    D, N = size(X)
    dX = X .- X[:, k]
    dz = z .- z[k]
    m = Model(optimizer)
    set_silent(m)
    @variable(m, a[1:D])
    # The gaps between the data points and the plane, which must be nonnegative
    r = [dz[l] - sum(a[i] * dX[i, l] for i ∈ 1:D) for l ∈ 1:N]
    @constraint(m, r .≥ 0)
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
function _reduced_order_pwa(
    f::FunctionEvaluations{D},
    p::Int,
    U,
    optimizer,
    μ = 1e-4;
    support = _envelope_support(f, optimizer),
) where {D}
    X = [pt[i] for i ∈ 1:D, pt ∈ f.points]
    z = f.values
    N = length(z)
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
)
    U = Uᴵ
    prev_obj = Inf64
    obj, pwa_red = _reduced_order_pwa(f, p, U, optimizer; support)
    while obj < prev_obj
        prev_obj = obj
        U = _update_allocation(f, pwa_red)
        obj, pwa_red = _reduced_order_pwa(f, p, U, optimizer; support)
    end
    return pwa_red
end

# Add an extra plane to the pwa by identifying the point lying
# further away and use the corresponding plane in the full order pwa
# to update the allocation of points to planes.
function _increase_order(f::FunctionEvaluations, pwa_red, pwa, used)
    s = [z - evaluate(pwa_red, x) for (x, z) ∈ f]
    imax = 0
    smax = -1
    for i ∈ 1:length(f)
        if s[i] > smax && !(i in used)
            imax = i
            smax = s[i]
        end
    end
    push!(pwa_red.planes, pwa.planes[imax])
    U = _update_allocation(f, pwa_red)
    return U, vcat(used, imax)
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

    # Increase the number of planes until the required tolerance is met
    support = _envelope_support(f, optimizer)
    pwa_red = _allocation_improvement(f, p, U, optimizer; support)
    used = Int[]
    while _tolerance_error(f, pwa_red, metric) > δᵗᵒˡ && p < length(f.points)
        U, used = _increase_order(f, pwa_red, pwa, used)
        p += 1
        pwa_red = _allocation_improvement(f, p, U, optimizer; support)
    end
    e = _tolerance_error(f, pwa_red, metric)
    @debug "Fitting finished, error = $(round(e; sigdigits = 3)), p = $p"
    return pwa_red
end

# Map algorithm structure to correct parameters
function _progressive_pwa(f::FunctionEvaluations, options::Progressive)
    return _progressive_pwa(
        f,
        options.optimizer,
        options.tolerance,
        options.metric,
    )
end
