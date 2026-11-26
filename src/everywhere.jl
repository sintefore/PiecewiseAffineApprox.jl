# Approximations of a known function that lie on one side of it over the whole
# domain, not only at the sampled points, and the largest deviations between
# an approximation and the function over the domain.
#
# Both build on the curvature of the function. For a convex function f and a
# convex piecewise affine p = maxₖ pₖ over a box:
#
# - f - p is convex on each piece of p, so its maximum over the box is at a
#   vertex of p (the vertices of the pieces, including where they meet the
#   boundary of the box). Evaluating f at these finitely many points gives the
#   maximum exactly.
# - p - f = maxₖ (pₖ - f), and each pₖ - f is concave on the box, so its
#   maximum is found by a local search.
#
# Concave functions are handled by negation.

# The value of f at a point given as a tuple, see _evaluate_points
_call(f, x::Tuple) = applicable(f, x) ? f(x) : f(x...)

"""
    error_bounds(pwa::PWAFunc, f::Function, bounds = domain(pwa))

The largest deviations between the approximation `pwa` and the function `f`
over the box `bounds = [(lo₁, hi₁), …]`, not only at sampled points. Returns a
named tuple with

- `below`: the maximum of `f - pwa`, the most the approximation lies below `f`,
- `above`: the maximum of `pwa - f`, the most it lies above `f`,
- `x_below`, `x_above`: points where these are attained.

A negative value means that the approximation is on the other side of `f`
everywhere, by at least that much. The largest absolute error over the box is
`max(below, above)`.

The function `f` must have the curvature of `pwa` on the box (convex for a
convex approximation); this is not checked. For a convex approximation, `below`
is then exact up to rounding: `f` is evaluated at the vertices of `pwa`.
`above` is found by a search for the maximum of
each plane minus `f`, which is concave: a golden section search in one
variable, and searches along the coordinates in several, which requires `f` to
be differentiable there. For a concave approximation the roles are
interchanged.
"""
function error_bounds(
    pwa::PWAFunc{Convex,D},
    f::Function,
    bounds = domain(pwa),
) where {D}
    isnothing(bounds) && throw(
        ArgumentError(
            "The approximation has no domain, give the bounds as " *
            "[(lo, hi), ...] with one tuple for each variable",
        ),
    )
    lo, hi = _check_bounds(bounds, D)
    x_below, below = _max_below(pwa, f, lo, hi)
    above, x_above = -Inf, Tuple(lo)
    for p ∈ pwa.planes
        x, v = _max_plane_above(p, f, lo, hi)
        if v > above
            above, x_above = v, x
        end
    end
    return (; below, above, x_below, x_above)
end

function error_bounds(
    pwa::PWAFunc{Concave,D},
    f::Function,
    bounds = domain(pwa),
) where {D}
    # In terms of the convex function -pwa, see `evaluate`
    r = error_bounds(_as_convex(pwa), x -> -_call(f, x), bounds)
    return (
        below = r.above,
        above = r.below,
        x_below = r.x_above,
        x_above = r.x_below,
    )
end

# The convex function max of the planes of a concave function, which is its
# negative
_as_convex(pwa::PWAFunc{Concave,D,T}) where {D,T} =
    PWAFunc{Convex,D,T}(pwa.planes, pwa.domain)
_as_concave(pwa::PWAFunc{Convex,D,T}) where {D,T} =
    PWAFunc{Concave,D,T}(pwa.planes, pwa.domain)

# The maximum of f - pwa over the box for convex f, attained at a vertex of
# pwa, and all vertices with their values of f - pwa
function _vertex_gaps(pwa::PWAFunc{Convex,D}, f, lo, hi) where {D}
    verts = _enumerate_vertices_box(pwa, lo, hi, 1e-9)
    return [(Tuple(x), _call(f, Tuple(x)) - z) for (x, z) ∈ verts]
end
function _max_below(pwa::PWAFunc{Convex,D}, f, lo, hi) where {D}
    gaps = _vertex_gaps(pwa, f, lo, hi)
    x, v = gaps[argmax(last.(gaps))]
    return x, v
end

# The maximum over the box of the plane minus the convex function f
function _max_plane_above(p::Plane{D}, f, lo, hi) where {D}
    φ(x) = evaluate(p, x) - _call(f, x)
    return _maximize_concave(φ, lo, hi)
end

# Maximize a concave function over a box by golden section searches along the
# coordinates, followed after each sweep by a search along the direction of
# the sweep. Exact in one variable. In several, the search can only stop short
# of the maximum where the function is not differentiable.
function _maximize_concave(φ, lo::Vector{Float64}, hi::Vector{Float64})
    D = length(lo)
    x = (lo .+ hi) ./ 2
    at(i, t) = Tuple(j == i ? t : x[j] for j ∈ 1:D)
    for i ∈ 1:D
        x[i] = _golden_section(t -> -φ(at(i, t)), lo[i], hi[i])
    end
    v = φ(Tuple(x))
    D == 1 && return Tuple(x), v
    scale = max(abs(v), 1.0)
    for _ ∈ 1:200
        x₀, v₀ = copy(x), v
        for i ∈ 1:D
            x[i] = _golden_section(t -> -φ(at(i, t)), lo[i], hi[i])
        end
        # Along the direction of the sweep, within the box
        d = x .- x₀
        if any(!iszero, d)
            t̲, t̄ = -Inf, Inf
            for i ∈ 1:D
                d[i] == 0 && continue
                a, b = (lo[i] - x[i]) / d[i], (hi[i] - x[i]) / d[i]
                t̲, t̄ = max(t̲, min(a, b)), min(t̄, max(a, b))
            end
            line(t) = Tuple(clamp.(x .+ t .* d, lo, hi))
            t = _golden_section(t -> -φ(line(t)), t̲, t̄)
            φ(line(t)) > φ(Tuple(x)) && (x = collect(line(t)))
        end
        v = φ(Tuple(x))
        v - v₀ ≤ 1e-12 * scale && break
    end
    return Tuple(x), v
end

# Fit the function with the algorithm such that the side given by its `strict`
# option holds over the whole box. The fit is made in terms of a convex
# function (concave functions are negated), for which :inner is above the
# function and :outer below it.
function _fit_everywhere(
    f::Function,
    bbox,
    data::FunctionEvaluations{D},
    c::Curvature,
    a::ApproximationAlgorithm,
    refinements::Int,
) where {D}
    c isa DC && throw(
        ArgumentError(
            "everywhere = true requires a convex or concave function, as " *
            "the side of a difference of convex approximation between the " *
            "data points cannot be determined from its vertices",
        ),
    )
    refinements ≥ 0 || throw(ArgumentError("refinements must be nonnegative"))
    side = strict(a)
    lo, hi = _check_bounds(bbox, D)
    σ = c isa Concave ? -1 : 1
    fc(x) = σ * _call(f, x)
    data = FunctionEvaluations(data.points, σ .* data.values)
    scale = max(maximum(data.values) - minimum(data.values), eps())
    tol = 1e-9 * scale

    local r, viol
    rounds = 0
    t = @elapsed begin
        r = fit(data, Convex(), a)
        viol = _side_violations(r.func, fc, lo, hi, side)
        while rounds < refinements && !isempty(viol.points) && viol.largest > tol
            # Refit with the points of the violations added to the data
            new = setdiff(viol.points, data.points)
            isempty(new) && break
            refined = FunctionEvaluations(
                [data.points; new],
                [data.values; [fc(x) for x ∈ new]],
            )
            @debug "Refinement $(rounds + 1): $(length(new)) points added, " *
                   "largest violation $(viol.largest)"
            r′ = fit(refined, Convex(), a)
            viol′ = _side_violations(r′.func, fc, lo, hi, side)
            rounds += 1
            # A refit that is not closer to the side is not used
            viol′.largest < viol.largest || break
            data, r, viol = refined, r′, viol′
        end
        func = _shift_to_side(r.func, viol, side, tol)
        bounds = error_bounds(func, fc, bbox)
    end
    func = _with_domain(func, _data_domain(data))
    data = FunctionEvaluations(data.points, σ .* data.values)
    if c isa Concave
        func = _as_concave(func)
        bounds = (
            below = bounds.above,
            above = bounds.below,
            x_below = bounds.x_above,
            x_above = bounds.x_below,
        )
    end
    details = (;
        r.details...,
        refinements = rounds,
        shift = maximum(viol.shifts; init = 0.0),
        points = length(data.points),
        bounds,
    )
    return _result(data, a, func, details, t)
end

# The violations of the side by a convex approximation of the convex function
# f over the box: the points to add to the data, the largest violation, and the
# shift of each plane that removes them
function _side_violations(pwa::PWAFunc{Convex,D}, f, lo, hi, side) where {D}
    K = length(pwa.planes)
    if side == :inner
        # Above f: the largest violation is at a vertex, and all planes are
        # shifted by it
        gaps = _vertex_gaps(pwa, f, lo, hi)
        largest = max(maximum(last, gaps), 0.0)
        points = [x for (x, v) ∈ gaps if v > 0]
        return (; points, largest, shifts = fill(largest, K))
    elseif side == :outer
        # Below f: each plane is shifted by its own largest violation
        found = [_max_plane_above(p, f, lo, hi) for p ∈ pwa.planes]
        shifts = [max(v, 0.0) for (_, v) ∈ found]
        points = unique([x for (x, v) ∈ found if v > 0])
        return (; points, largest = maximum(shifts), shifts)
    end
    return (; points = NTuple{D,Float64}[], largest = 0.0, shifts = zeros(K))
end

# The planes shifted to the side, with a margin for the rounding in the search
function _shift_to_side(pwa::PWAFunc{Convex,D,T}, viol, side, margin) where {D,T}
    (side == :none || viol.largest ≤ 0) && return pwa
    sgn = side == :inner ? 1 : -1
    planes = map(pwa.planes, viol.shifts) do p, s
        return s > 0 ? Plane{D,T}(p.α, p.β + sgn * (s + margin)) : p
    end
    return PWAFunc{Convex,D,T}(planes, pwa.domain)
end
