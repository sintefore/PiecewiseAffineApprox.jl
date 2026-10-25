# Approximation of general (non-convex) functions as a difference of convex
# piecewise linear functions, f ≈ g - h, represented as the sum of a convex and
# a concave part, see DCFunc.
#
# The generic approach works with any algorithm for convex approximation: the
# data is split into a convex and a concave part by a linear program, each part
# is fitted with the algorithm, and the two parts are refitted alternately.

# The default number of alternating refits for an algorithm: none for the
# algorithms that require convex data
_dc_iterations(a::ApproximationAlgorithm) = requires_convex_data(a) ? 0 : 5

# The metric used for comparing the alternating refits
function _dc_metric(a::ApproximationAlgorithm)
    m = metric(a)
    return isnothing(m) || m ∉ (:l1, :l2, :rms, :max) ? :l1 : m
end

"""
    approx(input::FunctionEvaluations{D}, c::DC, a::ApproximationAlgorithm)

Return a [`DCFunc`](@ref), approximating the `input` points in `D` dimensions
as the sum of a convex and a concave piecewise linear function. Any function
can be approximated in this way, so the points do not need to be convex or
concave.

The algorithm `a` is used for both parts, as described for
[`Alternating`](@ref), with `a.planes` planes in the convex part and
`a.concave_planes` in the concave part for the algorithms that have these
options. Use `Alternating` to fit the parts with different algorithms. The
`strict` option of the algorithm applies to the approximation as a whole: with
`:outer` it lies below the data and with `:inner` above it.
"""
function approx(
    input::FunctionEvaluations{D},
    c::DC,
    a::ApproximationAlgorithm,
) where {D}
    return _dc_approx(input, a, a, optimizer(a), _dc_iterations(a))
end

function approx(input::FunctionEvaluations{D}, c::DC, a::Alternating) where {D}
    opt = a.optimizer
    if isnothing(opt)
        opt = optimizer(a.convex)
    end
    if isnothing(opt)
        opt = optimizer(a.concave)
    end
    return _dc_approx(input, a.convex, a.concave, opt, a.iterations)
end

function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Alternating,
) where {D}
    return _alternating_error(c)
end
approx(::FunctionEvaluations, c::Concave, ::Alternating) = _alternating_error(c)
function _alternating_error(c)
    return throw(
        ArgumentError(
            "Alternating is only available for DC(), use its convex or " *
            "concave algorithm directly for $(typeof(c))()",
        ),
    )
end

function _default_samples(a::Alternating, D)
    return max(_default_samples(a.convex, D), _default_samples(a.concave, D))
end

# The algorithm for fitting the concave part: strict = :outer and :inner are
# interchanged, and the number of planes is concave_planes
function _dc_concave_algorithm(a::ApproximationAlgorithm)
    a = _dc_swap_strict(a)
    p = planes(a)
    (isnothing(p) || concave_planes(a) == p) && return a
    return with(a; planes = concave_planes(a))
end

# The algorithm with strict = :outer and :inner interchanged, for fitting the
# concave part. Algorithms without strictness are returned as they are.
function _dc_swap_strict(a::ApproximationAlgorithm)
    s = strict(a)
    swapped = s == :outer ? :inner : s == :inner ? :outer : s
    swapped == s && return a
    return with(a; strict = swapped)
end

# Fit the convex part with ag and the concave part with ah, alternately. The
# strictness of the algorithms refers to the sum of the parts: with :outer,
# the sum lies below the data. Each part is fitted to the data minus the other
# part, so the convex part must then lie below these values, and the concave
# part as well, which for a concave approximation is strict = :inner.
function _dc_approx(
    f::FunctionEvaluations{D},
    ag::ApproximationAlgorithm,
    ah::ApproximationAlgorithm,
    optimizer,
    iterations,
) where {D}
    ah = _dc_concave_algorithm(ah)
    metric = _dc_metric(ag)
    x = f.points
    z = f.values
    function residual(pwa)
        return FunctionEvaluations(
            x,
            [z[l] - evaluate(pwa, x[l]) for l ∈ eachindex(x)],
        )
    end
    if isnothing(optimizer) && !_dc_closed_form(f)
        # Without a decomposition, start from a convex fit to the data
        g = approx(f, Convex(), ag)
        h = approx(residual(g), Concave(), ah)
    else
        # The split uses the default metric (a linear program), independent
        # of the metric of the fits
        parts = dc_decomposition(f, optimizer)
        g = approx(parts.convex, Convex(), ag)
        h = approx(parts.concave, Concave(), ah)
    end
    best = DCFunc(g, h)
    eᵇᵉˢᵗ = _tolerance_error(f, best, metric)
    eᵖʳᵉᵛ = eᵇᵉˢᵗ
    for _ ∈ 1:iterations
        g = approx(residual(h), Convex(), ag)
        h = approx(residual(g), Concave(), ah)
        dc = DCFunc(g, h)
        e = _tolerance_error(f, dc, metric)
        if e < eᵇᵉˢᵗ
            best, eᵇᵉˢᵗ = dc, e
        end
        # Stop when the refits no longer change the error (a fixed point)
        abs(e - eᵖʳᵉᵛ) ≤ 1e-9 * max(e, eᵖʳᵉᵛ) && break
        eᵖʳᵉᵛ = e
    end
    return best
end

# Whether the decomposition can be computed without an optimizer: for
# one-dimensional data with distinct points
function _dc_closed_form(f::FunctionEvaluations{D}) where {D}
    return D == 1 && allunique(first.(f.points))
end

"""
    dc_decomposition(f::FunctionEvaluations, optimizer = nothing; metric = :l1)

Split the function values into a convex and a concave part, returning a named
tuple `(convex, concave)` of `FunctionEvaluations` at the same points with
`convex.values + concave.values == f.values`, such that the convex part can be
interpolated by a convex function and the concave part by a concave function.

Such a decomposition always exists for distinct points, but is not unique. The
concave part is chosen with as little curvature as possible, since it is the
expensive part to model, and its affine component is moved to the convex part,
so for convex data the concave part is zero.

In one dimension with distinct points, the concave part collects the decreases
of the slope between consecutive points, and is computed directly without an
`optimizer`. Otherwise, the concave part is chosen as close to affine as
possible, measured by `metric`: the sum of absolute deviations (`:l1`), the sum
of squared deviations (`:l2`, a quadratic program) or the maximum absolute
deviation (`:max`), by an optimization problem with the convexity constraints
generated lazily, as in [`enforce_curvature`](@ref).
"""
function dc_decomposition(
    f::FunctionEvaluations{D},
    optimizer = nothing;
    metric = :l1,
) where {D}
    if metric ∉ (:l1, :l2, :rms, :max)
        throw(
            ArgumentError(
                "Unsupported metric :$metric for dc_decomposition, " *
                "use :l1, :l2 or :max",
            ),
        )
    end
    _dc_closed_form(f) && return _dc_decomposition_1d(f)
    if isnothing(optimizer)
        throw(
            ArgumentError(
                "dc_decomposition requires an optimizer, except for " *
                "one-dimensional data with distinct points",
            ),
        )
    end
    # Work with normalized data for the solver tolerances to be independent of
    # the scale of the data
    X, _, _ = _normalize_points(f)
    z, z̄, s = _normalize_values(f.values, 0.0)
    N = length(z)

    m = Model(optimizer)
    set_silent(m)
    # The convex part z + vₕ and the concave part -vₕ, with vₕ convex
    @variable(m, vₕ[1:N])
    v₉ = z .+ vₕ
    lazy = !(D == 1 && allunique(X[1, :]))
    if lazy
        lc₉, lcₕ = _LazyConvexity(X), _LazyConvexity(X)
        a₉, b₉ = _add_convexity!(m, lc₉, v₉)
        aₕ, bₕ = _add_convexity!(m, lcₕ, vₕ)
    else
        _convexity_constraints_1d!(m, X, v₉)
        _convexity_constraints_1d!(m, X, vₕ)
    end
    # The deviation of vₕ from an affine function c⋅x + d
    @variable(m, c[1:D])
    @variable(m, d)
    r = [vₕ[l] - sum(c[i] * X[i, l] for i ∈ 1:D) - d for l ∈ 1:N]
    if metric == :l2 || metric == :rms
        @objective(m, Min, sum(r .^ 2))
    elseif metric == :max
        @variable(m, t)
        @constraint(m, t .≥ r)
        @constraint(m, t .≥ -r)
        @objective(m, Min, t)
    else
        @variable(m, t[1:N])
        @constraint(m, t .≥ r)
        @constraint(m, t .≥ -r)
        @objective(m, Min, sum(t))
    end

    while true
        optimize!(m)
        (lazy && is_solved_and_feasible(m; allow_local = true)) || break
        # All solution values are queried before the model is modified
        new₉ = _violated_convexity(lc₉, value.(v₉), value.(a₉), value.(b₉))
        newₕ = _violated_convexity(lcₕ, value.(vₕ), value.(aₕ), value.(bₕ))
        _add_convexity_pairs!(m, lc₉, v₉, a₉, b₉, new₉)
        _add_convexity_pairs!(m, lcₕ, vₕ, aₕ, bₕ, newₕ)
        length(new₉) + length(newₕ) == 0 && break
    end
    if termination_status(m) ∈ (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
        error(
            "No difference of convex decomposition exists for the data, " *
            "possibly due to repeated points with different values",
        )
    elseif !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end

    # Move the affine component of the concave part to the convex part, and
    # ignore deviations at the level of the solver tolerance
    ĉ, d̂ = value.(c), value(d)
    v̂ₕ = [value(vₕ[l]) - sum(ĉ[i] * X[i, l] for i ∈ 1:D) - d̂ for l ∈ 1:N]
    if maximum(abs, v̂ₕ) ≤ 1e-7
        v̂ₕ .= 0.0
    end
    convex = FunctionEvaluations(f.points, z̄ .+ s .* (z .+ v̂ₕ))
    concave = FunctionEvaluations(f.points, -s .* v̂ₕ)
    return (; convex, concave)
end

# For sorted distinct points, the data is convex if and only if the slopes
# between consecutive points are nondecreasing. The concave part is the convex
# function (negated) whose slope increases are the decreases of the slopes of
# the data, so that the convex part has slope changes max(Δ, 0). Its affine
# component (the least squares fit) is moved to the convex part.
function _dc_decomposition_1d(f::FunctionEvaluations)
    x = first.(f.points)
    z = f.values
    N = length(z)
    o = sortperm(x)
    vₕ = zeros(N)
    if N > 2
        xₒ, zₒ = x[o], z[o]
        slopes = diff(zₒ) ./ diff(xₒ)
        # Slope of vₕ on each segment, starting from zero
        t = cumsum([0.0; max.(.-diff(slopes), 0.0)])
        vₒ = cumsum([0.0; t .* diff(xₒ)])
        vₕ[o] = vₒ
        A = [x ones(N)]
        vₕ .-= A * (A \ vₕ)
    end
    convex = FunctionEvaluations(f.points, z .+ vₕ)
    concave = FunctionEvaluations(f.points, -vₕ)
    return (; convex, concave)
end
