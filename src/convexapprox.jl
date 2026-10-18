defaultmetric() = :l1
defaultplanes() = 4
defaulttimelimit() = 60

"""
    approx(input::FunctionEvaluations{D}, c::Curvature, a::Algorithm)

Return PWAFunc{Convex,D} or PWAFunc{Concave,D} depending on `c`, approximating the `input` points in `D` dimensions
"""
function approx(input, c::Concave, a::Algorithm)
    cv = approx(FunctionEvaluations(input.points, -input.values), Convex(), a;)
    return PWAFunc{Concave,dims(cv)}(cv.planes)
end

dims(pwa::PWAFunc{C,D}) where {C,D} = D

# Using dispatch for specializing on dimensions. If performance were a concern,
# maybe just do branching and call specialized function directly
function approx(input::FunctionEvaluations{D}, c::Convex, a;) where {D}
    return approx(input, c, a, Val(D))
end

# Specialized for 1D and interpolation
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Interpol,
    ::Val{1};
) where {D}
    return _convex_linearization_ipol(
        [i[1] for i ∈ input.points],
        input.values,
        a;
    )
end

function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Interpol,
    ::Val,
) where {D}
    throw(
        ArgumentError(
            "Interpol only supports 1D data, got $D dimensions. " *
            "Use another algorithm, e.g. Cluster.",
        ),
    )
end

"""
    approx(input::FunctionEvaluations, c::Convex, a::Cluster)

Approximate using a heuristic that works for general dimensions.

Additional keyword arguments:
- `trials`: number of restarts (default = 20)
- `itlim`: max refining iterations in each trial (default = 50),
"""
function approx(input::FunctionEvaluations{D}, c::Convex, a::Cluster;) where {D}
    x = [p[i] for i ∈ 1:D, p ∈ input.points]
    z = input.values
    return _convex_linearization_mb(x, z, a)
end

"""
    approx(input::FunctionEvaluations, c::Convex, a::Progressive)

Approximate using a progressive fitting heuristic that adds planes until
a specified error tolerance is met.

This algorithm requires that the data points provided are samples from
a convex function.
"""
function approx(input::FunctionEvaluations, c::Convex, a::Progressive;)
    return _progressive_pwa(input, a)
end

"""
    approx(input::FunctionEvaluations, c::Convex, a::FullOrder)

Approximate using full order fitting that adds planes for all sample points.

This algorithm requires that the data points provided are samples from
a convex function.
"""
function approx(input::FunctionEvaluations, c::Convex, a::FullOrder;)
    return _full_order_pwa(input, a.optimizer, a.metric)
end

# Optimal convex approximation using mixed integer optimization
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    options::MILP,
) where {D}
    _validate(options)
    x = input.points
    z = input.values
    𝒫 = eachindex(x)
    𝒦 = 1:options.planes
    ℐₚ = 1:D

    # A heuristic solution is used as a start for the solver, to bound the
    # deviation of the approximation from the data, and to make sure that the
    # big-M values do not cut it off
    start = options.warmstart ? _milp_start(input, options) : nothing

    # Big-M for each point, bounding the distance from the approximation down to
    # the planes that are not active in the point
    Mᵖ = if D ≤ 2
        2 .* _linear_big_M_points(input)
    elseif !isnothing(start)
        fill(2 * maximum(start.gap), length(x))
    else
        throw(
            ArgumentError(
                "MILP without warmstart only supports 1D and 2D, got $D dimensions",
            ),
        )
    end
    if !isnothing(start)
        Mᵖ = max.(Mᵖ, start.gap)
    end

    m = Model()

    set_time_limit_sec(m, options.maxtime)

    @variable(m, 𝑧̂[𝒫])
    @variable(m, a[ℐₚ, 𝒦])
    @variable(m, b[𝒦])

    @variable(m, 𝑢[𝒫, 𝒦], Bin)

    if options.metric == :l2
        @objective(m, Min, sum((z[p] - 𝑧̂[p])^2 for p ∈ 𝒫))
    elseif options.metric == :max
        𝑡 = m[:𝑡] = @variable(m)
        @objective(m, Min, 𝑡)
        for p ∈ 𝒫
            @constraint(m, 𝑡 ≥ (z[p] - 𝑧̂[p]))
            @constraint(m, 𝑡 ≥ (𝑧̂[p] - z[p]))
        end
    elseif options.metric == :l1
        𝑡 = m[:𝑡] = @variable(m, [𝒫])
        @objective(m, Min, sum(𝑡))
        for p ∈ 𝒫
            @constraint(m, 𝑡[p] ≥ (z[p] - 𝑧̂[p]))
            @constraint(m, 𝑡[p] ≥ (𝑧̂[p] - z[p]))
        end
    end

    for p ∈ 𝒫, k ∈ 𝒦
        @constraint(m, 𝑧̂[p] ≥ sum(a[j, k] * x[p][j] for j ∈ ℐₚ) + b[k])
        @constraint(
            m,
            𝑧̂[p] ≤
            sum(a[j, k] * x[p][j] for j ∈ ℐₚ) + b[k] + Mᵖ[p] * (1 - 𝑢[p, k])
        )
    end

    if isnothing(start)
        # Workaround for infeasibility (really unbounded?) when options.strict == :none.
        # The bound is relative to the range of z to be invariant to shifts in z.
        𝑧ᵐᵃˣ = maximum(z) + (maximum(z) - minimum(z))
        for p ∈ 𝒫
            @constraint(m, 𝑧̂[p] <= 𝑧ᵐᵃˣ)
        end
    else
        # An optimal solution deviates at most r from the data in each point
        for p ∈ 𝒫
            set_lower_bound(𝑧̂[p], z[p] - start.r)
            set_upper_bound(𝑧̂[p], z[p] + start.r)
        end
    end

    if options.strict == :outer
        # All planes below the data, i.e. the approximation is below the data
        for p ∈ 𝒫, k ∈ 𝒦
            @constraint(m, z[p] ≥ sum(a[j, k] * x[p][j] for j ∈ ℐₚ) + b[k])
        end
    elseif options.strict == :inner
        # The approximation is above the data
        for p ∈ 𝒫
            @constraint(m, 𝑧̂[p] ≥ z[p])
        end
    end

    for p ∈ 𝒫
        @constraint(m, sum(𝑢[p, k] for k ∈ 𝒦) == 1)
    end

    if D == 1
        # For one-dimensional approximation, keep planes sorted by gradient
        for k ∈ 𝒦, j ∈ ℐₚ
            if k > 1
                @constraint(m, a[j, k-1] ≤ a[j, k])
            end
        end
    end

    if !isnothing(start)
        _set_milp_start(m, start, options.metric)
    end

    set_optimizer(m, options.optimizer)
    optimize!(m)

    status = termination_status(m)
    if status ∉ (MOI.OPTIMAL, MOI.TIME_LIMIT)
        error("Optimization failed $status\n$(raw_status(m))")
    end
    if primal_status(m) != MOI.FEASIBLE_POINT
        error("No feasible solution found within the time limit")
    end
    if status == MOI.TIME_LIMIT
        @warn "MILP reached the time limit, the approximation may not be optimal" maxlog =
            1
    end

    aᴼᵖᵗ = value.(a)
    bᴼᵖᵗ = value.(b)

    return PWAFunc{Convex,D}([Plane(Tuple(aᴼᵖᵗ.data[:, k]), bᴼᵖᵗ[k]) for k ∈ 𝒦])
end

function _validate(options::MILP)
    valid(name, value, values) =
        value ∈ values || throw(
            ArgumentError(
                "Unsupported $name = :$value for MILP, use one of " *
                join(":" .* string.(values), ", "),
            ),
        )
    valid("metric", options.metric, (:l1, :l2, :max))
    valid("strict", options.strict, (:none, :inner, :outer))
    options.planes ≥ 1 || throw(ArgumentError("MILP requires planes ≥ 1"))
    options.maxtime > 0 || throw(ArgumentError("MILP requires maxtime > 0"))
    return options
end

# Heuristic start for the MILP from a least squares Cluster fit, adjusted to
# satisfy options.strict. Returns the plane coefficients a and b, the values of
# the approximation ẑ, the active plane in each point, the distance gap from the
# approximation down to the lowest plane in each point, the deviations e = ẑ - z,
# and the bound r on the deviation from the data in each point for an optimal
# solution.
function _milp_start(input::FunctionEvaluations{D}, options::MILP) where {D}
    x = input.points
    z = input.values
    K = options.planes
    X = [p[i] for i ∈ 1:D, p ∈ x]
    heuristic = Cluster(; planes = K, rng = Random.Xoshiro(42))
    planes = _convex_linearization_mb(X, z, heuristic).planes
    if D == 1
        sort!(planes; by = p -> p.α[1])
    end
    # Fill up with copies of the last plane if the heuristic found fewer planes
    while length(planes) < K
        push!(planes, planes[end])
    end
    a = [planes[k].α[j] for j ∈ 1:D, k ∈ 1:K]
    b = [planes[k].β for k ∈ 1:K]

    planevalues() = [
        sum(a[j, k] * x[p][j] for j ∈ 1:D) + b[k] for p ∈ eachindex(x), k ∈ 1:K
    ]
    if options.strict == :outer
        b .-= max.(0, vec(maximum(planevalues() .- z; dims = 1)))
    elseif options.strict == :inner
        b .+= max(0, maximum(z .- vec(maximum(planevalues(); dims = 2))))
    end

    v = planevalues()
    ẑ = vec(maximum(v; dims = 2))
    active = [argmax(v[p, :]) for p ∈ eachindex(x)]
    gap = ẑ .- vec(minimum(v; dims = 2))

    e = ẑ .- z
    r = if options.metric == :l1
        sum(abs, e)
    elseif options.metric == :l2
        sqrt(sum(abs2, e))
    else
        maximum(abs, e)
    end
    # Small margin to keep the start feasible despite rounding
    r = r * (1 + 1e-6) + 1e-9 * (1 + maximum(abs, z))

    return (; a, b, ẑ, active, gap, e, r)
end

function _set_milp_start(m, start, metric)
    a, b, 𝑧̂, 𝑢 = m[:a], m[:b], m[:𝑧̂], m[:𝑢]
    for k ∈ axes(start.a, 2)
        for j ∈ axes(start.a, 1)
            set_start_value(a[j, k], start.a[j, k])
        end
        set_start_value(b[k], start.b[k])
    end
    for (p, ẑₚ) ∈ enumerate(start.ẑ)
        set_start_value(𝑧̂[p], ẑₚ)
        for k ∈ axes(start.a, 2)
            set_start_value(𝑢[p, k], k == start.active[p] ? 1.0 : 0.0)
        end
    end
    if metric == :max
        set_start_value(m[:𝑡], maximum(abs, start.e))
    elseif metric == :l1
        for (p, eₚ) ∈ enumerate(start.e)
            set_start_value(m[:𝑡][p], abs(eₚ))
        end
    end
    return
end

# Sample the function on a uniform grid within the given bounding box using nsamples in each dimension
function _sample_uniform(f::Function, bbox::Vector{<:Tuple}, nsamples)
    dims = length(bbox)
    if dims == 1
        it = LinRange(bbox[dims][1], bbox[dims][2], nsamples)
        x = Tuple.(collect(it))
    else
        it = Iterators.product(
            (LinRange(bbox[d][1], bbox[d][2], nsamples) for d ∈ 1:dims)...,
        )
        x = vec(collect(it))
    end
    y = [f(xx) for xx ∈ x]
    return FunctionEvaluations(x, y)
end

"""
    approx(f::Function, bbox::Vector{<:Tuple}, c::Curvature, a::Algorithm;  kwargs...)

Approximate the function using a uniform sampling over the bounding box `bbox`

Additional keyword arguments:
- `nsample`: the number of points used in each dimension (default = 10)
"""
function approx(f::Function, bbox::Vector{<:Tuple}, c::Curvature, a::Algorithm;)
    samples = 3 * a.planes

    return approx(_sample_uniform(f, bbox, samples), c, a)
end

# Utility function to find the value of a hyperplane at the point x
function _plane_f(x, normal, d)
    return abs(normal[3]) > 1e-4 ?
           (-1 / normal[3]) * (normal[1] * x[1] + normal[2] * x[2] - d) : 0.0
end

function _sign(p1, p2, p3)
    return (p1[1] - p3[1]) * (p2[2] - p3[2]) - (p2[1] - p3[1]) * (p1[2] - p3[2])
end

# Returns true if the point pt is lying inside
# the triangle with vertices v1, v2 and v3
function _point_in_triangle(pt, v1, v2, v3)
    d1 = _sign(pt, v1, v2)
    d2 = _sign(pt, v2, v3)
    d3 = _sign(pt, v3, v1)

    has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0)
    has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0)

    return !(has_neg && has_pos)
end

# Measure on the elongation of the triangle defined
# by v2, v2 and v3
function _sliver_parameter(v1, v2, v3)
    area =
        0.5 * abs(
            v1[1] * (v2[2] - v3[2]) +
            v2[1] * (v3[2] - v1[2]) +
            v3[1] * (v1[2] - v2[2]),
        )
    perimeter = norm(v1 .- v2) + norm(v1 .- v3) + norm(v2 .- v3)
    return 2 * area / perimeter
end

_linear_big_M(input::FunctionEvaluations) = maximum(_linear_big_M_points(input))

# Big-M for each data point: the largest distance from the point down to the
# candidate planes
function _linear_big_M_points(input::FunctionEvaluations{D}) where {D}
    Mᵖ = zeros(length(input))
    if D == 1
        # For one dimensional problems calculate a big M by
        # checking all segments not having another data point in the interior
        x = [p[1] for p ∈ input.points]
        z = input.values
        x₊ = collect(zip(x, z))

        # Consider points that are closer than 1/100 of the average distance
        # to be so close that they are not considered as a segment
        δ = (maximum(x) - minimum(x)) / (100 * length(x))

        segments = collect(combinations(x₊, 2))
        for s ∈ segments
            x1, z1 = s[1]
            x2, z2 = s[2]
            # Avoid points that are too close
            if abs(x1 - x2) > δ
                # Check for datapoints in the interior
                empty = (
                    count(p -> p > min(x1, x2) + δ && p < max(x1, x2) - δ, x) == 0
                )
                if empty
                    a = (z2 - z1) / (x2 - x1)
                    b = z1 - a * x1
                    Mᵖ .= max.(Mᵖ, z .- a .* x .- b)
                end
            end
        end
    elseif D == 2
        # For two dimensions, chcck all planes that
        # is tangential to a triangle formed by three
        # of the points, avoiding triangles that have data points
        # in the interior or those that are very elongated.
        # The points are scaled to the unit square to make the elongation
        # and vertical plane thresholds independent of the scale of the data.
        # The scaling is affine, so the values of the planes at the points
        # are unchanged. The elongation threshold of 0.025 on the unit square
        # corresponds to 0.05 on [-1, 1]².
        lo = [minimum(p[i] for p ∈ input.points) for i ∈ 1:2]
        hi = [maximum(p[i] for p ∈ input.points) for i ∈ 1:2]
        w = [hi[i] > lo[i] ? hi[i] - lo[i] : 1.0 for i ∈ 1:2]
        x = [Tuple((p[i] - lo[i]) / w[i] for i ∈ 1:2) for p ∈ input.points]
        z = input.values
        x₊ = [(x[i]..., z[i]) for i ∈ 1:length(x)]
        triangles = collect(combinations(x₊, 3))
        for t ∈ triangles
            v1 = t[1][1:2]
            v2 = t[2][1:2]
            v3 = t[3][1:2]
            empty = true
            for pt ∈ x
                if !(pt in [v1, v2, v3]) && _point_in_triangle(pt, v1, v2, v3)
                    empty = false
                end
            end

            if empty && _sliver_parameter(v1, v2, v3) > 0.025
                normal = cross(collect(t[1] .- t[2]), collect(t[1] .- t[3]))
                d = dot(normal, t[1])
                Mᵖ .= max.(Mᵖ, [z[i] - _plane_f(x[i], normal, d) for i ∈ 1:length(x)])
            end
        end
    else
        error("linear big M only supported for 1D and 2D")
    end

    # Increase slightly to avoid rounding errors
    return 1.1 .* Mᵖ
end
