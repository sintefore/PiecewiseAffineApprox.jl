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
    # Solve on data normalized to the unit box, to make the solver tolerances,
    # big-M values and margins independent of the scale of the data
    lo = ntuple(i -> minimum(p[i] for p ∈ input.points), D)
    w = ntuple(D) do i
        hi = maximum(p[i] for p ∈ input.points)
        return hi > lo[i] ? hi - lo[i] : 1.0
    end
    zᵐⁱⁿ, zᵐᵃˣ = extrema(input.values)
    zʷ = zᵐᵃˣ > zᵐⁱⁿ ? zᵐᵃˣ - zᵐⁱⁿ : 1.0
    normalized = FunctionEvaluations(
        [ntuple(i -> (p[i] - lo[i]) / w[i], D) for p ∈ input.points],
        (input.values .- zᵐⁱⁿ) ./ zʷ,
    )

    a, b = _milp(normalized, options)

    # Transform the planes back: z = zᵐⁱⁿ + zʷ (Σⱼ aⱼ (xⱼ - loⱼ) / wⱼ + b)
    return PWAFunc{Convex,D}([
        Plane(
            ntuple(j -> zʷ * a[j, k] / w[j], D),
            zᵐⁱⁿ + zʷ * (b[k] - sum(a[j, k] * lo[j] / w[j] for j ∈ 1:D)),
        ) for k ∈ axes(a, 2)
    ])
end

# Maximum number of times the big-M values are doubled when they may cut off
# the optimal solution
const _MILP_BIG_M_ROUNDS = 3

# Solve the MILP for the (normalized) input, returning the plane coefficients
# as a D × K matrix a and a vector b. The initial big-M values Mᵖ and the
# maximum number of times they are doubled can be given for testing.
function _milp(
    input::FunctionEvaluations{D},
    options::MILP;
    Mᵖ = nothing,
    rounds = _MILP_BIG_M_ROUNDS,
) where {D}
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
    Mᵖ = if !isnothing(Mᵖ)
        copy(Mᵖ)
    elseif D ≤ 2
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
    end
    @constraint(
        m,
        bigM[p ∈ 𝒫, k ∈ 𝒦],
        𝑧̂[p] ≤
        sum(a[j, k] * x[p][j] for j ∈ ℐₚ) + b[k] + Mᵖ[p] * (1 - 𝑢[p, k])
    )

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

    # The big-M values are heuristic. If a big-M constraint is tight for a
    # plane that is not active, the big-M may cut off better solutions. Then
    # the big-M values are doubled and the problem is solved again, as long as
    # this improves the objective.
    t₀ = time()
    objᵖʳᵉᵛ = Inf
    for iteration ∈ 0:rounds
        set_time_limit_sec(m, max(options.maxtime - (time() - t₀), 0.1))
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
            break
        end

        # Stop when the last increase of the big-M values did not improve the
        # objective (beyond the default relative MIP gap)
        obj = objective_value(m)
        iteration == 0 || obj < objᵖʳᵉᵛ - 1e-4 * abs(objᵖʳᵉᵛ) - 1e-9 || break
        objᵖʳᵉᵛ = obj

        tight = _tight_big_M_points(m, x, Mᵖ)
        isempty(tight) && break
        if iteration == rounds || time() - t₀ ≥ options.maxtime
            @warn "The big-M values in MILP may cut off the optimal solution"
            break
        end

        # Increase the big-M values for all points, since a cut off solution
        # usually needs larger values also in the points that are not (yet)
        # tight. Start from the current solution, which is kept first since it
        # is invalidated by modifying the model.
        solution = Dict(v => value(v) for v ∈ all_variables(m))
        for p ∈ 𝒫
            Mᵖ[p] = max(2 * Mᵖ[p], 1.0)
            for k ∈ 𝒦
                set_normalized_coefficient(bigM[p, k], 𝑢[p, k], Mᵖ[p])
                set_normalized_rhs(bigM[p, k], Mᵖ[p])
            end
        end
        for (v, val) ∈ solution
            set_start_value(v, val)
        end
    end

    return [value(a[j, k]) for j ∈ ℐₚ, k ∈ 𝒦], [value(b[k]) for k ∈ 𝒦]
end

# Points where a big-M constraint is tight for a plane that is not active
function _tight_big_M_points(m, x, Mᵖ; tol = 1e-6)
    a, b, 𝑧̂, 𝑢 = value.(m[:a]), value.(m[:b]), value.(m[:𝑧̂]), value.(m[:𝑢])
    D, K = size(a)
    return [
        p for p ∈ eachindex(x) if any(
            𝑢[p, k] < 0.5 &&
            𝑧̂[p] - sum(a[j, k] * x[p][j] for j ∈ 1:D) - b[k] ≥ Mᵖ[p] - tol for
            k ∈ 1:K
        )
    ]
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
    Mᵖ = if D == 1
        _linear_big_M_points_1d(input)
    elseif D == 2
        _linear_big_M_points_2d(input)
    else
        error("linear big M only supported for 1D and 2D")
    end

    # Increase slightly to avoid rounding errors
    return 1.1 .* Mᵖ
end

# For one dimensional problems the candidates are the lines through all
# segments not having another data point in the interior
function _linear_big_M_points_1d(input)
    x = [p[1] for p ∈ input.points]
    z = input.values
    N = length(x)
    Mᵖ = zeros(N)

    # Consider points that are closer than 1/100 of the average distance
    # to be so close that they are not considered as a segment
    δ = (maximum(x) - minimum(x)) / (100 * N)

    # With the points sorted, the segments from point i without a point in the
    # interior (more than δ from the end points) go to the points j with
    # x[i] + δ < x[j] ≤ x[j₀] + δ, where j₀ is the first point beyond x[i] + δ
    order = sortperm(x)
    xˢ = x[order]
    j₀ = 1
    for i ∈ 1:N
        while j₀ ≤ N && !(xˢ[j₀] > xˢ[i] + δ)
            j₀ += 1
        end
        j = j₀
        while j ≤ N && !(xˢ[j₀] < xˢ[j] - δ)
            if xˢ[j] - xˢ[i] > δ
                x1, z1 = xˢ[i], z[order[i]]
                x2, z2 = xˢ[j], z[order[j]]
                a = (z2 - z1) / (x2 - x1)
                b = z1 - a * x1
                for p ∈ 1:N
                    Mᵖ[p] = max(Mᵖ[p], z[p] - a * x[p] - b)
                end
            end
            j += 1
        end
    end
    return Mᵖ
end

# For two dimensions, the candidates are the planes through all triangles
# formed by three of the points, avoiding triangles that have data points in
# the interior or on the boundary, and those that are very elongated.
# The points are scaled to the unit square to make the elongation and vertical
# plane thresholds independent of the scale of the data. The scaling is affine,
# so the values of the planes at the points are unchanged. The elongation
# threshold of 0.025 on the unit square corresponds to 0.05 on [-1, 1]².
function _linear_big_M_points_2d(input)
    lo = [minimum(p[i] for p ∈ input.points) for i ∈ 1:2]
    hi = [maximum(p[i] for p ∈ input.points) for i ∈ 1:2]
    w = [hi[i] > lo[i] ? hi[i] - lo[i] : 1.0 for i ∈ 1:2]
    x = [((p[1] - lo[1]) / w[1], (p[2] - lo[2]) / w[2]) for p ∈ input.points]
    z = input.values
    N = length(x)
    Mᵖ = zeros(N)

    # Whether there is a point, other than the vertices, in the triangle
    function nonempty(v1, v2, v3)
        xᵐⁱⁿ, xᵐᵃˣ = extrema((v1[1], v2[1], v3[1]))
        yᵐⁱⁿ, yᵐᵃˣ = extrema((v1[2], v2[2], v3[2]))
        for pt ∈ x
            xᵐⁱⁿ ≤ pt[1] ≤ xᵐᵃˣ && yᵐⁱⁿ ≤ pt[2] ≤ yᵐᵃˣ || continue
            (pt == v1 || pt == v2 || pt == v3) && continue
            _point_in_triangle(pt, v1, v2, v3) && return true
        end
        return false
    end

    for i ∈ 1:N, j ∈ i+1:N, k ∈ j+1:N
        v1, v2, v3 = x[i], x[j], x[k]
        _sliver_parameter(v1, v2, v3) > 0.025 || continue
        nonempty(v1, v2, v3) && continue
        t1, t2, t3 = (v1..., z[i]), (v2..., z[j]), (v3..., z[k])
        normal = cross(collect(t1 .- t2), collect(t1 .- t3))
        d = dot(normal, collect(t1))
        for p ∈ 1:N
            Mᵖ[p] = max(Mᵖ[p], z[p] - _plane_f(x[p], normal, d))
        end
    end
    return Mᵖ
end
