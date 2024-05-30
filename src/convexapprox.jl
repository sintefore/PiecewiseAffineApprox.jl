defaultpenalty() = :l1
defaultplanes() = 4
defaulttimelimit() = 60

"""
    approx(input::FunctionEvaluations{D}, c::Curvature, a::Algorithm)

Return PWLFunc{Convex,D} or PWLFunc{Concave,D} depending on `c`, approximating the `input` points in `D` dimensions
"""
function approx(input, c::Concave, a::Algorithm)
    cv = approx(FunctionEvaluations(input.points, -input.values), Convex(), a;)
    return PWLFunc{Concave,dims(cv)}(cv.planes)
end

dims(pwa::PWLFunc{C,D}) where {C,D} = D

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

"""
    approx(input::FunctionEvaluations{D}, c::Convex, a::Heuristic; kwargs...) where D

Approximate using a heuristic that works for general dimensions.

Additional keyword arguments:
- `trials`: number of restarts (default = 20)
- `itlim`: max refining iterations in each trial (default = 50),
"""
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Heuristic;
) where {D}
    x = [p[i] for i ∈ 1:D, p ∈ input.points]
    z = input.values
    return _convex_linearization_mb(x, z, a)
end

# Optimal convex approximation using mixed integer optimization
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    options::Optimized,
) where {D}
    𝒫 = input.points
    z = input.values
    zᵖ = Dict(zip(𝒫, z))
    𝒦 = 1:options.planes
    ℐₚ = 1:D

    Mᵇⁱᵍ = _linear_big_M(input)

    m = Model()

    set_time_limit_sec(m, options.maxtime)

    @variable(m, 𝑧̂[𝒫])
    @variable(m, a[ℐₚ, 𝒦])
    @variable(m, b[𝒦])

    @variable(m, 𝑢[𝒫, 𝒦], Bin)

    if options.pen == :l2
        @objective(m, Min, sum((zᵖ[p] - 𝑧̂[p])^2 for p ∈ 𝒫))
    elseif options.pen == :max
        𝑡 = @variable(m)
        @objective(m, Min, 𝑡)
        for p ∈ 𝒫
            @constraint(m, 𝑡 ≥ (zᵖ[p] - 𝑧̂[p]))
            @constraint(m, 𝑡 ≥ (𝑧̂[p] - zᵖ[p]))
        end
    elseif options.pen == :l1
        𝑡 = @variable(m, [𝒫])
        @objective(m, Min, sum(𝑡))
        for p ∈ 𝒫
            @constraint(m, 𝑡[p] ≥ (zᵖ[p] - 𝑧̂[p]))
            @constraint(m, 𝑡[p] ≥ (𝑧̂[p] - zᵖ[p]))
        end
    else
        error("Unrecognized/unsupported penalty type $(options.pen)")
    end

    for p ∈ 𝒫, k ∈ 𝒦
        @constraint(m, 𝑧̂[p] ≥ sum(a[j, k] * p[j] for j ∈ ℐₚ) + b[k])
        @constraint(
            m,
            𝑧̂[p] ≤
            sum(a[j, k] * p[j] for j ∈ ℐₚ) + b[k] + 2 * Mᵇⁱᵍ * (1 - 𝑢[p, k])
        )
    end

    # Workaround for infeasibility (really unbounded?) when options.strict == :none
    for p ∈ 𝒫
        @constraint(m, 𝑧̂[p] <= maximum(z) * 1.1)
    end

    if options.strict == :outer
        for p ∈ 𝒫, k ∈ 𝒦
            @constraint(m, zᵖ[p] ≥ sum(a[j, k] * p[j] for j ∈ ℐₚ) + b[k])
        end
    elseif options.strict == :inner
        for p ∈ 𝒫, k ∈ 𝒦
            @constraint(m, zᵖ[p] ≤ sum(a[j, k] * p[j] for j ∈ ℐₚ) + b[k])
        end
    end

    for p ∈ 𝒫
        @constraint(m, sum(𝑢[p, k] for k ∈ 𝒦) ≥ 1)
    end

    if D == 1
        # For one-dimensional approximation, keep planes sorted by gradient
        for k ∈ 𝒦, j ∈ ℐₚ
            if k > 1
                @constraint(m, a[j, k-1] ≤ a[j, k])
            end
        end
    end

    set_optimizer(m, options.optimizer)
    optimize!(m)

    if termination_status(m) ∉ [MOI.OPTIMAL, MOI.TIME_LIMIT]
        error("Optimization failed $(termination_status(m))\n$(raw_status(m))")
    end

    aᴼᵖᵗ = value.(a)
    bᴼᵖᵗ = value.(b)

    return PWLFunc{Convex,D}([Plane(Tuple(aᴼᵖᵗ.data[:, k]), bᴼᵖᵗ[k]) for k ∈ 𝒦])
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

function _linear_big_M(input::FunctionEvaluations{D}) where {D}
    Mᵇⁱᵍ = 0
    if D == 1
        # For one dimensional problems calculate a big M by
        # checking all segments not having another data point in the interior
        x = [p[1] for p in input.points]
        z = input.values
        x₊ = collect(zip(x, z))

        # Consider points that are closer than 1/100 of the average distance
        # to be so close that they are not considered as a segment
        δ = (maximum(x) - minimum(x)) / (100 * length(x))

        segments = collect(combinations(x₊, 2))
        for s in segments
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
                    M = maximum(z .- a .* x .- b)
                    Mᵇⁱᵍ = max(M, Mᵇⁱᵍ)
                end
            end
        end
    elseif D == 2
        # For two dimensions, chcck all planes that
        # is tangential to a triangle formed by three
        # of the points, avoiding triangles that have data points
        # in the interior or those that are very elongated.
        x = input.points
        z = input.values
        x₊ = [(x[i]..., z[i]) for i ∈ 1:length(x)]
        triangles = collect(combinations(x₊, 3))
        for t in triangles
            v1 = t[1][1:2]
            v2 = t[2][1:2]
            v3 = t[3][1:2]
            empty = true
            for pt in x
                if !(pt in [v1, v2, v3]) && _point_in_triangle(pt, v1, v2, v3)
                    empty = false
                end
            end

            if empty && _sliver_parameter(v1, v2, v3) > 0.05
                normal = cross(collect(t[1] .- t[2]), collect(t[1] .- t[3]))
                d = dot(normal, t[1])
                M = maximum(
                    z[i] - _plane_f(x[i], normal, d) for i ∈ 1:length(x)
                )
                Mᵇⁱᵍ = max(M, Mᵇⁱᵍ)
            end
        end
    else
        error("linear big M only supported for 1D and 2D")
    end

    # Increase slightly to avoid rounding errors
    return 1.1 * Mᵇⁱᵍ
end
