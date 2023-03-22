defaultpenalty() = :l1
defaultpenalty2D() = :l2
defaultplanes() = 4
defaultbigM() = :linear_big_M

"""
    approx(input::FunctionEvaluations{D}, c::Curvature, a::Algorithm; kwargs...)

Return ConvexPWLFunc{D} or ConcavePWLFunc{D} depending on `c`, approximating the `input` points in `D` dimensions

Accepted keyword arguments currently include:
- `optimizer`: JuMP Optimizer
- `planes`: number of (hyper)planes to use for approximation
- `strict`: (TODO: Better name?) `strict ∈ (:none, :over, :under)`
- `pen`:  the metric used to measure deviation `pen ∈ (:l1,:l2)`
- `show_res`: TODO:Remove this 
"""
function approx(input, c::Concave, a; kwargs...)
    cv = approx(
        FunctionEvaluations(input.points, -input.values),
        Convex(),
        a;
        kwargs...,
    )
    # TODO: Generalize on D (hard coded to 2 for now)
    return PWLFunc{Concave,dims(cv)}(cv.planes)
end

dims(pwl::PWLFunc{C,D}) where {C,D} = D

# Using dispatch for specializing on dimensions. If performance were a concern,
# maybe just do branching and call specialized function directly
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a;
    kwargs...,
) where {D}
    return approx(input, c, a, Val(D); kwargs...)
end

# Specialized for 1D
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Interpol,
    ::Val{1};
    kwargs...,
) where {D}
    defaults = (
        planes = defaultplanes(),
        pen = defaultpenalty2D(),
        strict = :none,
        show_res = false,
    )
    options = merge(defaults, kwargs)
    # Wrap for now, TODO: move here
    return _convex_linearization_ipol(
        [i[1] for i ∈ input.points],
        input.values,
        options.optimizer;
        kwargs...,
    )
end

function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Optimized,
    ::Val{1};
    kwargs...,
) where {D}
    defaults = (
        planes = defaultplanes(),
        pen = defaultpenalty2D(),
        bigM = defaultbigM(),
        strict = :none,
        show_res = false,
    )
    options = merge(defaults, kwargs)
    # Wrap until big M issue is solved generally
    # TODO: move here   
    return convex_linearization_fit(
        [i[1] for i ∈ input.points],
        input.values,
        options.optimizer;
        kwargs...,
    )
end

"""
    approx(input::FunctionEvaluations{D}, c::Convex, a::Heuristic; kwargs...) where D

Approximate using heuristic for general dimension
Additional keyword arguments:
- `trials`=20
- `itlim`=50,
"""
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Heuristic;
    kwargs...,
) where {D}
    x = [p[i] for i ∈ 1:D, p ∈ input.points]
    z = input.values
    return _convex_linearization_mb(x, z; kwargs...)
end

# General D
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::Optimized,
    dims;
    kwargs...,
) where {D}
    defaults = (
        planes = defaultplanes(),
        pen = defaultpenalty2D(),
        strict = :none,
        show_res = false,
    )
    options = merge(defaults, kwargs)

    𝒫 = input.points
    z = input.values
    zᵖ = Dict(zip(𝒫, z))
    𝒦 = 1:options.planes
    ℐₚ = 1:length(𝒫[1])

    Mᵇⁱᵍ = linear_big_M(𝒫, z)

    m = Model()
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
            sum(a[j, k] * p[j] for j ∈ ℐₚ) + b[k] + Mᵇⁱᵍ * (1 - 𝑢[p, k])
        )
    end

    # Workaround for infeasibility (really unbounded?) when options.strict == :none
    for p ∈ 𝒫
        @constraint(m, 𝑧̂[p] <= maximum(z) * 1.1)
    end

    if options.strict == :above
        for p ∈ 𝒫, k ∈ 𝒦
            @constraint(m, zᵖ[p] ≥ sum(a[j, k] * p[j] for j ∈ ℐₚ) + b[k])
        end
    elseif options.strict == :below
        for p ∈ 𝒫, k ∈ 𝒦
            @constraint(m, zᵖ[p] ≤ sum(a[j, k] * p[j] for j ∈ ℐₚ) + b[k])
        end
    end

    for p ∈ 𝒫
        @constraint(m, sum(𝑢[p, k] for k ∈ 𝒦) ≥ 1)
    end

    set_optimizer(m, options.optimizer)
    optimize!(m)

    if termination_status(m) ∉ [MOI.OPTIMAL, MOI.TIME_LIMIT]
        error("Optimization failed $(termination_status(m))\n$(raw_status(m))")
    end

    if options.show_res
        println("Optimize succeed for $(options.pen)")
        val = objective_value(m)
        println("Objective value = $val")
    end

    aᴼᵖᵗ = value.(a)
    bᴼᵖᵗ = value.(b)

    return PWLFunc{Convex,D}([Plane(Tuple(aᴼᵖᵗ.data[:, k]), bᴼᵖᵗ[k]) for k ∈ 𝒦])
end

# Sample the function on a uniform grid within the given bounding box using nsamples in each dimension
function sample_uniform(f::Function, bbox::Vector{<:Tuple}, nsamples)
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
- `nsample`=10
"""
function approx(
    f::Function,
    bbox::Vector{<:Tuple},
    c::Curvature,
    a::Algorithm;
    kwargs...,
)
    defaults = (nsample = 10, planes = defaultplanes())
    options = merge(defaults, kwargs)

    samples = max(options.nsample, 3 * options.planes)

    return approx(sample_uniform(f, bbox, samples), c, a; kwargs...)
end

#=
    Original interface starts here
=#

"""
    convex_linearization(x, z, optimizer; kwargs...)

Computes a piecewise linear function that approximates the measurements given by `x` and `z`.

# Arguments
- `method::Symbol:=fit`: the method used for approximation
- `dimensions::Integer:=2`: the number of dimensions of the function domain
- `planes::Integer=5`: the number of segments to use 
- `planes::Integer=4`: the number of planes to use in 2D PWL functions
- `strict::Symbol=:none`: defines it is a general approximation, or an overestimation or underestimation
- `pen::Symbol=:l1`: the metric used to measure deviation


"""
function convex_linearization(x, z, optimizer; kwargs...)
    defaults = (; method = :fit, dimensions = 1)
    options = merge(defaults, kwargs)

    method = options.method
    dimensions = options.dimensions

    if dimensions == 1
        @assert(length(x) == length(z))
        if method == :fit
            return convex_linearization_fit(x, z, optimizer; kwargs...)
        elseif method == :ipol
            return convex_linearization_ipol(x, z, optimizer; kwargs...)
        else
            error("Unrecognized method $method")
        end
    elseif dimensions >= 2
        if method == :fit
            return convex_ND_linearization_fit(x, z, optimizer; kwargs...)
        end
    else
        error("Unrecognized number of dimensions $dimensions")
    end
end
function conv_linear_big_M(x, z)
    N = length(x)
    cᵉˢᵗ = (z[N] - z[N-1]) / (x[N] - x[N-1])
    return 2 * cᵉˢᵗ * (last(x) - first(x)) - maximum(z)
end

function conv_linear_big_M_ND(x, z)
    return 2 * maximum(z)
end

#convex_linearization(x, z, optimizer; kwargs...)  = 
#    convex_linearization([xx for xx in x], [zz for zz in z], optimizer; kwargs...)

function convex_linearization_fit(x::Vector, z::Vector, optimizer; kwargs...)
    defaults = (
        planes = defaultplanes(),
        pen = defaultpenalty(),
        bigM = defaultbigM(),
        strict = false,
        start_origin = false,
        show_res = false,
    )
    options = merge(defaults, kwargs)

    N = length(x)
    𝒩 = 1:N
    𝒦 = 1:options.planes

    Mᵇⁱᵍ = big_M(x, z, options.bigM)

    m = Model()
    @variable(m, 𝑧̂[𝒩])
    @variable(m, 𝑐[𝒦])
    @variable(m, 𝑑[𝒦])
    @variable(m, 𝑢[𝒩, 𝒦], Bin)

    if options.pen == :l2
        @objective(m, Min, sum((z[i] - 𝑧̂[i])^2 for i ∈ 𝒩))
    elseif options.pen == :max
        @variable(m, 𝑡)
        @objective(m, Min, 𝑡)
        for i ∈ 𝒩
            @constraint(m, 𝑡 ≥ (z[i] - 𝑧̂[i]))
            @constraint(m, 𝑡 ≥ (𝑧̂[i] - z[i]))
        end
    elseif options.pen == :l1
        @variable(m, 𝑡[𝒩])
        @objective(m, Min, sum(𝑡))
        for i ∈ 𝒩
            @constraint(m, 𝑡[i] ≥ (z[i] - 𝑧̂[i]))
            @constraint(m, 𝑡[i] ≥ (𝑧̂[i] - z[i]))
        end
    else
        error("Unrecognized/unsupported penalty type $(options.pen)")
    end

    for i ∈ 𝒩, k ∈ 𝒦
        @constraint(m, 𝑧̂[i] ≥ 𝑐[k] * x[i] + 𝑑[k])

        @constraint(m, 𝑧̂[i] ≤ 𝑐[k] * x[i] + 𝑑[k] + Mᵇⁱᵍ * (1 - 𝑢[i, k]))
    end

    if options.strict
        for i ∈ 𝒩, k ∈ 𝒦
            @constraint(m, z[i] ≥ 𝑐[k] * x[i] + 𝑑[k])
        end
    end

    if options.start_origin
        @constraint(m, 𝑑[1] == 0.0)
    end

    for i ∈ 𝒩
        @constraint(m, sum(𝑢[i, k] for k ∈ 𝒦) ≥ 1)
    end

    for k ∈ 𝒦
        if k > 1
            @constraint(m, 𝑐[k-1] ≤ 𝑐[k])
        end
    end

    set_optimizer(m, optimizer)
    optimize!(m)

    if termination_status(m) ∉ [MOI.OPTIMAL, MOI.TIME_LIMIT]
        error("Optimization failed $(termination_status(m))")
    end

    if options.show_res
        println("Optimize succeed for $(options.pen)")
        val = objective_value(m)
        println("Objective value = $val")
    end

    𝑐ᴼᵖᵗ = value.(𝑐)
    𝑑ᴼᵖᵗ = value.(𝑑)

    return PWLFunc{Convex,1}([Plane(Tuple(𝑐ᴼᵖᵗ[k]), 𝑑ᴼᵖᵗ[k]) for k ∈ 𝒦])
    # return ConvexPWLFunction([𝑐ᴼᵖᵗ[k] for k ∈ 𝒦], [𝑑ᴼᵖᵗ[k] for k ∈ 𝒦], minimum(x), maximum(x))
end

function convex_linearization(f::Function, xmin, xmax, optimizer; kwargs...)
    @assert(xmin < xmax)

    defaults = (nsample = 10, planes = defaultplanes())
    options = merge(defaults, kwargs)

    samples = max(options.nsample, 3 * options.planes)

    step = (xmax - xmin) / samples
    x = [i for i ∈ xmin:step:xmax]
    y = [f(xx) for xx ∈ x]
    return convex_linearization(x, y, optimizer; kwargs...)
end

function concave_linearization(x, z, optimizer; kwargs...)
    if (size(x, 2) < 2) && (length(first(x)) < 2)
        return ConcavePWLFunction(
            convex_linearization(x, -z, optimizer; kwargs...),
        )
    else
        return ConcavePWLFunctionND(
            convex_linearization(x, -z, optimizer; kwargs...),
        )
    end
end

function concave_linearization(f::Function, xmin, xmax, optimizer; kwargs...)
    if (size(x, 2) < 2) && (length(first(X)) < 2)
        return ConcavePWLFunction(
            convex_linearization(x -> -f(x), xmin, xmax, optimizer; kwargs...),
        )
    else
        return ConcavePWLFunctionND(
            convex_linearization(x -> -f(x), xmin, xmax, optimizer; kwargs...),
        )
    end
end

function mat2tuples(x::Matrix)
    return collect(Tuple(x'[:, i]) for i ∈ 1:size(x', 2))
end
# mat2tuples2(x::Matrix) = ((Tuple(x[i,:]) for i in 1:size(x,1)))

function tuples2mat(𝒫::Vector{Tuple{Float64,Float64}})
    return reduce(hcat, getindex.(𝒫, i) for i ∈ eachindex(𝒫[1]))
end

function convex_ND_linearization_fit(
    x::Matrix{Float64},
    z,
    optimizer;
    kwargs...,
)
    return convex_ND_linearization_fit(mat2tuples(x), z, optimizer; kwargs...)
end

function big_M(x, z, bigM)
    if bigM == :linear_big_M
        return linear_big_M(x, z)
    elseif bigM == :conv_linear_big_M_ND
        return conv_linear_big_M_ND(x, z)
    elseif bigM == :conv_linear_big_M
        return conv_linear_big_M(x, z)
    else
        error(
            "Big-M function not available. Options are :linear_big_M, :conv_linear_big_M_ND, and :conv_linear_big_M",
        )
    end
end

function linear_big_M(x, z)
    dims = length(x[1])

    ## calculate extremum points
    xₗ = if dims > 2
        minimum.(x)
    else
        minimum(x)
    end
    xᵤ = if dims > 2
        maximum.(x)
    else
        maximum(x)
    end

    x₊ = [Tuple([x[i]..., z[i]]) for i ∈ 1:length(x)]

    if dims <= 2 # works only for up to 3 dimensions for now
        ℋ = collect(combinations(x₊, dims + 1))
        Mᵇⁱᵍ = zeros(length(ℋ))

        k = 1
        for ℎ ∈ ℋ
            if dims == 1
                a = (ℎ[2][2] - ℎ[1][2]) / (ℎ[2][1] - ℎ[1][1])
                b = ℎ[1][2] - a * ℎ[1][1]

                line_f(x, a, b) = a .* x + b
                # extrapolation for corner points
                extremum_values = [line_f(xₗ, a, b), line_f(xᵤ, a, b)]
                Mᵇⁱᵍ[k] =
                    abs(maximum(extremum_values) - minimum(extremum_values))
            elseif dims == 2
                normal = cross(collect(ℎ[1] .- ℎ[2]), collect(ℎ[1] .- ℎ[3]))
                d = dot(normal, ℎ[1])
                if sum(abs.(normal) .> 1e-4) == dims + 1
                    function plane_f(x, normal, d)
                        return abs(normal[3]) > 1e-4 ?
                               (-1 / normal[3]) * normal[1] .* x[1] +
                               normal[2] .* x[2] +
                               d : 0.0
                    end

                    # extrapolation for corner points
                    extremum_values = [
                        plane_f((xₗ[1], xₗ[2]), normal, d),
                        plane_f((xₗ[1], xᵤ[2]), normal, d),
                        plane_f((xᵤ[1], xₗ[2]), normal, d),
                        plane_f((xᵤ[1], xᵤ[2]), normal, d),
                    ]
                    Mᵇⁱᵍ[k] =
                        abs(maximum(extremum_values) - minimum(extremum_values))
                end
            else
                error("Big-M not defined for dimensions > 3")
            end

            k = k + 1
        end
        return maximum(Mᵇⁱᵍ)
    else
        error("Big-M calculation works only for 2 and 3 dimensions.")
    end
end
