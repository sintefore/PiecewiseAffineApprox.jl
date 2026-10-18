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
