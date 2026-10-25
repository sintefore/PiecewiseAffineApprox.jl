defaultmetric() = :l1
defaultplanes() = 4
defaulttimelimit() = 60

"""
    approx(input::FunctionEvaluations{D}, c::Curvature, a::ApproximationAlgorithm)

Return a `PWAFunc{Convex,D}` or `PWAFunc{Concave,D}` depending on `c`,
approximating the `input` points in `D` dimensions with the algorithm `a`, or
a [`DCFunc`](@ref) for `c = DC()`. Concave data is approximated by negating
the values, approximating with a convex function, and negating the result.
"""
function approx(
    input::FunctionEvaluations{D},
    c::Convex,
    a::ApproximationAlgorithm,
) where {D}
    if D ∉ supported_dimensions(a)
        throw(
            ArgumentError(
                "$(nameof(typeof(a))) supports functions of " *
                "$(_dimensions_text(supported_dimensions(a))), got $D variables. " *
                "Use another algorithm, e.g. Cluster.",
            ),
        )
    end
    isempty(input.points) && throw(ArgumentError("The data has no points"))
    return _fit(a, input)
end

function approx(
    input::FunctionEvaluations,
    c::Concave,
    a::ApproximationAlgorithm,
)
    cv = approx(FunctionEvaluations(input.points, -input.values), Convex(), a)
    return PWAFunc{Concave,dims(cv)}(cv.planes)
end

dims(pwa::PWAFunc{C,D}) where {C,D} = D

function _dimensions_text(r)
    return length(r) == 1 ? "$(first(r)) variable" :
           "$(first(r)) to $(last(r)) variables"
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
    y = _evaluate_points(f, x)
    return FunctionEvaluations(x, y)
end

# The values of f at the points, given as tuples: f is called with the tuple
# if it accepts one, and otherwise with the coordinates as separate arguments
function _evaluate_points(f::Function, x)
    isempty(x) && return Float64[]
    if applicable(f, first(x))
        return [f(xx) for xx ∈ x]
    end
    return [f(xx...) for xx ∈ x]
end

"""
    approx(f::Function, bbox::Vector{<:Tuple}, c::Curvature, a::ApproximationAlgorithm; nsamples)

Approximate the function using a uniform sampling over the bounding box `bbox`,
with one `(lower, upper)` tuple for each dimension. The function is called with
the point as a tuple, e.g. `x -> x[1]^2 + x[2]^2`, if it accepts one, and
otherwise with the coordinates as separate arguments, e.g. `sin` or
`(x, y) -> x^2 + y^2`.

Additional keyword arguments:
- `nsamples`: the number of points used in each dimension. The default depends
  on the algorithm and the number of dimensions D:
  - [`MILP`](@ref): about `3 * D * planes` points in total, i.e.
    `ceil((3 * D * planes)^(1 / D))` in each dimension, since the size of
    the optimization problem grows with the number of points,
  - [`Cluster`](@ref) and [`Interpol`](@ref): `3 * planes` in each dimension,
  - other algorithms: about 100 points in total.
"""
function approx(
    f::Function,
    bbox::Vector{<:Tuple},
    c::Curvature,
    a::ApproximationAlgorithm;
    nsamples = _default_samples(a, length(bbox)),
)
    nsamples ≥ 2 || throw(ArgumentError("nsamples must be at least 2"))
    return approx(_sample_uniform(f, bbox, nsamples), c, a)
end

# Default number of samples in each of D dimensions
function _default_samples(a::ApproximationAlgorithm, D)
    p = planes(a)
    return isnothing(p) ? max(3, round(Int, 100^(1 / D))) : 3p
end
_default_samples(a::MILP, D) = max(2, ceil(Int, (3 * D * a.planes)^(1 / D)))
