# Input types
abstract type Curvature end
struct Concave <: Curvature end
struct Convex <: Curvature end

abstract type Algorithm end
"""
    Cluster
Compute affine approximation using the method proposed by Mangani & Boyd.

Note that this algoritm computes multiple approximations and selects the best.
If julia is started with multiple threads, these are computed in parallel. Consider
 how many threads will be beneficial, particularly when using a commercial solver where
 the license may restrict the number of simultanous solver instances.
"""
@kwdef struct Cluster{T} <: Algorithm
    planes::Int = defaultplanes()
    metric::Symbol = defaultmetric()
    trials::Int = 20
    itlim::Int = 50
    strict::Symbol = :none
    optimizer::T
end

"""
    Interpol

Compute affine approximation by method proposed by Flatberg. Only available for 1D.

The points do not need to be sorted, and points with the same x-value are merged
using the average of their function values. The number of `planes` must be less
than the number of distinct points.

The knots of the interpolant are chosen among the data points to minimize the
deviation from the data, measured by `metric`: the sum of absolute deviations
(`:l1`), the sum of squared deviations (`:l2`) or the maximum absolute deviation
(`:max`). Convexity of the interpolant is ensured by `convexify`:

- `:dp`: the best convex interpolant through (a subset of) the data points,
  found by dynamic programming. Exact, and does not need an `optimizer`, but
  may give a poor fit for data that is far from convex.
- `:data`: the data is first made convex by changing the function values as
  little as possible (measured by `metric`), and then interpolated. Handles
  non-convex and noisy data well. Requires an `optimizer`.
- `:knots`: the data is interpolated, and the function values at the knots
  are then made convex. Requires an `optimizer`.
- `:auto` (default): uses `:dp` if the data is convex or no `optimizer` is
  given. Otherwise both `:dp` and `:data` are computed, and the one with the
  smallest deviation from the data is returned.

The number of planes is given by `planes`, or selected automatically by
`selection`, with `planes` as the maximum number of planes considered (capped at
the number of distinct points minus one):

- `:none` (default): use `planes` planes.
- `:tolerance`: the fewest planes with an error of at most `tolerance`, where the
  error is the mean absolute deviation (`:l1`), the root mean square deviation
  (`:l2`) or the maximum absolute deviation (`:max`).
- `:aic`, `:bic`: minimize the Akaike or Bayesian information criterion,
  assuming normally (`:l2`) or Laplace (`:l1`) distributed errors, with two
  parameters per plane. Not available for `:max`.
- `:cv`: `folds`-fold cross validation, with the points assigned to folds by
  interleaving them in sorted order. Uses the one standard error rule: the
  fewest planes with an error on the held-out points within one standard error
  of the minimum.
"""
@kwdef struct Interpol{T} <: Algorithm
    planes::Int = defaultplanes()
    metric = defaultmetric()
    optimizer::T = nothing
    convexify::Symbol = :auto
    selection::Symbol = :none
    tolerance::Float64 = 0.0
    folds::Int = 5
end

"""
    MILP
Compute affine approximation using a variation of the method proposed by Toriello & Vielma.

Note that the resulting approximation is sensitive to the selection of the Big-M value used when
solving the optimization problem.
"""
@kwdef struct MILP{T} <: Algorithm
    planes::Int = defaultplanes()
    metric::Symbol = defaultmetric()
    strict::Symbol = :none
    maxtime::Int = defaulttimelimit()
    optimizer::T
end

"""
    Progressive

Compute affine approximation based on a variation of the method of Kazda and Li (2024)
specialized to convex approximations.
"""
@kwdef struct Progressive{T} <: Algorithm
    tolerance::Float64 = 0.01
    metric::Symbol = defaultmetric()
    optimizer::T
end

"""
    FullOrder

Compute affine approximation based on a variation of the method of Kazda and Li (2024)
specialized to convex approximations with full order approximation (no reduction of number of planes).
"""
@kwdef struct FullOrder{T} <: Algorithm
    metric::Symbol = defaultmetric()
    optimizer::T
end

"""
    FunctionEvaluations{D}

A structure holding a set of points and the associated
function values for a function f:ℜᴰ → ℜ.
"""
struct FunctionEvaluations{D,P<:NTuple{D,Number},V<:Number}
    points::Vector{P}
    values::Vector{V}
end
Base.length(feval::FunctionEvaluations) = length(feval.points)
function Base.iterate(feval::FunctionEvaluations, state = 0)
    state == length(feval) && return nothing
    return (feval.points[state+1], feval.values[state+1]), state + 1
end
Base.eltype(_::FunctionEvaluations{D,P,V}) where {D,P,V} = Tuple{P,V}

function Base.summary(io::IO, feval::FunctionEvaluations)
    n = length(feval)
    return print(io, "$(typeof(feval)) with $n data points")
end
function _write_points(io::IO, datapoints)
    for p ∈ datapoints
        println(io, " $(p[1]) ⟶ $(p[2])")
    end
end
function Base.show(
    io::IO,
    ::MIME"text/plain",
    feval::FunctionEvaluations{D},
) where {D}
    summary(io, feval)
    println(io, ":")
    fv = collect(feval)
    if length(feval) > 10
        _write_points(io, first(fv, 5))
        print(io, " ⋮\n")
        _write_points(io, last(fv, 5))
    else
        _write_points(io, fv)
    end
end

"""
    Plane{D}

A D-dimensional hyperplane (αᵀx = β).
"""
struct Plane{D}
    α::NTuple{D,Number}
    β::Number
end
Plane(a::NTuple{D}, b) where {D} = Plane{D}(a, b)
Plane(a::Vector, b) = Plane(Tuple(a), b)
Plane(a::Number, b::Number) = Plane{1}(Tuple(a), b)
evaluate(p::Plane, x, c = Convex) = dot(p.α, x) + p.β # Planes defined for convex functions
evaluate(p::Plane, x, c::Type{Concave}) = -evaluate(p, x, Convex) # Flip for concave functions

"""
    PWAFunc{C<:Curvature,D}

A piecewise linear function that can be either convex and concave
represented by a set of (hyper)planes.

If the curvature of the function is convex, the piecewise linear Function
is defined as the pointwise maximum over the planes. A concave pwa function
is handled by storing its negative and negating when evaluting.
"""
struct PWAFunc{C<:Curvature,D}
    planes::Vector{Plane{D}}
end
function PWAFunc(planes::Vector{Plane{D}}, C::Curvature) where {D}
    return PWAFunc{typeof(C),D}(planes)
end
PWAFunc{C,D}() where {C,D} = PWAFunc{C,D}(Vector{Plane{D}}())

function _write_planes(io::IO, C::Curvature, planes::Vector{Plane{D}}) where {D}
    for p ∈ planes
        println(io, " ", _write_plane(C, p))
    end
end
_write_plane(c::Convex, plane) = "z ≥ $(_xepr(plane.α))$(plane.β)"
_write_plane(c::Concave, plane) = "z ≤ $(_xepr(-1 .* plane.α))$(-plane.β)"
function _xepr(α)
    x = ["x₁", "x₂", "x₃", "x₄", "x₅"]
    expr = ""
    for i ∈ 1:min(5, length(α))
        expr *= "$(α[i]) $(x[i]) + "
    end
    if length(α) > 5
        expr *= "⋯ + "
    end
    return expr
end
function Base.summary(io::IO, pwa::PWAFunc)
    n = _planes(pwa)
    return print(io, "$(typeof(pwa)) with $n planes")
end
function Base.show(io::IO, ::MIME"text/plain", pwa::PWAFunc{C,D}) where {C,D}
    summary(io, pwa)
    println(io, ":")
    if length(pwa.planes) > 10
        _write_planes(io, C(), first(pwa.planes, 5))
        print(io, " ⋮\n")
        _write_planes(io, C(), last(pwa.planes, 5))
    else
        _write_planes(io, C(), pwa.planes)
    end
end

"""
    evaluate(pwa::PWAFunc{Convex,D}, x) where {D}

Evaluate the convex piecewise linear function at the point x.
"""
function evaluate(pwa::PWAFunc{Convex,D}, x) where {D}
    return maximum(dot(p.α, x) + p.β for p ∈ pwa.planes)
end

"""
    evaluate(pwa::PWAFunc{Concave,D}, x) where {D}

Evaluate the concave piecewise linear function at the point x.
"""
function evaluate(pwa::PWAFunc{Concave,D}, x) where {D}
    return -evaluate(PWAFunc{Convex,D}(pwa.planes), x)
end

# The number of planes defining the piecewise linar function
_planes(pwa) = length(pwa.planes)

_addplane!(pwa::PWAFunc{C,D}, p::Plane{D}) where {C,D} = push!(pwa.planes, p)
_addplane!(pwa::PWAFunc, α, β) = push!(pwa.planes, Plane(α, β))
