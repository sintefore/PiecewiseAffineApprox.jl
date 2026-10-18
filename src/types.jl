# Input types
abstract type Curvature end
struct Concave <: Curvature end
struct Convex <: Curvature end

abstract type Algorithm end
"""
    Cluster
Compute affine approximation using the method proposed by Magnani & Boyd.

Note that this algoritm computes multiple approximations and selects the best.
If julia is started with multiple threads, these are computed in parallel. Consider
 how many threads will be beneficial, particularly when using a commercial solver where
 the license may restrict the number of simultanous solver instances.

The random initial partitions are generated from `rng`. Pass a seeded generator,
e.g. `rng = StableRNG(123)` from StableRNGs.jl, to get reproducible results that are
independent of the number of threads. Note that `rng` is advanced by each call to
`approx`, so reusing the same `Cluster` object gives a different (but still
reproducible) sequence of results.

With `metric = :l2` (the default, or `:rms`) and `strict = :none`, each plane is
fitted with closed form least squares as in the original method. This is much
faster and does not require an `optimizer`. Other metrics, and
`strict = :inner`/`:outer`, solve an optimization problem for each fit and require
an `optimizer`.

!!! warning "Use a robust QP solver for `:l2` with `strict`"
    With `metric = :l2` and `strict = :inner` or `:outer`, each fit is a convex
    quadratic program. The QP solver in HiGHS (tested with HiGHS 1.15.1) fails on
    some of these problems: it either reports the problem as non-convex or does not
    terminate. Use another QP solver for this combination, e.g. Ipopt or Xpress
    (both tested), and consider setting a time limit. HiGHS works well for the
    linear programs solved with `metric = :l1` and `:max`.

Each trial returns the best approximation found during its iterations, and the best
over all trials is selected, both as measured by `metric`. For the fits that solve
an optimization problem, each trial first refines its partition of the points with
closed form least squares, and then for at most five iterations with the requested
fit, which is both faster and better than the requested fit alone.
"""
@kwdef struct Cluster{T,R<:AbstractRNG} <: Algorithm
    planes::Int = defaultplanes()
    metric::Symbol = :l2
    trials::Int = 20
    itlim::Int = 50
    strict::Symbol = :none
    optimizer::T = nothing
    rng::R = Random.default_rng()
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

The knots of the interpolant are data points, or adjusted after fitting by
`adjust`:

- `:none` (default): keep the knots.
- `:values`: keep the breakpoints, and refit the function values at the knots to
  minimize the deviation from the data, subject to convexity.
- `:breakpoints`: also move the interior breakpoints by a local search, refitting
  the function values for each candidate position.

Adjusting never increases the deviation from the data. For `:l2` it is solved as
a nonnegative least squares problem, while `:l1` and `:max` require an
`optimizer`. When the number of planes is selected, the selection is based on
the unadjusted fits.
"""
@kwdef struct Interpol{T} <: Algorithm
    planes::Int = defaultplanes()
    metric::Symbol = defaultmetric()
    optimizer::T = nothing
    convexify::Symbol = :auto
    selection::Symbol = :none
    tolerance::Float64 = 0.0
    folds::Int = 5
    adjust::Symbol = :none
end

"""
    MILP
Compute affine approximation using a variation of the method proposed by Toriello & Vielma.

The deviation from the data is measured by `metric` (`:l1`, `:l2` or `:max`), and
`strict` restricts the approximation to lie below the data (`:outer`), above the
data (`:inner`) or neither (`:none`, default). For concave approximations, below
and above are interchanged.

With `warmstart = true` (default), a heuristic solution from [`Cluster`](@ref) is
passed to the solver as a start. It is also used to bound the deviation from
the data in each point, and to make sure the big-M values do not cut it off.
The warm start is required for more than two dimensions.

The problem is solved on data scaled to the unit box, so the result does not
depend on the units of the data. The big-M values in the formulation are
heuristic. If they may cut off the optimal solution, i.e. a big-M constraint is
tight for a plane that is not active, they are doubled and the problem is solved
again as long as this improves the objective, within the total time `maxtime`.
"""
@kwdef struct MILP{T} <: Algorithm
    planes::Int = defaultplanes()
    metric::Symbol = defaultmetric()
    strict::Symbol = :none
    maxtime::Int = defaulttimelimit()
    optimizer::T
    warmstart::Bool = true
end

"""
    Progressive

Compute affine approximation based on a variation of the method of Kazda and Li (2024)
specialized to convex approximations.

Planes are added until the approximation meets the `tolerance`, in one of two
ways depending on `strict`:

- `:outer` (default): the planes lie below all data points, and planes are added
  until the error per data point is at most `tolerance`, measured by `metric`:
  the mean absolute deviation (`:l1`), the root mean square deviation (`:l2`) or
  the maximum absolute deviation (`:max`). Each plane is fitted to its points
  in the same metric (for `:l2` with two or more inputs, the absolute deviation
  is used for the fit). The data points must be convex, i.e.
  possible to interpolate by a convex function. Use [`enforce_curvature`](@ref)
  first for noisy or non-convex data.
- `:none`: the approximation lies within `tolerance` of all data points (above
  or below), as in the original method of Kazda and Li. The data does not need
  to be convex, as long as some convex function lies within the tolerance, so
  noise up to the tolerance is handled. `metric` is not used.
- `:inner`: as `:none`, but the approximation lies between the data and
  `tolerance` above it.

For concave approximations, above and below are interchanged.
"""
@kwdef struct Progressive{T} <: Algorithm
    tolerance::Float64 = 0.01
    metric::Symbol = defaultmetric()
    optimizer::T
    strict::Symbol = :outer
end

"""
    FullOrder

Compute affine approximation based on a variation of the method of Kazda and Li (2024)
specialized to convex approximations with full order approximation (no reduction of number of planes).

Each data point gets a plane through the point that lies below all other data
points. Among such planes, the one minimizing the deviation from the other points
is chosen, measured by `metric` (`:l1`, `:l2`, `:max` or `:none`).

The data points must be convex, i.e. possible to interpolate by a convex
function. Use [`enforce_curvature`](@ref) first for noisy or non-convex data.
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
    Plane{D,T<:Real}

A D-dimensional affine function z = αᵀx + β, with coefficients of type `T`.
"""
struct Plane{D,T<:Real}
    α::NTuple{D,T}
    β::T
end
function Plane(a::Tuple{Vararg{Real,D}}, b::Real) where {D}
    T = promote_type(typeof(b), map(typeof, a)...)
    return Plane{D,T}(a, b)
end
Plane{D}(a::Tuple{Vararg{Real,D}}, b::Real) where {D} = Plane(a, b)
Plane(a::AbstractVector{<:Real}, b::Real) = Plane(Tuple(a), b)
Plane(a::Real, b::Real) = Plane((a,), b)

Base.convert(::Type{Plane{D,T}}, p::Plane{D}) where {D,T} = Plane{D,T}(p.α, p.β)
function Base.promote_rule(::Type{Plane{D,S}}, ::Type{Plane{D,T}}) where {D,S,T}
    return Plane{D,promote_type(S, T)}
end

evaluate(p::Plane, x, c = Convex) = dot(p.α, x) + p.β # Planes defined for convex functions
evaluate(p::Plane, x, c::Type{Concave}) = -evaluate(p, x, Convex) # Flip for concave functions

"""
    PWAFunc{C<:Curvature,D,T<:Real}

A piecewise linear function that can be either convex and concave
represented by a set of (hyper)planes with coefficients of type `T`.

If the curvature of the function is convex, the piecewise linear Function
is defined as the pointwise maximum over the planes. A concave pwa function
is handled by storing its negative and negating when evaluting.

`PWAFunc{C,D}(planes)` infers `T` from the planes, and `PWAFunc{C,D}()` creates
an empty function with `T = Float64`.
"""
struct PWAFunc{C<:Curvature,D,T<:Real}
    planes::Vector{Plane{D,T}}
end
function PWAFunc{C,D}(planes::AbstractVector{<:Plane{D}}) where {C<:Curvature,D}
    T = isempty(planes) ? _numtype(eltype(planes)) :
        mapreduce(p -> typeof(p.β), promote_type, planes)
    return PWAFunc{C,D,T}(convert(Vector{Plane{D,T}}, planes))
end
PWAFunc{C,D}() where {C<:Curvature,D} = PWAFunc{C,D,Float64}(Plane{D,Float64}[])
function PWAFunc(planes::AbstractVector{<:Plane{D}}, C::Curvature) where {D}
    return PWAFunc{typeof(C),D}(planes)
end

_numtype(::Type{Plane{D,T}}) where {D,T} = T
_numtype(::Type) = Float64

function _write_planes(io::IO, C::Curvature, planes::AbstractVector{<:Plane})
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
function Base.show(io::IO, ::MIME"text/plain", pwa::PWAFunc{C}) where {C}
    summary(io, pwa)
    println(io, ":")
    if length(pwa.planes) > 10
        _write_planes(io, C(), first(pwa.planes, 5))
        print(io, " ⋮
")
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
    return -maximum(dot(p.α, x) + p.β for p ∈ pwa.planes)
end

# The number of planes defining the piecewise linar function
_planes(pwa) = length(pwa.planes)

_addplane!(pwa::PWAFunc{C,D}, p::Plane{D}) where {C,D} = push!(pwa.planes, p)
_addplane!(pwa::PWAFunc, α, β) = push!(pwa.planes, Plane(α, β))
