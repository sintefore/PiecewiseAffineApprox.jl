# Input types
abstract type Curvature end
struct Concave <: Curvature end
struct Convex <: Curvature end

abstract type Algorithm end
"""
    Heuristic
Compute affine approximation using the method proposed by Mangani & Boyd.

Note that this algoritm computes multiple approximations and selects the best.
If julia is started with multiple threads, these are computed in parallel. Consider
 how many threads will be beneficial, particularly when using a commercial solver where
 the license may restrict the number of simultanous solver instances.
"""
@kwdef struct Heuristic{T} <: Algorithm
    planes::Int = defaultplanes()
    pen::Symbol = defaultpenalty()
    trials::Int = 20
    itlim::Int = 50
    strict::Symbol = :none
    optimizer::T
end
@kwdef struct Interpol{T} <: Algorithm
    planes::Int = defaultplanes()
    pen = defaultpenalty()
    optimizer::T
end
"""
    Optimized
Compute affine approximation using a variation of the method proposed by Toriello & Vielma.

Note that the resulting approximation is sensitive to the selection of the Big-M value used when
solving the optimization problem.
"""
@kwdef struct Optimized{T} <: Algorithm
    planes::Int = defaultplanes()
    pen::Symbol = defaultpenalty()
    strict::Symbol = :none
    maxtime::Int = defaulttimelimit()
    optimizer::T
end

"""
    ProgressiveFitting

Compute affine approximation based on a variation of the method of Kazda and Li (2024)
specialized to convex approximations.
"""
@kwdef struct ProgressiveFitting{T} <: Algorithm
    tolerance::Float64 = 0.01
    pen::Symbol = defaultpenalty()
    optimizer::T
end

"""
    FunctionEvaluations{D}

A structure holding a set of points and the associated
function values for a function f:ℜᴰ → ℜ.
"""
struct FunctionEvaluations{D}
    points::Vector{<:NTuple{D,Number}}
    values::Vector{<:Number}
end

point_vals(f::FunctionEvaluations) = collect(zip(f.points, f.values))

"""
    Plane{D}

A D-dimensional hyperplane (αᵀx = β).
"""
struct Plane{D}
    α::NTuple{D,Number}
    β::Number
end
Plane(a::NTuple{N}, b) where {N} = Plane{N}(a, b)
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
