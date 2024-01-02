# Input types
abstract type Curvature end
struct Concave <: Curvature end
struct Convex <: Curvature end

abstract type Algorithm end
struct Heuristic <: Algorithm end
struct Interpol <: Algorithm end
struct Optimized <: Algorithm end

"""
    FunctionEvaluations{D}

A structure holding a set of points and the associated
function values for a function f:ℜᴰ → ℜ.
"""
struct FunctionEvaluations{D}
    points::Vector{<:NTuple{D,Number}}
    values::Vector{<:Number}
end

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
evaluate(p::Plane, x, c=Convex) = dot(p.α, x) + p.β # Planes defined for convex functions
evaluate(p::Plane, x, c::Type{Concave}) = -evaluate(p, x, Convex) # Flip for concave functions 

"""
    PWLFunc{C<:Curvature,D}

A piecewise linear function that can be either convex and concave
represented by a set of (hyper)planes.

If the curvature of the function is convex, the piecewise linear Function
is defined as the pointwise maximum over the planes. A concave pwl function
is handled by storing its negative and negating when evaluting.
"""
struct PWLFunc{C<:Curvature,D}
    planes::Vector{Plane{D}}
end
function PWLFunc(planes::Vector{Plane{D}}, C::Curvature) where {D}
    return PWLFunc{typeof(C),D}(planes)
end
PWLFunc{C,D}() where {C,D} = PWLFunc{C,D}(Vector{Plane{D}}())

"""
    evaluate(pwl::PWLFunc{Convex,D}, x) where {D}

Evaluate the convex piecewise linear function at the point x.
"""
function evaluate(pwl::PWLFunc{Convex,D}, x) where {D}
    return maximum(dot(p.α, x) + p.β for p ∈ pwl.planes)
end

"""
    evaluate(pwl::PWLFunc{Concave,D}, x) where {D}

Evaluate the concave piecewise linear function at the point x.
"""
function evaluate(pwl::PWLFunc{Concave,D}, x) where {D}
    return -evaluate(PWLFunc{Convex,D}(pwl.planes), x)
end

# The number of planes defining the piecewise linar function
_planes(pwl) = length(pwl.planes)

_addplane!(pwl::PWLFunc{C,D}, p::Plane{D}) where {C,D} = push!(pwl.planes, p)
_addplane!(pwl::PWLFunc, α, β) = push!(pwl.planes, Plane(α, β))
