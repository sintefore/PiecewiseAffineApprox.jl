module JSONExt

using PiecewiseAffineApprox
using JSON

# Reading and writing of PWAFunc and DCFunc as JSON. A function is written as
# an object with its curvature ("convex", "concave" or "dc"), its planes (for
# a DCFunc those of the convex and of the concave part) with coefficients α and
# constant β, and the lower and upper bounds of its domain, which are empty for
# a function without a domain and may be absent in a file.

_plane(p::Plane) = (α = collect(p.α), β = p.β)
_plane(d::AbstractDict) = Plane(Float64.(d["α"]), Float64(d["β"]))

# The domain as lower and upper vectors, and back
function _bounds(f)
    isnothing(f.domain) && return (lower = Float64[], upper = Float64[])
    return (lower = collect(first.(f.domain)), upper = collect(last.(f.domain)))
end
function _domain(d::AbstractDict)
    lower = Float64.(get(d, "lower", Float64[]))
    upper = Float64.(get(d, "upper", Float64[]))
    return isempty(lower) ? nothing : collect(zip(lower, upper))
end

_curvature(::PWAFunc{Convex}) = "convex"
_curvature(::PWAFunc{Concave}) = "concave"

function _serializable(f::PWAFunc)
    return (; curvature = _curvature(f), planes = _plane.(f.planes), _bounds(f)...)
end
function _serializable(f::DCFunc)
    return (;
        curvature = "dc",
        convex = _plane.(f.convex.planes),
        concave = _plane.(f.concave.planes),
        _bounds(f)...,
    )
end

function _pwafunc(d::AbstractDict)
    planes = _plane.(d["planes"])
    isempty(planes) && error("PWAFunc without planes")
    D = length(first(planes).α)
    c = d["curvature"]
    C = c == "convex" ? Convex : c == "concave" ? Concave : error(c)
    return PWAFunc{C,D}(planes; domain = _domain(d))
end

function _dcfunc(d::AbstractDict)
    d["curvature"] == "dc" || error(d["curvature"])
    convex, concave = _plane.(d["convex"]), _plane.(d["concave"])
    planes = isempty(convex) ? concave : convex
    isempty(planes) && error("DCFunc without planes")
    D = length(first(planes).α)
    domain = _domain(d)
    # A part without planes has an empty list of the same type of planes
    P = eltype(planes)
    part(C, list) = PWAFunc{C,D}(P[p for p ∈ list]; domain)
    return DCFunc(part(Convex, convex), part(Concave, concave); domain)
end

JSON.json(fn::String, f::PWAFunc) = JSON.json(fn, _serializable(f))
JSON.json(fn::String, f::DCFunc) = JSON.json(fn, _serializable(f))
function JSON.parsefile(s::AbstractString, ::Type{PWAFunc})
    return _pwafunc(JSON.parsefile(s))
end
JSON.parsefile(s::AbstractString, ::Type{DCFunc}) = _dcfunc(JSON.parsefile(s))

end
