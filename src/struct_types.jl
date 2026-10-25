# Non-parametric representation of planes for suppport from StructTypes
struct STPlane
    α::Vector{Float64}
    β::Float64
end
STPlane(p::Plane) = STPlane([i for i ∈ p.α], p.β)
StructTypes.StructType(::Type{STPlane}) = StructTypes.Struct()

# The domain as lower and upper vectors (empty without a domain), and back
function _st_domain(f)
    isnothing(f.domain) && return (Float64[], Float64[])
    return (first.(f.domain), last.(f.domain))
end
_domain_from_st(lower, upper) = isempty(lower) ? nothing : collect(zip(lower, upper))

# Non-parametric representation of PWAFunc for support from StructTypes. The
# struct is mutable with defaults, so that files without a domain can be read.
mutable struct STFunc
    curvature::Any
    planes::Vector{STPlane}
    lower::Vector{Float64}
    upper::Vector{Float64}
    STFunc() = new("", STPlane[], Float64[], Float64[])
    STFunc(curvature, planes, lower, upper) = new(curvature, planes, lower, upper)
end
function STFunc(pwaf::PWAFunc{Convex,N}) where {N}
    return STFunc(:convex, STPlane.(pwaf.planes), _st_domain(pwaf)...)
end
function STFunc(pwaf::PWAFunc{Concave,N}) where {N}
    return STFunc(:concave, STPlane.(pwaf.planes), _st_domain(pwaf)...)
end
StructTypes.StructType(::Type{STFunc}) = StructTypes.Mutable()

# Constructor to cast back to PWAFunc from STFunc
function PWAFunc(stf::STFunc)
    D = length(first(stf.planes).α)
    domain = _domain_from_st(stf.lower, stf.upper)
    if stf.curvature == "convex"
        return PWAFunc{Convex,D}(Plane.(stf.planes); domain)
    elseif stf.curvature == "concave"
        return PWAFunc{Concave,D}(Plane.(stf.planes); domain)
    else
        error(stf.curvature)
    end
end

Plane(fp::STPlane) = Plane(fp.α, fp.β)

# Non-parametric representation of DCFunc for support from StructTypes
mutable struct STDCFunc
    curvature::Any
    convex::Vector{STPlane}
    concave::Vector{STPlane}
    lower::Vector{Float64}
    upper::Vector{Float64}
    STDCFunc() = new("", STPlane[], STPlane[], Float64[], Float64[])
    function STDCFunc(curvature, convex, concave, lower, upper)
        return new(curvature, convex, concave, lower, upper)
    end
end
function STDCFunc(f::DCFunc)
    return STDCFunc(
        :dc,
        STPlane.(f.convex.planes),
        STPlane.(f.concave.planes),
        _st_domain(f)...,
    )
end
StructTypes.StructType(::Type{STDCFunc}) = StructTypes.Mutable()

# Constructor to cast back to DCFunc from STDCFunc
function DCFunc(st::STDCFunc)
    st.curvature == "dc" || error(st.curvature)
    planes = isempty(st.convex) ? st.concave : st.convex
    isempty(planes) && error("DCFunc without planes")
    D = length(first(planes).α)
    domain = _domain_from_st(st.lower, st.upper)
    return DCFunc(
        PWAFunc{Convex,D}(Plane.(st.convex); domain),
        PWAFunc{Concave,D}(Plane.(st.concave); domain);
        domain,
    )
end

# Construction from parsed JSON (dictionaries with string keys), for readers
# that do not fill in absent fields, such as JSON.jl
_st_plane(p) = STPlane(Float64.(p["α"]), Float64(p["β"]))
_st_vector(d, key) = Float64.(get(d, key, Float64[]))
function STFunc(d::AbstractDict)
    return STFunc(
        d["curvature"],
        [_st_plane(p) for p ∈ d["planes"]],
        _st_vector(d, "lower"),
        _st_vector(d, "upper"),
    )
end
function STDCFunc(d::AbstractDict)
    return STDCFunc(
        d["curvature"],
        [_st_plane(p) for p ∈ d["convex"]],
        [_st_plane(p) for p ∈ d["concave"]],
        _st_vector(d, "lower"),
        _st_vector(d, "upper"),
    )
end
