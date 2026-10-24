# Non-parametric representation of planes for suppport from StructTypes
struct STPlane
    α::Vector{Float64}
    β::Float64
end
STPlane(p::Plane) = STPlane([i for i ∈ p.α], p.β)
StructTypes.StructType(::Type{STPlane}) = StructTypes.Struct()

# Non-parametric representation of PWAFunc for support from StructTypes
struct STFunc
    curvature::Any
    planes::Vector{STPlane}
end
function STFunc(pwaf::PWAFunc{Convex,N}) where {N}
    return STFunc(:convex, STPlane.(pwaf.planes))
end
function STFunc(pwaf::PWAFunc{Concave,N}) where {N}
    return STFunc(:concave, STPlane.(pwaf.planes))
end
StructTypes.StructType(::Type{STFunc}) = StructTypes.Struct()

# Constructor to cast back to PWAFunc from STFunc
function PWAFunc(stf::STFunc)
    D = length(first(stf.planes).α)
    if stf.curvature == "convex"
        return PWAFunc{Convex,D}(Plane.(stf.planes))
    elseif stf.curvature == "concave"
        return PWAFunc{Concave,D}(Plane.(stf.planes))
    else
        error(stf.curvature)
    end
end

Plane(fp::STPlane) = Plane(fp.α, fp.β)

# Non-parametric representation of DCFunc for support from StructTypes
struct STDCFunc
    curvature::Any
    convex::Vector{STPlane}
    concave::Vector{STPlane}
end
function STDCFunc(f::DCFunc)
    return STDCFunc(:dc, STPlane.(f.convex.planes), STPlane.(f.concave.planes))
end
StructTypes.StructType(::Type{STDCFunc}) = StructTypes.Struct()

# Constructor to cast back to DCFunc from STDCFunc
function DCFunc(st::STDCFunc)
    st.curvature == "dc" || error(st.curvature)
    planes = isempty(st.convex) ? st.concave : st.convex
    isempty(planes) && error("DCFunc without planes")
    D = length(first(planes).α)
    return DCFunc(
        PWAFunc{Convex,D}(Plane.(st.convex)),
        PWAFunc{Concave,D}(Plane.(st.concave)),
    )
end
