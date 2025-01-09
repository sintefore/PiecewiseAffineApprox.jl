# Non-parametric representation of planes for suppport from StructTypes
struct STPlane
    α::Vector{Float64}
    β::Float64
end
STPlane(p::Plane) = STPlane([i for i in p.α],p.β)
StructTypes.StructType(::Type{STPlane}) = StructTypes.Struct()

# Non-parametric representation of PWAFunc for support from StructTypes
struct STFunc
    curvature
    planes::Vector{STPlane}
end
STFunc(pwaf::PWAFunc{Convex,N}) where {N} = STFunc(:convex, STPlane.(pwaf.planes))
STFunc(pwaf::PWAFunc{Concave,N}) where{N} = STFunc(:concave, STPlane,(pwaf.planes))
StructTypes.StructType(::Type{STFunc}) = StructTypes.Struct()

# Constructor to cast back to PWAFunc from STFunc
function PWAFunc(stf::STFunc)
	D = length(first(stf.planes).α)
	if stf.curvature == "convex"
		return PWAFunc{Convex, D}(Plane.(stf.planes))
	elseif stf.curvature == "concave"
		return PWAFunc{Concave, D}(Plane.(stf.planes))
	else
		error(stf.curvature)
	end
end

Plane(fp::STPlane) = Plane(fp.α,fp.β)

# Overload read/write of JSON3 to use ST-compatible intermediate representation
# (TODO: move to extension)
JSON3.write(fn::String, pwaf::PWAFunc) = JSON3.write(fn, STFunc(pwaf))
JSON3.read(s::AbstractString, ::Type{PWAFunc}) = PWAFunc(JSON3.read(s, STFunc))