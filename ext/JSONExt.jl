module JSONExt

using PiecewiseAffineApprox
using JSON

# Overload read/write of JSON to use ST-compatible intermediate representation that also stores curvature
function JSON.json(fn::String, pwaf::PWAFunc)
    return JSON.json(fn, PiecewiseAffineApprox.STFunc(pwaf))
end
function JSON.parsefile(s::AbstractString, ::Type{PWAFunc})
    return PWAFunc(PiecewiseAffineApprox.STFunc(JSON.parsefile(s)))
end
function JSON.json(fn::String, f::DCFunc)
    return JSON.json(fn, PiecewiseAffineApprox.STDCFunc(f))
end
function JSON.parsefile(s::AbstractString, ::Type{DCFunc})
    return DCFunc(PiecewiseAffineApprox.STDCFunc(JSON.parsefile(s)))
end

end
