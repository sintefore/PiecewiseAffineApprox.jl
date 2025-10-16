module JSONExt

using PiecewiseAffineApprox
using JSON

# Overload read/write of JSON to use ST-compatible intermediate representation that also stores curvature
function JSON.json(fn::String, pwaf::PWAFunc)
    return JSON.json(fn, PiecewiseAffineApprox.STFunc(pwaf))
end
function JSON.parsefile(s::AbstractString, ::Type{PWAFunc})
    return PWAFunc(JSON.parsefile(s, PiecewiseAffineApprox.STFunc))
end

end
