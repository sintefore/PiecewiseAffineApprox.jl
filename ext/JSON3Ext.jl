module JSON3Ext

using PiecewiseAffineApprox
using JSON3

# Overload read/write of JSON3 to use ST-compatible intermediate representation
function JSON3.write(fn::String, pwaf::PWAFunc)
    return JSON3.write(fn, PiecewiseAffineApprox.STFunc(pwaf))
end
function JSON3.read(s::AbstractString, ::Type{PWAFunc})
    return PWAFunc(JSON3.read(s, PiecewiseAffineApprox.STFunc))
end

end
