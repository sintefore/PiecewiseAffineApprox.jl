module JSON3Ext

using PiecewiseAffineApprox
using JSON3

# Overload read/write of JSON3 to use ST-compatible intermediate representation
JSON3.write(fn::String, pwaf::PWAFunc) = JSON3.write(fn, PiecewiseAffineApprox.STFunc(pwaf))
JSON3.read(s::AbstractString, ::Type{PWAFunc}) = PWAFunc(JSON3.read(s, PiecewiseAffineApprox.STFunc))


end