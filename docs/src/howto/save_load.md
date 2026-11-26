# Save and load

Fitting can take time, so approximations are often computed once and stored.
When [JSON](https://github.com/JuliaIO/JSON.jl) is loaded, a [`PWAFunc`](@ref)
or [`DCFunc`](@ref) is written with `JSON.json(filename, f)` and read with
`JSON.parsefile(filename, PWAFunc)` or `JSON.parsefile(filename, DCFunc)`:

```jldoctest save
using PiecewiseAffineApprox, JSON

pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 4))
file = joinpath(mktempdir(), "pwa.json")
JSON.json(file, pwa)
JSON.parsefile(file, PWAFunc) == pwa

# output

true
```

The file holds the curvature, the planes with their coefficients `α` and
constant `β`, and the bounds of the domain. For `|x|` on `[-1, 1]`, built from
its two planes:

```jldoctest save
absval = PWAFunc{Convex,1}([Plane((-1.0,), 0.0), Plane((1.0,), 0.0)]; domain = [(-1.0, 1.0)])
JSON.json(file, absval)
print(read(file, String))

# output

{"curvature":"convex","planes":[{"α":[-1.0],"β":0.0},{"α":[1.0],"β":0.0}],"lower":[-1.0],"upper":[1.0]}
```

A [`DCFunc`](@ref) is written with the planes of its `convex` and `concave`
parts. Files written with the JSON3 extension of versions before 0.7 have the
same format and are read in the same way.
