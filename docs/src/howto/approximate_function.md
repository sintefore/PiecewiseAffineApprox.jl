# Approximate a function

A function is approximated by sampling it on a uniform grid over a box and
fitting the samples. [`approx`](@ref) takes the function, the box with one
`(lower, upper)` tuple for each variable, the shape of the approximation and
the algorithm:

```jldoctest function
using PiecewiseAffineApprox

pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 5))
length(pwa.planes)

# output

5
```

The result is a [`PWAFunc`](@ref): the maximum of its planes for a convex
approximation, and the minimum for a concave one. It is evaluated by calling
it with numbers, and added to a JuMP model by calling it with variables, see
[Add to a JuMP model](@ref).

```jldoctest function
round(pwa(0.5); digits = 4), round(evaluate(pwa, (0.5,)); digits = 4)

# output

(0.2755, 0.2755)
```

## Functions of several variables

The function is called with the point as a tuple if it accepts one, e.g.
`x -> x[1]^2 + x[2]^2`, and otherwise with the coordinates as separate
arguments, e.g. `sin` or `(x, y) -> x^2 + y^2`. [`Interpol`](@ref) handles one
variable only, while [`Cluster`](@ref), [`MILP`](@ref),
[`Progressive`](@ref) and [`FullOrder`](@ref) handle any number.

!!! tip "Functions of one variable"
    A function with an untyped argument, e.g. `f(x) = exp(x) - 2x`, accepts a
    tuple, so it is called with the tuple `(x,)` and fails. Either index the
    point, `x -> exp(x[1]) - 2x[1]`, or restrict the argument to numbers,
    `f(x::Real) = exp(x) - 2x`.


```jldoctest function
using StableRNGs

f(x, y) = x^2 + y^2 + x * y
pwa = approx(f, [(-1, 1), (-1, 1)], Convex(), Cluster(planes = 8, rng = StableRNG(1)))
length(pwa.planes), domain(pwa)

# output

(8, [(-1.0, 1.0), (-1.0, 1.0)])
```

The approximation records the box as its [`domain`](@ref), which is used as
the default bounds of the inputs when it is added to a model, and for plotting.

## The number of samples

The keyword `nsamples` sets the number of points in each dimension. By default
it depends on the algorithm and the number of variables `D`:

- [`MILP`](@ref): about `3 * D * planes` points in total, since the size of
  the optimization problem grows with the number of points,
- [`Cluster`](@ref) and [`Interpol`](@ref): `3 * planes` in each dimension,
- other algorithms: about 100 points in total.

More samples give a better picture of the function, at the cost of a longer
fit:

```jldoctest function
x = range(0, 2; length = 1000)
test = FunctionEvaluations(x, exp.(x))
errors = map((5, 50)) do n
    pwa = approx(x -> exp(x[1]), [(0, 2)], Convex(), Interpol(planes = 4); nsamples = n)
    round(approximation_error(pwa, test, :max); digits = 3)
end

# output

(0.181, 0.104)
```

Here [`approximation_error`](@ref) measures the largest deviation from the
function on a fine grid.

## Concave and non-convex functions

A concave function is approximated in the same way with [`Concave`](@ref),
and a function that is neither convex nor concave with [`DC`](@ref), as the
sum of a convex and a concave part, see [Non-convex functions](@ref):

```jldoctest function
g = approx(x -> log(x[1]), [(0.5, 4)], Concave(), Interpol(planes = 4))
h = approx(sin, [(0, 2π)], DC(), Interpol(planes = 4, concave_planes = 3))
summary(g), summary(h)

# output

("PWAFunc{Concave, 1, Float64} with 4 planes on [0.5, 4.0]", "DCFunc{1, Float64} with 4 convex and 3 concave planes on [0.0, 6.283]")
```

## Diagnostics

[`fit`](@ref) takes the same arguments as `approx`, and returns a
[`FitResult`](@ref) with the approximation `func`, its `error` from the
samples, the run time and `details` from the algorithm:

```jldoctest function
r = fit(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 5))
round(r.error; digits = 4), r.metric

# output

(0.0231, :l1)
```

To make the approximation lie on one side of the function, at the samples or
over the whole box, see [Control the error](@ref).
