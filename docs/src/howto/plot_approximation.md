# Plot approximations

To keep the dependencies light, plotting is provided by package extensions,
available when [Makie](https://docs.makie.org/stable/) (e.g. through
CairoMakie or GLMakie) or [Plots](https://docs.juliaplots.org/stable/) is
loaded.

## Makie

A function of one variable is plotted together with the data, with its planes
dashed:

```@example plot
using PiecewiseAffineApprox, CairoMakie

x = range(0, 2; length = 20)
f(x::Real) = exp(x) - 2x
pwa = approx(f, [(0, 2)], Convex(), Interpol(planes = 4))
plot(x, f.(x), pwa)
```

An approximation fitted with [`approx`](@ref) records its domain, so it can
also be plotted on its own, `plot(pwa)`.

A function of two variables is shown as a surface, together with the error at
each data point and the points allocated to each plane:

```@example plot
using StableRNGs

rng = StableRNG(1)
points = [Tuple(2 .* rand(rng, 2) .- 1) for _ ∈ 1:200]
data = FunctionEvaluations(points, [x^2 + y^2 for (x, y) ∈ points])
pwa2 = approx(data, Convex(), Cluster(planes = 9, rng = StableRNG(2)))
plot(data, pwa2)
```

A [`DCFunc`](@ref) is shown with its convex and concave parts in separate
axes, together with their planes; pass `parts = false` to show only the
approximation and the data:

```@example plot
x = range(-1, 1; length = 40)
g(x::Real) = sin(3x)
dc = approx(g, [(-1, 1)], DC(), Interpol(planes = 4, concave_planes = 3))
plot(x, g.(x), dc)
```

## Plots

With Plots, a function of one variable is plotted over an interval,
`plot(pwa, (lo, hi))`, and added to an existing plot with
`plot!(p, pwa, (lo, hi))`:

```julia
using Plots
p = scatter(x, f.(x); label = "data")
plot!(p, pwa, (0, 2))
```

A [`DCFunc`](@ref) is plotted in the same way, with its parts unless
`parts = false`.

## Animations

The accuracy when adding more planes:

```@raw html
<video autoplay loop muted playsinline style="max-width: 100%">
<source src="../../assets/approxanim.mp4" type="video/mp4" />
</video>
```

An approximation of a function of two variables:

```@raw html
<video autoplay loop muted playsinline style="max-width: 100%">
<source src="../../assets/rotation.mp4" type="video/mp4" />
</video>
```
