# Approximate data

Measured or precomputed data is wrapped in [`FunctionEvaluations`](@ref), the
points together with their function values, and passed to [`approx`](@ref).
The points can be given

- as a vector of numbers, for a function of one variable,
- as a vector of tuples, one tuple per point,
- as a matrix with one point per row (`dims = 1`, the default) or per column
  (`dims = 2`).

```jldoctest points
using PiecewiseAffineApprox

x = collect(range(-1, 1; length = 30))
z = x .^ 2 .+ 0.02 .* sin.(40 .* x)       # convex with some noise
data = FunctionEvaluations(x, z)

pwa = approx(data, Convex(), Interpol(planes = 4))
length(pwa.planes), domain(pwa)

# output

(4, [(-1.0, 1.0)])
```

The range of the data is recorded as the [`domain`](@ref) of the
approximation. The points do not need to be sorted, and repeated measurements
at the same point are allowed.

## Several variables

For data in several variables, use [`Cluster`](@ref) for many points,
[`MILP`](@ref) for the best fit of a small data set, or [`Progressive`](@ref)
to meet a tolerance. Here 200 random points of a convex function of two
variables, given as a matrix with one point per row:

```jldoctest points
using StableRNGs

rng = StableRNG(1)
X = 2 .* rand(rng, 200, 2) .- 1
z = [x^2 + y^2 + x * y for (x, y) ∈ eachrow(X)]
data2 = FunctionEvaluations(X, z)

pwa = approx(data2, Convex(), Cluster(planes = 6, rng = StableRNG(2)))
length(pwa.planes)

# output

6
```

[`Cluster`](@ref) starts from random partitions of the points. Pass a seeded
generator such as `StableRNG(2)` from StableRNGs.jl as `rng` for reproducible
results.

## Concave data

Concave data is approximated with [`Concave`](@ref), as the minimum of planes:

```jldoctest points
pwa = approx(FunctionEvaluations(x, -x .^ 2), Concave(), Interpol(planes = 4))
pwa(0.0) ≤ 0

# output

true
```

## Noisy and non-convex data

The data does not have to be exactly convex:

- [`MILP`](@ref), [`Cluster`](@ref) and [`Interpol`](@ref) fit any data with
  a convex function. Interpol can make the data convex first with
  `convexify = :data`, see [Interpol](@ref).
- [`Progressive`](@ref) with `strict = :none` or `:inner` handles noise up to
  its tolerance.
- [`Progressive`](@ref) with the default `strict = :outer`, and
  [`FullOrder`](@ref), need data that can be interpolated by a convex function.
  [`enforce_curvature`](@ref) changes the values as little as possible to make
  it so:

```jldoctest points
using HiGHS, JuMP
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

convex = enforce_curvature(data, Convex(), optimizer)
pwa = approx(convex, Convex(), Progressive(; optimizer, tolerance = 0.02))
length(pwa.planes)

# output

4
```

Data that is clearly neither convex nor concave is approximated with
[`DC`](@ref), see [Non-convex functions](@ref).

## Check the fit

[`approximation_error`](@ref) measures the deviation of an approximation from
data, and [`fit`](@ref) reports it together with diagnostics from the
algorithm, see [Control the error](@ref):

```jldoctest points
r = fit(data, Convex(), Interpol(planes = 4, metric = :l2))
round(r.error; digits = 4), round(approximation_error(r.func, data, :max); digits = 4)

# output

(0.0384, 0.0884)
```
