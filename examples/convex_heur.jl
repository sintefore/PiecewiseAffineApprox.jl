using GLMakie
using HiGHS

using PiecewiseAffineApprox

optimizer = HiGHS.Optimizer

# Function values on a uniform grid with n points in each dimension
function sample_grid(f, box, n)
    grid = Iterators.product((LinRange(lo, hi, n) for (lo, hi) ∈ box)...)
    points = vec(collect(grid))
    return FunctionEvaluations(points, f.(points))
end

# Quadratic with random sampling

# 1D
I = 100
x = 2 * rand(I) .- 1
z = x .^ 2
vals = FunctionEvaluations(Tuple.(x), z)
pwa = approx(vals, Convex(), Cluster(; optimizer = optimizer, planes = 5))

# 2D
xmat = 2 * rand(2, I) .- 1
x = [Tuple(col) for col ∈ eachcol(xmat)]
z = [p[1]^2 + p[2]^2 for p ∈ x]
vals = FunctionEvaluations(x, z)
pwa = approx(
    vals,
    Convex(),
    Cluster(; optimizer = optimizer, planes = 9, strict = :none),
)
plot(vals, pwa)

pwa = approx(vals, Convex(), FullOrder(; optimizer = optimizer, metric = :max))
plot(vals, pwa)

pwa = approx(
    vals,
    Convex(),
    Progressive(; optimizer = optimizer, tolerance = 0.5, metric = :l2),
)
plot(vals, pwa)

# Uniform sampling for various convex functions in 2D

# Quadratic
f(x) = x[1]^2 + x[2]^2
vals = sample_grid(f, [(-1, 1), (-1, 1)], 10)
pwa = approx(
    vals,
    Convex(),
    Cluster(; optimizer = optimizer, planes = 10, metric = :l2),
)
plot(vals, pwa)

# Concave quadratic
f(x) = 4 - x[1]^2 - x[2]^2
vals = sample_grid(f, [(-1, 1), (-1, 1)], 10)
pwa = approx(
    vals,
    Concave(),
    Cluster(; optimizer = optimizer, planes = 10, metric = :l1),
)

# Log sum exp
h(x) = log(exp(x[1]) + exp(x[2]))
vals = sample_grid(h, [(-1, 1), (-1, 1)], 10)
pwa = approx(vals, Convex(), Cluster(; optimizer = optimizer, planes = 5))
plot(vals, pwa)

# Quadratic over linear
g(x) = x[1]^2 / x[2]
vals = sample_grid(g, [(-1, 1), (0.1, 1)], 10)
pwa = approx(vals, Convex(), Cluster(; optimizer = optimizer, planes = 10))
plot(vals, pwa)

# Geometric mean
f(x) = -sqrt(x[1] * x[2])
vals = sample_grid(f, [(0.1, 1), (0.1, 1)], 10)
pwa = approx(vals, Convex(), Cluster(; optimizer = optimizer, planes = 10))
plot(vals, pwa)

# Non differentiable
f(x) = max(x[1]^2, x[2]^2)
vals = sample_grid(f, [(-1, 1), (-1, 1)], 10)
pwa = approx(
    vals,
    Convex(),
    Cluster(; optimizer = optimizer, planes = 8, metric = :l1),
)
plot(vals, pwa)

# Uniform sampling 3D
h(x) = log(exp(x[1]) + exp(x[2]) + exp(x[3]))
vals = sample_grid(h, [(-5, 5), (-5, 5), (-5, 5)], 10)
pwa = approx(
    vals,
    Convex(),
    Cluster(; optimizer = optimizer, planes = 8, metric = :rms),
)
