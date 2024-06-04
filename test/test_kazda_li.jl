
using PiecewiseAffineApprox

using JuMP
using HiGHS
using GLMakie

const optimizer =
    optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

# 1D
x = collect(range(-1, 1; length = 30))
z = x .^ 2

f = FunctionEvaluations(tuple.(x), z)
pwa_red = approx(
    f,
    Convex(),
    ProgressiveFitting(optimizer = optimizer, tolerance = 0.05, pen = :max),
)
plot(x, z, pwa_red)

# 2D
g(x) = x[1]^2 + x[2]^2
vals = PiecewiseAffineApprox._sample_uniform(g, [(-1, 1), (-1, 1)], 10)

pwa_red = approx(
    vals,
    Convex(),
    ProgressiveFitting(optimizer = optimizer, tolerance = 0.5, pen = :max),
)

plot(vals, pwl_red)
