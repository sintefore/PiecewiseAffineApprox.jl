
using PiecewiseAffineApprox

using JuMP
using HiGHS
using GLMakie

const optimizer =
    optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

# 1D
x = collect(range(-1, 1; length = 30))
z = x.^2

f = FunctionEvaluations(tuple.(x), z)
pwl_red = PiecewiseAffineApprox._progressive_pwa(f, optimizer, 0.05, :max)
plot(x, z, pwl_red)


# 2D
g(x) = x[1]^2 + x[2]^2
vals = PiecewiseAffineApprox._sample_uniform(g, [(-1, 1), (-1, 1)], 10)

pwl_red = PiecewiseAffineApprox._progressive_pwa(vals, optimizer, 0.5, :max)


pwl_red = approx(vals, Convex(), ProgressiveFitting(optimizer = optimizer, tolerance = 2.5, pen = :l1))


plot(vals, pwl_red)
