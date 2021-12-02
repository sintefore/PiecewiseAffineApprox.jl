using GLMakie
using Xpress

using PiecewiseLinearApprox

const PWL = PiecewiseLinearApprox

I = 100

x = rand(2,I)
z = [x[1,i]^2 + x[2,i]^2 for i in 1:I]

pwl = convex_linearization_mb(x, z, GLPK.Optimizer; nseg=4)

x = 2*rand(1,I) .- 1
z = x.^2

pwl = convex_linearization_mb(x, z, GLPK.Optimizer; nseg=5)

# Quadratic with random sampling
I = 100

# 1D
x = 2 * rand(I) .- 1
sort!(x)
z = x.^2
pwl = PWL.convex_linearization_mb(x', z, Xpress.Optimizer; nseg=5, pen=:l1, trials=50)
pwl_opt = PWL.convex_linearization(x, z, Xpress.Optimizer; nseg=5, pen=:l1)

# 2D
x = 2 * rand(2,I) .- 1
z = [x[1,i]^2 + x[2,i]^2 for i in 1:I]
pwl = PWL.convex_linearization_mb(x, z, Xpress.Optimizer; nseg=7, strict=:below)
PWL.plotconv2D(x,z,pwl,-1,1,-1,1)

# Uniform sampling for various convex functions in 2D

# Quadratic
f(x) = x[1]^2 + x[2]^2
vals = PWL.sample_uniform(f, [(-1,1), (-1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), solver=Xpress.Optimizer, nseg=10)
PWL.plotconv2D(vals, pwl)

# Log sum exp
h(x) = log(exp(x[1]) + exp(x[2]))
vals = PWL.sample_uniform(h, [(-1,1), (-1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), solver=Xpress.Optimizer, nseg=5)
PWL.plotconv2D(vals, pwl)

# Quadratic over linear
g(x) = x[1]^2 / x[2]
vals = PWL.sample_uniform(g, [(-1,1), (0.1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), solver=Xpress.Optimizer, nseg=10)
PWL.plotconv2D(vals, pwl)

# Geometric mean
f(x) = -sqrt(x[1]*x[2])
vals = PWL.sample_uniform(f, [(0.1,1), (0.1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), solver=Xpress.Optimizer, nseg=10)
PWL.plotconv2D(vals, pwl)

# Non differentiable
f(x) = max(x[1]^2,x[2]^2)
vals = PWL.sample_uniform(f, [(-1,1), (-1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), solver=Xpress.Optimizer, nseg=8, pen=:max)
PWL.plotconv2D(vals, pwl)


# Uniform sampling 3D
h(x) = log(exp(x[1]) + exp(x[2]) + exp(x[3]))
vals = PWL.sample_uniform(h, [(-5,5), (-5,5), (-5,5)], 10)
pwl = approx(vals, Convex(), Heuristic(), solver=Xpress.Optimizer, nseg=8, pen=:rms)