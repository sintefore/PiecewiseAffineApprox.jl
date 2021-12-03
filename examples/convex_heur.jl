using GLMakie
using Xpress

using PiecewiseLinearApprox

const PWL = PiecewiseLinearApprox

# Quadratic with random sampling

# 1D
I = 100
x = 2 * rand(I) .- 1
z = x.^2
vals = FunctionEvaluations(Tuple.(x), z)
pwl = approx(vals, Convex(), Heuristic(); optimizer=Xpress.Optimizer, planes=5)

# 2D
xmat = 2 * rand(2,I) .- 1
x = [Tuple(xmat[:,i]) for i in 1:size(xmat,2)]
z = [p[1]^2 + p[2]^2 for p in x]
vals = FunctionEvaluations(x, z)
pwl = approx(vals, Convex(), Heuristic(); optimizer=Xpress.Optimizer, planes=7, strict=:below)
PWL.plotconv2D(vals,pwl)

# Uniform sampling for various convex functions in 2D

# Quadratic
f(x) = x[1]^2 + x[2]^2
vals = PWL.sample_uniform(f, [(-1,1), (-1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), optimizer=Xpress.Optimizer, planes=10)
PWL.plotconv2D(vals, pwl)

# Log sum exp
h(x) = log(exp(x[1]) + exp(x[2]))
vals = PWL.sample_uniform(h, [(-1,1), (-1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), optimizer=Xpress.Optimizer, planes=5)
PWL.plotconv2D(vals, pwl)

# Quadratic over linear
g(x) = x[1]^2 / x[2]
vals = PWL.sample_uniform(g, [(-1,1), (0.1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), optimizer=Xpress.Optimizer, planes=10)
PWL.plotconv2D(vals, pwl)

# Geometric mean
f(x) = -sqrt(x[1]*x[2])
vals = PWL.sample_uniform(f, [(0.1,1), (0.1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), optimizer=Xpress.Optimizer, planes=10)
PWL.plotconv2D(vals, pwl)

# Non differentiable
f(x) = max(x[1]^2,x[2]^2)
vals = PWL.sample_uniform(f, [(-1,1), (-1,1)], 10)
pwl = approx(vals, Convex(), Heuristic(), optimizer=Xpress.Optimizer, planes=8, pen=:max)
PWL.plotconv2D(vals, pwl)


# Uniform sampling 3D
h(x) = log(exp(x[1]) + exp(x[2]) + exp(x[3]))
vals = PWL.sample_uniform(h, [(-5,5), (-5,5), (-5,5)], 10)
pwl = approx(vals, Convex(), Heuristic(), optimizer=Xpress.Optimizer, planes=8, pen=:rms)