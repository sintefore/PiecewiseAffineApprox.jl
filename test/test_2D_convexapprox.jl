PWL = PiecewiseLinearApprox

xg = [i for i in -1:0.1:1]
yg = [j for j in -1:0.1:1]
x = [xg, yg]

X = ndgrid(x[1], x[2])
z = X[1].^2 + X[2].^2 + X[1].*X[2]

# Use Xpress without output as the solver
Opt = optimizer_with_attributes(Xpress.Optimizer)

np = 8

pwl1 = convex_linearization(x, z, Opt; nplanes=np, dimensions=2)
@test length(pwl1.c) == np
@test isapprox(PWL.evaluate(pwl1, 0.5, 0.5), 0.75, atol=0.01)
