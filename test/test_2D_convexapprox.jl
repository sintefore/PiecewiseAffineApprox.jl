PWL = PiecewiseLinearApprox

xg = [i for i in -1:0.025:1]
yg = [j for j in -1:0.025:1]
x = [xg, yg]

X = ndgrid(x[1], x[2])
z = X[1].^2 + X[2].^2

# Use Xpress without output as the solver
Opt = optimizer_with_attributes(Xpress.Optimizer)

pwl1 = convex_linearization(x, z, Opt; nplanes=4, dimensions=2)
@test length(pwl1.c) == 8
@test isapprox(PWL.evaluate(pwl1, 0.25, 0.25), 0.5, atol=0.01)

