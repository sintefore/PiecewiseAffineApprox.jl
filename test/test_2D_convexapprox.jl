PWL = PiecewiseLinearApprox

xg = [i for i in -1:0.5:1]
yg = [j for j in -1:0.5:1]

X = [repeat(xg, inner=[size(yg,1)]) repeat(yg, outer=[size(xg,1)])]

z = X[:,1].^2 + X[:,2].^2

# Use Xpress without output as the solver
Opt = optimizer_with_attributes(Xpress.Optimizer)

np = 4

pwl1 = convex_linearization(X, z, Opt; nplanes=np, dimensions=2, strict=:above, pen=:l2)
@test length(pwl1.a) == np
@test isapprox(PWL.evaluate(pwl1, (0.5, 0.5)), 0.5, atol=0.1)
