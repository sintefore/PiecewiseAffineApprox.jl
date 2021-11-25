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

# Test with constraints added using already existing y-variable
m = Model()
@variable(m, xvar)
@variable(m, yvar)
@variable(m, test_f)

t = (xvar, yvar)

y = PiecewiseLinearApprox.convex_pwlinear(m,t,X,z,Opt;nplanes=np, dimensions=2, strict=:above, pen=:l2, z=test_f)
@objective(m, Min, y)
set_optimizer(m,Opt)
@constraint(m, xvar >= 0.3)
optimize!(m)

xval = JuMP.value(m[:xvar])
yval = JuMP.value(m[:yvar])
fval = JuMP.value(m[:test_f])

@test isapprox(value(m[:test_f]), -0.15, atol=0.015)







