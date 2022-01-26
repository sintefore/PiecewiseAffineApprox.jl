xg = [i for i in -1:0.5:1]
yg = [j for j in -1:0.5:1]

X = [repeat(xg, inner=[size(yg,1)]) repeat(yg, outer=[size(xg,1)])]

z = X[:,1].^2 + X[:,2].^2
z_concave = z.*-1

np = 5

pwl1 = approx(FunctionEvaluations(PWL.mat2tuples(X), z), Convex(), Optimized() ;optimizer=quadopt(1.7), planes=np, dimensions=2, strict=:above, pen=:l2)
pwl2 = approx(FunctionEvaluations(PWL.mat2tuples(X), z_concave), Concave(), Optimized(); optimizer=quadopt(1.7), planes=np, dimensions=2, strict=:above, pen=:l2)

# @test length(pwl1.a) == np
@test isapprox(PWL.evaluate(pwl1, (0.5, 0.5)), 0.5, atol=0.1)
@test isapprox(PWL.evaluate(pwl1, (-0.3, 0.4)), -PWL.evaluate(pwl2, (-0.3, 0.4)), atol=0.01)


# Test with constraints added using already existing variables (tuple with xvar and yvar)
m = Model()
@variable(m, xvar)
@variable(m, yvar)
@variable(m, test_f)

tuple_var = (xvar, yvar)

# TODO: FIX/Update (Optimization fails with most recent Xpress release)
#y = PiecewiseLinearApprox.convex_pwlinear(m,tuple_var,X,z,quadopt(1.7);planes=4, dimensions=2, strict=:above, pen=:l2, z=test_f)

# @objective(m, Min, y)
# set_optimizer(m,quadopt())
# @constraint(m, xvar >= 0.3)
# optimize!(m)

# xval = JuMP.value(m[:xvar])
# yval = JuMP.value(m[:yvar])
# fval = JuMP.value(m[:test_f])

@test_broken isapprox(value(m[:test_f]), -0.36, atol=0.04)

# m = Model()
# @variable(m, xvar_conc)
# @variable(m, yvar_conc)
# @variable(m, f_conc)

# tuple_var_conc = (xvar_conc, yvar_conc)

# y_concave = PWL.pwlinear(m,tuple_var_conc,FunctionEvaluations(PWL.mat2tuples(X),z_concave),Concave(),Optimized();optimizer=quadopt(1.7),planes=np, dimensions=2, strict=:above, pen=:l2, z=f_conc)

# @objective(m, Max, y_concave)
# set_optimizer(m,quadopt(2.2))
# @constraint(m, xvar_conc >= 0.3)
# optimize!(m)

# xval_conc = JuMP.value(m[:xvar_conc])
# yval_conc = JuMP.value(m[:yvar_conc])
# fval_conc = JuMP.value(m[:f_conc])

@test_broken isapprox(fval, -fval_conc, atol=0.01)
@test_broken isapprox(xval, xval_conc, atol=0.01)
@test_broken isapprox(yval, yval_conc, atol=0.01)
