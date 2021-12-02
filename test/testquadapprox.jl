# Points to approximate
points = [i=>i^2 for i in 0:0.1:1]

# Test with constraints added to existing model
m = Model()
@variable(m, x)
y = PiecewiseLinearApprox.pwlinear(m, tuple(x),
    FunctionEvaluations([tuple(i.first) for i ∈ points], [i.second for i ∈ points]),
     Convex(), Optimized(); optimizer, nseg=5)
@objective(m, Min, y)
set_optimizer(m, optimizer)
@constraint(m, x >= 0.3)
optimize!(m)
@test isapprox(value(y), 0.09, atol=0.015)
@constraint(m, x>= 0.9)
optimize!(m)
@test isapprox(value(y), 0.81, atol=0.1)

# Test with constraints added using already existing y-variable
m = Model()
@variable(m, x)
@variable(m, test_y)
y = PiecewiseLinearApprox.pwlinear(m, tuple(x),
    FunctionEvaluations([tuple(i.first) for i ∈ points],[i.second for i ∈ points]),
    Convex(), Optimized(); optimizer, nseg=5, z=test_y)
@objective(m, Min, y)
set_optimizer(m, optimizer)
@constraint(m, x >= 0.3)
optimize!(m)
@test isapprox(value(m[:test_y]), 0.09, atol=0.015)
