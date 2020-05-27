# Points to approximate
points = [i=>i^2 for i in 0:0.1:1]

# Use Cbc without output as the solver
Opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)


# Test with constraints added to existing model
m = Model()
@variable(m, x)
y = PiecewiseLinearApprox.convex_pwlinear(m,x,[i.first for i ∈ points],[i.second for i ∈ points],Opt; nseg=5)
@objective(m, Min, y)
set_optimizer(m,Opt)
@constraint(m, x >= 0.3)
optimize!(m)
@test isapprox(value(y), 0.09, atol=0.01)
@constraint(m, x>= 0.9)
optimize!(m)
@test isapprox(value(y), 0.81, atol=0.01)

# Test with constraints added using already existing y-variable
m = Model()
@variable(m, x)
@variable(m, test_y)
y = PiecewiseLinearApprox.convex_pwlinear(m,x,[i.first for i ∈ points],[i.second for i ∈ points],Opt;nseg=5,z=test_y)
@objective(m, Min, y)
set_optimizer(m,Opt)
@constraint(m, x >= 0.3)
optimize!(m)
@test isapprox(value(m[:test_y]), 0.09, atol=0.01)
