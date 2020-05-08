# Points to approximate
points = [i=>i^2 for i in 0:0.1:1]

# Use Cbc without output as the solver
Opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

# Check that generated points are as expected
prs =  PiecewiseLinearApprox.optimize_points(points,Opt;numpoints=3);

sort!(prs) # Sort to avoid symmetrical solutions
@test prs[1][1] ≈ -0.47115384615384626
@test prs[1][2] ≈  1.4134615384615385

# Test with constraints added to existing model
m = Model()
@variable(m, x)
y = PiecewiseLinearApprox.convexlinearization(m,x,[i.first for i ∈ points],[i.second for i ∈ points],opt=Opt,K=5)
@objective(m, Min, y)
set_optimizer(m,Opt)
@constraint(m, x >= 0.3)
optimize!(m)
@test isapprox(value(m[:y]), 0.09, atol=0.01)
@constraint(m, x>= 0.9)
optimize!(m)
@test isapprox(value(m[:y]), 0.81, atol=0.01)

# Test with constraints added using already existing y-variable
m = Model()
@variable(m, x)
@variable(m, test_y)
y = PiecewiseLinearApprox.convexlinearization(m,x,[i.first for i ∈ points],[i.second for i ∈ points],opt=Opt,K=5,y=test_y)
@objective(m, Min, y)
set_optimizer(m,Opt)
@constraint(m, x >= 0.3)
optimize!(m)
@test isapprox(value(m[:test_y]), 0.09, atol=0.01)
