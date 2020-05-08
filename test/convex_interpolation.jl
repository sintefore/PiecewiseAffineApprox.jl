# Points to approximate
points = [i => i^2 for i in -2:0.1:2]

# Use Cbc without output as the solver
Opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

# Test interpolation routine
prs =  PiecewiseLinearApprox.interpolatepw(points, Opt; nseg=5);

@test length(prs) == 6
@test prs[1][2] == prs[6][2]
@test prs[2][2] == prs[5][2]
@test prs[3][2] == prs[4][2]

# Test convexification
rpoints = [i => i^2 - 0.5 + rand() for i in -2:0.1:2]
@test !PiecewiseLinearApprox.isconvex(rpoints)
cpoints = PiecewiseLinearApprox.convexify(rpoints, Opt)
@test PiecewiseLinearApprox.isconvex(cpoints)

# Test with constraints added to existing model
d = first.(points)
fd = last.(points)
m = Model()
@variable(m, x)
y = PiecewiseLinearApprox.convex_pwlinear_interpolate(m, x, d, fd, Opt; nseg=11)
@objective(m, Min, y)
set_optimizer(m,Opt)
@constraint(m, x >= 0.9)
optimize!(m)
@test isapprox(value(m[:y]), 0.81, atol=0.05)
