PWL = PiecewiseLinearApprox

x = [i for i in -1:0.1:1]
z = x.^2

# Use Cbc without output as the solver
Opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

pwl1 = convex_linearization(x, z, Opt; nseg=5)
@test length(pwl1.c) == 5
@test length(pwl1.x) == 6
@test issorted(pwl1.c)
@test isapprox(PWL.evaluate(pwl1, 0.4), 0.16, atol=0.001)


pwl2 = convex_linearization(x -> x^2, -1, 1, Opt; nseg=5)
@test length(pwl2.c) == 5
@test length(pwl2.x) == 6
@test isapprox(PWL.evaluate(pwl2, 0.4), 0.1733, atol=0.001)

pwl3 = convex_linearization(x -> x^2, -1, 1, Opt; nseg=4, strict=true)
@test length(pwl2.c) == 5
@test length(pwl2.x) == 6
@test isapprox(PWL.evaluate(pwl3, 0.4), 0.0666, atol=0.001)

pwl4 = convex_linearization(x -> x^2, -1, 1, Opt; nseg=5, method=:ipol)
@test isapprox(PWL.evaluate(pwl4, 0.4), 0.2, atol=0.001)


