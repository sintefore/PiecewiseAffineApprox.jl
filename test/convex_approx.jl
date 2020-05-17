PWL = PiecewiseLinearApprox

x = [i for i in -1:0.1:1]
z = x.^2

# Use Cbc without output as the solver
Opt = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

pwl1 = convex_linearization(x, z, Opt; nseg=5)
@test length(pwl1.c) == 5
@test length(pwl1.x) == 6
@test issorted(pwl1.c)
@test PWL.evaluate(pwl1, 0.4) ≈ 0.18


pwl2 = convex_linearization(x -> x^2, -1, 1, Opt; nseg=5)
@test length(pwl2.c) == 5
@test length(pwl2.x) == 6
@test PWL.evaluate(pwl2, 0.4) ≈ 0.164

pwl3 = convex_linearization(x -> x^2, -1, 1, Opt; nseg=4, strict=true)
@test length(pwl2.c) == 5
@test length(pwl2.x) == 6
@test PWL.evaluate(pwl2, 0.4) ≈ 0.144

pwl4 = convex_linearization(x -> x^2, -1, 1, Opt; nseg=5, method=:ipol)
@test PWL.evaluate(pwl2, 0.4) ≈ 0.2


