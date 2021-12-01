PWL = PiecewiseLinearApprox

x = [i for i in -1:0.1:1]
z = x.^2

pwl1 = approx(FunctionEvaluations(Tuple.(x), z), Convex(), Optimized(), ;optimizer, nplanes=5)
@test length(pwl1.planes) == 5
@test issorted((p.α for p in pwl1.planes))
@test isapprox(PWL.evaluate(pwl1, 0.4), 0.16, atol=0.03)

pwl2 = approx(x -> x^2, -1, 1, Convex(), Optimized(); optimizer, nplanes=5)
@test length(pwl2.planes) == 5
@test isapprox(PWL.evaluate(pwl2, 0.4), 0.1733, atol=0.015)

pwl3 = approx(x -> x^2, -1, 1, Convex(), Optimized(); optimizer, nplanes=5, strict=true)
@test length(pwl2.planes) == 5
@test isapprox(PWL.evaluate(pwl3, 0.4), 0.16, atol=0.015)

# pwl4 = convex_linearization(x -> x^2, -1, 1, optimizer; nseg=5, method=:ipol)
pwl4 = approx(x -> x^2, -1, 1, Convex(), Interpol(); optimizer, nseg=5)
@test isapprox(PWL.evaluate(pwl4, 0.4), 0.2, atol=0.015)


