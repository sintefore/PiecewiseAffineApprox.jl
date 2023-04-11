
x = collect(range(-1, 1; length = 10))
z = x .^ 2

pwl1 = approx(
    FunctionEvaluations(Tuple.(x), z),
    Convex(),
    Optimized(),
    ;
    optimizer,
    pen = :l1,
    planes = 5,
)
@test length(pwl1.planes) == 5
@test issorted((p.α for p ∈ pwl1.planes))
@test isapprox(PWA.evaluate(pwl1, 0.4), 0.16, atol = 0.035)

pwl2 = approx(
    x -> x[1]^2,
    [(-1, 1)],
    Convex(),
    Optimized();
    optimizer,
    pen = :l1,
    planes = 5,
)
@test length(pwl2.planes) == 5
@test isapprox(PWA.evaluate(pwl2, 0.4), 0.1733, atol = 0.015)

pwl3 = approx(
    x -> x[1]^2,
    [(-1, 1)],
    Convex(),
    Optimized();
    optimizer,
    pen = :l1,
    planes = 5,
    strict = true,
)
@test length(pwl2.planes) == 5
@test isapprox(PWA.evaluate(pwl3, 0.4), 0.16, atol = 0.03)

# pwl4 = convex_linearization(x -> x^2, -1, 1, optimizer; planes=5, method=:ipol)
pwl4 =
    approx(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(); optimizer, planes = 10)
@test isapprox(PWA.evaluate(pwl4, 0.4), 0.16, atol = 0.015)

#Check approximation for flat function with releatively high values:

I = 10
xmat = 2 * rand(rng, 2, I) .- 1
x = [Tuple(xmat[:, i]) for i ∈ 1:size(xmat, 2)]
z = [10_000 + 0.001 * p[1]^2 + 0.001 * p[2]^2 for p ∈ x]
vals = FunctionEvaluations(x, z)
pwl = approx(
    vals,
    Convex(),
    Optimized();
    optimizer = optimizer,
    planes = 4,
    pen = :l1,
)
@test PWA.evaluate(pwl, (0, 0)) ≈ 10_000.0 atol = 0.1
