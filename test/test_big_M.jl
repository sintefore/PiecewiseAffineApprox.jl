# big-M for a 1D function
x = collect(range(-1, 1; length = 10))
z = x .^ 2

M₁ = PWA.conv_linear_big_M(x, z)
M₁⁺ = PWA.linear_big_M(x, z)

pwl₁ = approx(
    FunctionEvaluations(Tuple.(x), z),
    Convex(),
    Optimized(optimizer, planes = 5, bigM = :conv_linear_big_M),
)
pwl₁⁺ = approx(
    FunctionEvaluations(Tuple.(x), z),
    Convex(),
    Optimized(optimizer, planes = 5, bigM = :linear_big_M),
)

@test isapprox(PWA.evaluate(pwl₁, 0.4), 0.16, atol = 0.015)
@test isapprox(PWA.evaluate(pwl₁⁺, 0.4), 0.16, atol = 0.015)

@test isless(M₁⁺, M₁)  # new big-M is smaller than previous but gives a better approximation

# big-M for a 2D function
xg = collect(range(-1, 1; length = 4))
yg = collect(range(-1, 1; length = 4))

X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]

z = X[:, 1] .^ 2 + X[:, 2] .^ 2
np = 5

pwl₂ = approx(
    FunctionEvaluations(PWA.mat2tuples(X), z),
    Convex(),
    Optimized(
        optimizer = optimizer,
        planes = np,
        dimensions = 2,
        strict = :above,
        pen = :l1,
    )bigM = :conv_linear_big_M_ND,
) # old estimate for big-M

pwl₂⁺ = approx(
    FunctionEvaluations(PWA.mat2tuples(X), z),
    Convex(),
    Optimized(
        optimizer = optimizer,
        planes = np,
        dimensions = 2,
        strict = :above,
        pen = :l1,
    )bigM = :linear_big_M,
)

@test isapprox(
    PWA.evaluate(pwl₂, (0.5, 0.5)),
    PWA.evaluate(pwl₂⁺, (0.5, 0.5)),
    atol = 0.01,
)

M₂⁺ = PWA.linear_big_M(x, z)
M₂ = PWA.conv_linear_big_M_ND(x, z)  # old estimate for big-M

@test isless(M₂, M₂⁺)

@testset "Duplicate inputs" begin
    x = [302.3633312955, 300.6680553405, 274.281686739, 274.281686739]
    z = [-184.9368247, -183.9371662, -167.94263, -167.94263]

    @test_broken !isnan(PWA.linear_big_M(x, z))
    @test PWA.linear_big_M(unique(x), unique(z)) ≈ 17.02 atol = 0.01
end
