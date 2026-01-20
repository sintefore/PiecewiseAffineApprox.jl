function calculate_X(x1, x2, x3)
    X = hcat(
        repeat(x1, inner = [length(x2) * length(x3)]),
        repeat(x2, inner = [length(x3)], outer = [length(x1)]),
        repeat(x3, outer = [length(x1) * length(x2)]),
    )
    valid_indices = X[:, 1] .^ 2 .> X[:, 2] .^ 2
    X = X[valid_indices, :]
    return X
end

function calculate_z(x1, x2, x3)
    M1 = 16.042
    M2 = 2.016
    M3 = 28.96
    constant = 0.1375

    return sqrt(constant * (x1^2 - x2^2) * (M3 / (M1 * (1 - x3) + M2 * x3)))
end

function compute_concave(strict = :outer)
    x1 = [i for i ∈ 145:5:175]
    x2 = [i for i ∈ 145:5:175]
    x3 = [i for i ∈ 0.00:0.01:0.2]

    X = calculate_X(x1, x2, x3)
    z = calculate_z.(X[:, 1], X[:, 2], X[:, 3])

    pwa = approx(
        FunctionEvaluations(collect(zip(X[:, 1], X[:, 2], X[:, 3])), z),
        Concave(),
        Cluster(; optimizer, planes = 10, strict = strict, metric = :l1),
    )
    return (; X, z, pwa)
end

@testset "Outer approx for concave functions" begin
    _, _, pwa = compute_concave(:outer)
    p_in = 175
    fraction = 0.05
    total = 0
    violations = 0
    for p_out ∈ 140:175
        actual = PiecewiseAffineApprox.evaluate(pwa, (p_in, p_out, fraction))
        expected = calculate_z(p_in, p_out, fraction)
        total += 1
        if actual < expected
            violations += 1
        end
    end

    # Test that at least 90% of points satisfy the outer approximation
    @test violations / total ≤ 0.10
end
