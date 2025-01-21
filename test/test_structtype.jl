@testset "Serialization" begin
    x = 1:10
    y = x .^ 2
    tmp = mktempdir()
    fn = joinpath(tmp, "test.json")
    pwa = PiecewiseAffineApprox.approx(
        FunctionEvaluations(Tuple.(x), y),
        Convex(),
        MILP(; optimizer = HiGHS.Optimizer, planes = 5),
    )
    pwa_concave = PiecewiseAffineApprox.approx(
        FunctionEvaluations(Tuple.(x), -y),
        Concave(),
        MILP(; optimizer = HiGHS.Optimizer, planes = 5),
    )
    for p ∈ (pwa, pwa_concave)
        JSON3.write(fn, p)
        read_back = JSON3.read(fn, PWAFunc)
        @test isa(read_back, PWAFunc)
        @test length(read_back.planes) == 5
    end
end
