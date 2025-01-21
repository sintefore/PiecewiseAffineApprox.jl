@testset "Serialization" begin
    x = 1:10
    y = x .^ 2
    tmp = mktempdir()
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
    for (i,p) ∈ enumerate((pwa, pwa_concave))
        fn = joinpath(tmp, "test_$i.json")
        JSON3.write(fn, p)
        read_back = JSON3.read(fn, PWAFunc)
        @test isa(read_back, PWAFunc)
        @test length(read_back.planes) == 5
    end
end
