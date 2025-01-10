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
    JSON3.write(fn, pwa)
    read_back = JSON3.read(fn, PWAFunc)
    @test isa(read_back, PWAFunc)
    @test length(read_back.planes) == 5
end
