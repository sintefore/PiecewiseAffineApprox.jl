@testset "Files in the format written with JSON3 are read" begin
    # JSON3 was supported until version 0.6. It wrote the same format as the
    # JSON extension; these are files as it wrote them.
    tmp = mktempdir()
    function parsed(name, text, T)
        fn = joinpath(tmp, name)
        write(fn, text)
        return JSON.parsefile(fn, T)
    end
    p = parsed(
        "convex.json",
        """{"curvature":"convex","planes":[{"α":[1.0,2.0],"β":1.0},{"α":[-3.0,1.0],"β":-2.0}],"lower":[0.0,-1.0],"upper":[1.0,2.0]}""",
        PWAFunc,
    )
    @test p == PWAFunc{Convex,2}(
        [Plane((1.0, 2.0), 1.0), Plane((-3.0, 1.0), -2.0)];
        domain = [(0.0, 1.0), (-1.0, 2.0)],
    )
    q = parsed(
        "concave.json",
        """{"curvature":"concave","planes":[{"α":[0.5],"β":0.25}],"lower":[],"upper":[]}""",
        PWAFunc,
    )
    @test q == PWAFunc{Concave,1}([Plane((0.5,), 0.25)])
    # Files from before the domain was stored
    q = parsed(
        "old.json",
        """{"curvature":"concave","planes":[{"α":[0.5],"β":0.25}]}""",
        PWAFunc,
    )
    @test q == PWAFunc{Concave,1}([Plane((0.5,), 0.25)])
    f = parsed(
        "dc.json",
        """{"curvature":"dc","convex":[{"α":[1.0,2.0],"β":1.0}],"concave":[{"α":[0.5,0.5],"β":0.0}],"lower":[],"upper":[]}""",
        DCFunc,
    )
    @test f == DCFunc(
        PWAFunc{Convex,2}([Plane((1.0, 2.0), 1.0)]),
        PWAFunc{Concave,2}([Plane((0.5, 0.5), 0.0)]),
    )
end

@testset "JSON Serialization" begin
    x = 1:10
    y = x .^ 2
    tmp = mktempdir()
    pwa = PiecewiseAffineApprox.approx(
        FunctionEvaluations(Tuple.(x), y),
        Convex(),
        MILP(; optimizer, planes = 5),
    )
    pwa_concave = PiecewiseAffineApprox.approx(
        FunctionEvaluations(Tuple.(x), -y),
        Concave(),
        MILP(; optimizer, planes = 5),
    )
    for (i, p) ∈ enumerate((pwa, pwa_concave))
        fn = joinpath(tmp, "test_$i.json")
        JSON.json(fn, p)
        read_back = JSON.parsefile(fn, PWAFunc)
        @test isa(read_back, PWAFunc)
        @test length(read_back.planes) == 5
    end
end
