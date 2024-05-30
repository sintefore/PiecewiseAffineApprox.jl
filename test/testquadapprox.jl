@testset "Quad Approx" begin
    # Points to approximate
    points = [i => i^2 for i ∈ 0:0.1:1]

    @testset "L1" begin
        # Test with constraints added to existing model
        m = Model()
        @variable(m, x)
        y = PWA.pwaffine(
            m,
            tuple(x),
            FunctionEvaluations(
                [tuple(i.first) for i ∈ points],
                [i.second for i ∈ points],
            ),
            Convex(),
            Optimized(optimizer = optimizer, pen = :l1, planes = 5),
        )
        @objective(m, Min, y)
        set_optimizer(m, optimizer)
        @constraint(m, x >= 0.3)
        optimize!(m)
        @test isapprox(value(y), 0.09, atol = 0.015)
        @constraint(m, x >= 0.9)
        optimize!(m)
        @test isapprox(value(y), 0.81, atol = 0.1)

        # Test with constraints added using already existing y-variable
        m = Model()
        @variable(m, x)
        @variable(m, test_y)
        y = PWA.pwaffine(
            m,
            tuple(x),
            FunctionEvaluations(
                [tuple(i.first) for i ∈ points],
                [i.second for i ∈ points],
            ),
            Convex(),
            Optimized(optimizer = optimizer, pen = :l1, planes = 5),
            ;
            z = test_y,
        )
        @objective(m, Min, y)
        set_optimizer(m, optimizer)
        @constraint(m, x >= 0.3)
        optimize!(m)
        @test isapprox(value(m[:test_y]), 0.09, atol = 0.015)
    end

    @testset "L2" begin
        # Test with constraints added to existing model
        for met in [
            Heuristic(optimizer = qp_optimizer, pen = :l2, planes = 5),
            Optimized(optimizer = qp_optimizer, pen = :l2, planes = 5),
        ]
            m = Model()
            @variable(m, x)
            y = PWA.pwaffine(
                m,
                tuple(x),
                FunctionEvaluations(
                    [tuple(i.first) for i ∈ points],
                    [i.second for i ∈ points],
                ),
                Convex(),
                met;
            )
            @objective(m, Min, y)
            set_optimizer(m, optimizer)
            @constraint(m, x >= 0.3)
            optimize!(m)
            @test isapprox(value(y), 0.09, atol = 0.015)
            @constraint(m, x >= 0.9)
            optimize!(m)
            @test isapprox(value(y), 0.81, atol = 0.1)
        end
    end
end
