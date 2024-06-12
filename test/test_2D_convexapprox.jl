@testset "2D Convex Approx" begin
    xg = [i for i ∈ -1:0.5:1]
    yg = [j for j ∈ -1:0.5:1]

    X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]

    z = X[:, 1] .^ 2 + X[:, 2] .^ 2
    z_concave = z .* -1

    np = 17

    function mat2tuples(x::Matrix)
        return collect(Tuple(x'[:, i]) for i ∈ 1:size(x', 2))
    end

    pwa1 = approx(
        FunctionEvaluations(mat2tuples(X), z),
        Convex(),
        MILP(optimizer = optimizer, planes = np, strict = :outer, pen = :l1),
    )
    pwa2 = approx(
        FunctionEvaluations(mat2tuples(X), z_concave),
        Concave(),
        MILP(optimizer = optimizer, planes = np, strict = :outer, pen = :l1),
    )

    # @test length(pwa1.a) == np
    @test isapprox(PWA.evaluate(pwa1, (0.5, 0.5)), 0.5, atol = 0.1)
    @test isapprox(
        PWA.evaluate(pwa1, (-0.3, 0.4)),
        -PWA.evaluate(pwa2, (-0.3, 0.4)),
        atol = 0.01,
    )

    # Test with constraints added using already existing variables (tuple with xvar and yvar)
    m = Model()
    @variable(m, xvar)
    @variable(m, yvar)
    @variable(m, test_f)

    tuple_var = (xvar, yvar)

    y = PWA.pwaffine(
        m,
        tuple_var,
        approx(
            FunctionEvaluations(mat2tuples(X), z),
            Convex(),
            MILP(
                optimizer = optimizer,
                planes = np,
                strict = :outer,
                pen = :l1,
            ),
        );
        z = test_f,
    )

    @objective(m, Min, y)
    set_optimizer(m, optimizer)
    @constraint(m, xvar == √0.5)
    @constraint(m, yvar == √0.5)
    optimize!(m)

    xval = JuMP.value(m[:xvar])
    yval = JuMP.value(m[:yvar])
    fval = JuMP.value(m[:test_f])

    @test isapprox(value(m[:test_f]), 1.0, rtol = 0.12)

    m = Model()
    @variable(m, xvar_conc)
    @variable(m, yvar_conc)
    @variable(m, f_conc)

    tuple_var_conc = (xvar_conc, yvar_conc)

    y_concave = PWA.pwaffine(
        m,
        tuple_var_conc,
        FunctionEvaluations(mat2tuples(X), z_concave),
        Concave(),
        MILP(optimizer = optimizer, planes = np, strict = :outer, pen = :l1);
        z = f_conc,
    )

    @objective(m, Max, y_concave)
    set_optimizer(m, optimizer)
    @constraint(m, xvar_conc == √0.5)
    @constraint(m, yvar_conc == √0.5)
    optimize!(m)

    xval_conc = JuMP.value(m[:xvar_conc])
    yval_conc = JuMP.value(m[:yvar_conc])
    fval_conc = JuMP.value(m[:f_conc])

    @test isapprox(fval, -fval_conc, atol = 0.01)
    @test isapprox(xval, xval_conc, atol = 0.01)
    @test isapprox(yval, yval_conc, atol = 0.01)
end
