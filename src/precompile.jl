# Workload for precompilation, covering the parts that do not need an
# optimizer: fitting with Interpol and Cluster (convex, concave and DC), and
# adding the fitted functions to a JuMP model with the different formulations.
# The algorithms that solve optimization problems are compiled together with
# the optimizer on first use.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    rng = StableRNG(1)
    @compile_workload begin
        # One-dimensional fits, from a function and from data
        f₁ = approx(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 4))
        g₁ = approx(x -> -x[1]^2, [(-1, 1)], Concave(), Interpol(planes = 4))
        dc₁ = approx(sin, [(0, 2π)], DC(), Interpol(planes = 3))
        x = collect(range(-1, 1; length = 20))
        data = FunctionEvaluations(tuple.(x), x .^ 2)
        r = fit(data, Convex(), Cluster(; planes = 3, rng))
        approximation_error(r.func, data, :max)
        sprint(show, r)

        # Two-dimensional fits
        box = [(-1, 1), (-1, 1)]
        f₂ = approx((x, y) -> x^2 + y^2, box, Convex(), Cluster(; planes = 4, rng))
        g₂ = approx((x, y) -> -x^2 - y^2, box, Concave(), Cluster(; planes = 4, rng))
        dc₂ = approx((x, y) -> x^2 - y^2, box, DC(), Cluster(; planes = 3, rng))

        evaluate(f₁, (0.5,))
        evaluate(dc₁, (0.5,))
        evaluate(f₂, (0.5, 0.5))
        evaluate(dc₂, (0.5, 0.5))

        # Adding the functions to a model
        m = Model()
        @variable(m, -1 ≤ u ≤ 1)
        @variable(m, 0 ≤ w ≤ 2π)
        @variable(m, -1 ≤ v[1:2] ≤ 1)
        f₁(u)
        g₁(u)
        for F ∈ (
            PsiFormulation(),
            LambdaFormulation(),
            DeltaFormulation(),
            PhiFormulation(),
        )
            f₁(u; formulation = F)
        end
        f₂(v)
        g₂(v)
        f₂(v; formulation = LambdaFormulation())
        dc₂(v)
        dc₁(w)
        for F ∈ (PsiFormulation(), DeltaFormulation(), PhiFormulation())
            dc₁(w; formulation = F)
        end
        z, record = add_pwaformulation(m, u, f₁)
        delete(m, record)
    end
end
