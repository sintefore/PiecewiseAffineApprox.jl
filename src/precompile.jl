# Workload for precompilation, covering the parts that do not need an
# optimizer: fitting with Interpol and Cluster (convex, concave and DC), and
# adding the fitted functions to a JuMP model with the default formulations.
# The algorithms that solve optimization problems, and the other formulations,
# are compiled on first use. Most of the time goes to the first fit with each
# algorithm, which compiles all its options, as they are runtime values.

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    rng = StableRNG(1)
    @compile_workload begin
        # One-dimensional fits
        f₁ = approx(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 4))
        g₁ = approx(x -> -x[1]^2, [(-1, 1)], Concave(), Interpol(planes = 4))
        dc₁ = approx(sin, [(0, 2π)], DC(), Interpol(planes = 3))

        # Two-dimensional fits
        box = [(-1, 1), (-1, 1)]
        f₂ = approx((x, y) -> x^2 + y^2, box, Convex(), Cluster(; planes = 4, rng))
        g₂ = approx((x, y) -> -x^2 - y^2, box, Concave(), Cluster(; planes = 4, rng))
        dc₂ = approx((x, y) -> x^2 - y^2, box, DC(), Cluster(; planes = 3, rng))

        # A fit from data, which shares most of its code with the fits above
        x = collect(range(-1, 1; length = 20))
        data = FunctionEvaluations(tuple.(x), x .^ 2)
        r = fit(data, Convex(), Cluster(; planes = 3, rng))
        approximation_error(r.func, data, :max)
        sprint(show, r)

        evaluate(f₁, (0.5,))
        evaluate(dc₁, (0.5,))
        evaluate(f₂, (0.5, 0.5))
        evaluate(dc₂, (0.5, 0.5))

        # Adding the functions to a model with the default formulations
        m = Model()
        @variable(m, -1 ≤ u ≤ 1)
        @variable(m, 0 ≤ w ≤ 2π)
        @variable(m, -1 ≤ v[1:2] ≤ 1)
        f₁(u)
        g₁(u)
        dc₁(w)
        f₂(v)
        dc₂(v)
    end
end
