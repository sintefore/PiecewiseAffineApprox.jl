@testset "Linear optimization 1D" begin

    # 8 plane approximation of z = (x-1)^2 on [-1, 3]
    pwa = JSON.parsefile(joinpath(@__DIR__, "pwa_ex", "pwa_quad_1d_8planes.json"), PWAFunc)
    v = PWA.vertices_by_subsets(pwa)
    x_min, z_min = argmin(v -> v[2], v)
    # Test all 4 possible formulations in 1D, with bounds taken from the variable
    for f in [λ_Formulation, Ψ_Formulation, Δ_Formulation, Φ_Formulation]
        m = JuMP.Model(HiGHS.Optimizer)
        set_silent(m)
        @variable(m, -1 ≤ x ≤ 3)
        @variable(m, z)
        @objective(m, Min, z)
        pwaffine(m, x, pwa; z = z, formulation = f())
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value(x) ≈ x_min[1]
        @test objective_value(m) ≈ z_min

        # Restrict x-values
        c = @constraint(m, x >= 1.5)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value(x) ≈ 1.5
        @test objective_value(m) ≈ evaluate(pwa, 1.5)

        # Beyond the last breakpoint (≈ 2.71), but inside the bounds
        delete(m, c)
        @constraint(m, x >= 2.9)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value(x) ≈ 2.9
        @test objective_value(m) ≈ evaluate(pwa, 2.9)
    end
end

@testset "Linear optimization 2D" begin
    # Minimize approximation of x₁² + x₂² on [-1, 1]² subject to x₁ + 3x₂ ≥ 2
    # (the exact function has its minimum at (0.2, 0.6) with objective 0.4)
    pwa = JSON.parsefile(joinpath(@__DIR__, "pwa_ex", "pwa_quad_2d_12planes.json"), PWAFunc)
    v = PWA.vertices_by_subsets(pwa)
    x_min, z_min = argmin(v -> v[2], v)

    # Test both possible formulations in 2D, with explicit bounds
    for f in [λ_Formulation, Ψ_Formulation]
        m = JuMP.Model(HiGHS.Optimizer)
        set_silent(m)
        @variable(m, x[1:2])
        @variable(m, z)
        pwaffine(m, x, pwa; z = z, formulation = f(), bounds = [(-1, 1), (-1, 1)])
        @objective(m, Min, z)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value.(x) ≈ x_min
        @test objective_value(m) ≈ z_min

        # Restrict x-values
        @constraint(m, x[1] + 3*x[2] ≥ 2)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value.(x) ≈ [0.29, 0.57] atol=1e-6
        @test objective_value(m) ≈ 0.33739 atol=1e-5
        @test value(z) ≈ evaluate(pwa, value.(x)) atol=1e-5
    end
end

# Solve min/max z for the pwa function, with x restricted to [xlo, ∞)
function _solve_pwa(pwa::PWAFunc{C}, f, xlo; bounds) where {C}
    m = JuMP.Model(HiGHS.Optimizer)
    set_silent(m)
    @variable(m, x ≥ xlo)
    z = pwaffine(m, x, pwa; formulation = f(), bounds = bounds)
    C == Convex ? @objective(m, Min, z) : @objective(m, Max, z)
    optimize!(m)
    return termination_status(m), value(x), value(z)
end

@testset "Linear optimization formulations - edge cases" begin
    formulations = [λ_Formulation, Ψ_Formulation, Δ_Formulation, Φ_Formulation]
    # Tangents to x² at ±0.25 and ±0.75, with breakpoints at -0.5, 0 and 0.5
    planes = [Plane((2t,), -t^2) for t in (-0.75, -0.25, 0.25, 0.75)]
    cases = [
        # Between the last breakpoint and the upper bound
        (PWAFunc{Convex,1}(planes), 0.9, 0.7875),
        (PWAFunc{Concave,1}(planes), 0.9, -0.7875),
        # Single plane z = x and two planes z = |x|
        (PWAFunc{Convex,1}([Plane((1,), 0)]), 1.0, 1.0),
        (PWAFunc{Convex,1}([Plane((1,), 0), Plane((-1,), 0)]), 1.0, 1.0),
        # Three planes meeting at x = 0
        (PWAFunc{Convex,1}([Plane((1,), 0), Plane((-1,), 0), Plane((0,), 0), Plane((2,), -1)]), 1.5, 2.0),
    ]
    for (pwa, xlo, zopt) in cases, f in formulations
        status, xval, zval = _solve_pwa(pwa, f, xlo; bounds = [(-2, 2)])
        @test status == OPTIMAL
        @test xval ≈ xlo atol = 1e-6
        @test zval ≈ zopt atol = 1e-6
    end

    # Bounds are required for λ, Δ and Φ when x has no finite bounds
    pwa = PWAFunc{Convex,1}(planes)
    for f in [λ_Formulation, Δ_Formulation, Φ_Formulation]
        m = JuMP.Model()
        @variable(m, x ≥ 0)
        @test_throws ErrorException pwaffine(m, x, pwa; formulation = f())
    end

    # Δ and Φ are only available in 1D
    pwa2d = PWAFunc{Convex,2}([Plane([1, 0], 0), Plane([0, 1], 0), Plane([0, 0], 0)])
    for f in [Δ_Formulation, Φ_Formulation]
        m = JuMP.Model()
        @variable(m, x[1:2])
        @test_throws ErrorException pwaffine(m, x, pwa2d; formulation = f(), bounds = [(-1, 1), (-1, 1)])
    end

    # f(x,y) = max(x, y, 0) has a single vertex at the origin, but with bounds λ covers the box
    for f in [λ_Formulation, Ψ_Formulation]
        m = JuMP.Model(HiGHS.Optimizer)
        set_silent(m)
        @variable(m, -1 ≤ x[1:2] ≤ 1)
        z = pwaffine(m, x, pwa2d; formulation = f())
        @constraint(m, x[1] + x[2] ≥ 1.5)
        @objective(m, Min, z)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test objective_value(m) ≈ 0.75 atol = 1e-6
    end
end
