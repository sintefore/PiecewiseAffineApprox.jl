@testset "Linear optimization 1D" begin

    # 8 plane approximation of z = (x-1)^2
    pwa = JSON.parsefile(joinpath(@__DIR__, "pwa_ex", "pwa_quad_1d_8planes.json"), PWAFunc)
    v = PWA.vertices_by_subsets(pwa)
    x_min, z_min = argmin(v -> v[2], v)
    # Test all 4 possible formulations in 1D
    for f in [λ_Formulation, Ψ_Formulation, Δ_Formulation, Φ_Formulation]
        m = JuMP.Model(HiGHS.Optimizer)
        set_silent(m)
        @variable(m, x)
        @variable(m, z)
        @objective(m, Min, z)
        pwaffine(m, x, pwa; z = z, formulation = f())
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value(x) ≈ x_min[1]
        @test objective_value(m) ≈ z_min

        # Restrict x-values
        @constraint(m, x >= 1.5)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value(x) ≈ 1.5
        @test objective_value(m) ≈ evaluate(pwa, 1.5)
    end
end

@testset "Linear optimization 2D" begin
    # Test 1: Minimize x₁² + x₂² subject to x₁ + 3x₂ ≥ 2
    # Analytical solution: minimize ||x||² with x₁ + 3x₂ = 2
    # Solution: x = (2/10, 6/10) = (0.2, 0.6) with objective = 0.4
    pwa = JSON.parsefile(joinpath(@__DIR__, "pwa_ex", "pwa_quad_2d_12planes.json"), PWAFunc)
    v = PWA.vertices_by_subsets(pwa)
    x_min, z_min = argmin(v -> v[2], v)

    # Test both possible formulations in 2D
    for f in [λ_Formulation, Ψ_Formulation]
        m = JuMP.Model(HiGHS.Optimizer)
        set_silent(m)
        @variable(m, x[1:2])
        @variable(m, z)
        pwaffine(m, x, pwa; z = z, formulation = f())
        @objective(m, Min, z)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value.(x) ≈ x_min
        @test objective_value(m) ≈ z_min

        # Restrict x-values
        @constraint(m, x[1] + 3*x[2] ≥ 2)
        optimize!(m)
        @test termination_status(m) == OPTIMAL
        @test value.(x) ≈ [0.29, 0.57]
        @test objective_value(m) ≈ 0.33739 atol=1e-5
        @test value(z) ≈ evaluate(pwa, value.(x)) atol=1e-5
    end
end
