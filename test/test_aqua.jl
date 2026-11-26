using Aqua

@testset "Aqua" begin
    # Only ambiguities within this package, not those of JuMP and other deps
    Aqua.test_all(PiecewiseAffineApprox; ambiguities = (recursive = false,))
end
