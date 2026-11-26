using Aqua

@testset "Aqua" begin
    # Only ambiguities within this package, not those of JuMP and other deps.
    # The check for persistent tasks precompiles the package again in a new
    # process, which takes most of the time of these tests, and is left out.
    Aqua.test_all(
        PiecewiseAffineApprox;
        ambiguities = (recursive = false,),
        persistent_tasks = false,
    )
end
