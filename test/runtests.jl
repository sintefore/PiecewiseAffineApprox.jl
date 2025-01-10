using COSMO
using HiGHS
using JuMP
using JSON3
using Logging
using Pajarito
using PiecewiseAffineApprox
using StableRNGs
using Test
rng = StableRNG(123)

const PWA = PiecewiseAffineApprox
const optimizer =
    optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
const qp_optimizer = optimizer_with_attributes(
    Pajarito.Optimizer,
    "oa_solver" => optimizer_with_attributes(
        HiGHS.Optimizer,
        MOI.Silent() => true,
        "mip_feasibility_tolerance" => 1e-8,
        "mip_rel_gap" => 1e-6,
    ),
    "conic_solver" =>
        optimizer_with_attributes(COSMO.Optimizer, MOI.Silent() => true),
    MOI.Silent() => true,
)

@testset "PiecewiseAffineApprox" begin
    nologger = ConsoleLogger(devnull, Logging.Debug)
    with_logger(nologger) do
        include("testquadapprox.jl")
        include("convex_interpolation.jl")
        include("testconvexapprox.jl")
        include("test_2D_convexapprox.jl")
        include("test_big_M.jl")
        include("test_structtype.jl")
        return include("test_kazda_li.jl")
    end
end
