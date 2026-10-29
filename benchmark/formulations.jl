# Benchmark of the formulations of add_pwaformulation in a model with many piecewise
# affine functions: a dispatch problem with one cost function per unit,
# minimize the total cost subject to meeting a demand, with box bounds on
# the units.
#
# Usage (from the package root):
#
#     julia --project=benchmark benchmark/formulations.jl [--units=50,200,1000] [--repeats=2]
#
# Two families of cost functions are used: convex costs, which every
# formulation models as a linear program, and non-convex costs (a convex cost
# with a sinusoidal perturbation), which are difference of convex functions
# modelled with binary variables, by the Ψ-formulation with big-M selection
# and by the Δ- and Φ-formulations with binaries at the concave kinks only.
# Reports the build time, the solve time and the size of the model for each
# formulation, as the median over the repeats, to benchmark/results. The
# solves are limited to TIME_LIMIT seconds each; the status column shows
# which did not finish.

using PiecewiseAffineApprox
using HiGHS, JuMP, StableRNGs, Statistics, Printf, Dates

const TIME_LIMIT = 120.0
const OPTIMIZER = optimizer_with_attributes(
    HiGHS.Optimizer,
    MOI.Silent() => true,
    "time_limit" => TIME_LIMIT,
)

function parse_args(args)
    value(key, default) = begin
        a = findfirst(startswith("--$key="), args)
        isnothing(a) ? default : split(args[a], '='; limit = 2)[2]
    end
    units = parse.(Int, split(value("units", "50,200,1000"), ','))
    repeats = parse(Int, value("repeats", "2"))
    return (; units, repeats)
end

# The cost functions of the units: a quadratic-plus-reciprocal convex cost
# like the specific fuel consumption, with random coefficients, and the same
# with a sinusoidal perturbation for the non-convex case
function costs(rng, n, nonconvex)
    fs = map(1:n) do _
        a, b, c = 100 + 100rand(rng), 5 + 10rand(rng), 20 + 80rand(rng)
        w, h = 5 + 5rand(rng), nonconvex ? 10 + 10rand(rng) : 0.0
        return x -> a + b / x[1] + c * (x[1] - 0.7)^2 + h * sin(w * x[1])
    end
    return fs
end

function fit_costs(fs, nonconvex)
    return map(fs) do f
        if nonconvex
            approx(
                f,
                [(0.2, 1.0)],
                DC(),
                Interpol(; planes = 5, optimizer = OPTIMIZER);
                nsamples = 60,
            )
        else
            approx(
                f,
                [(0.2, 1.0)],
                Convex(),
                Interpol(; planes = 5);
                nsamples = 60,
            )
        end
    end
end

# Build and solve the dispatch model, returning the times and the model size
function dispatch(pwas, formulation, direction, demand)
    t₀ = time()
    m = Model(OPTIMIZER)
    n = length(pwas)
    @variable(m, 0.2 ≤ x[1:n] ≤ 1.0)
    z = [pwas[i](x[i]; formulation, direction) for i ∈ 1:n]
    @constraint(m, sum(x) ≥ demand)
    @objective(m, Min, sum(z))
    build = time() - t₀
    optimize!(m)
    solve = solve_time(m)
    status = termination_status(m)
    return (;
        build,
        solve,
        status,
        objective = primal_status(m) == MOI.FEASIBLE_POINT ? objective_value(m) : NaN,
        vars = num_variables(m),
        bins = count(is_binary, all_variables(m)),
        cons = sum(
            num_constraints(m, F, S) for (F, S) ∈ list_of_constraint_types(m)
        ),
    )
end

function main(args = ARGS)
    opts = parse_args(args)
    rows = []
    for nonconvex ∈ (false, true), n ∈ opts.units
        rng = StableRNG(7)
        pwas = fit_costs(costs(rng, n, nonconvex), nonconvex)
        demand = 0.6n
        formulations = if nonconvex
            [
                ("Psi-bigM", PsiFormulation()),
                ("Delta", DeltaFormulation()),
                ("Phi", PhiFormulation()),
            ]
        else
            [
                ("Psi", PsiFormulation()),
                ("Lambda", LambdaFormulation()),
                ("Delta", DeltaFormulation()),
                ("Phi", PhiFormulation()),
            ]
        end
        for (name, formulation) ∈ formulations
            direction = PiecewiseAffineApprox.Epigraph()
            # Compile first
            dispatch(pwas[1:min(n, 3)], formulation, direction, 0.6 * min(n, 3))
            runs = [
                dispatch(pwas, formulation, direction, demand) for
                _ ∈ 1:opts.repeats
            ]
            r = first(runs)
            row = (;
                family = nonconvex ? "nonconvex" : "convex",
                units = n,
                formulation = name,
                build = median(x.build for x ∈ runs),
                solve = median(x.solve for x ∈ runs),
                vars = r.vars,
                bins = r.bins,
                cons = r.cons,
                objective = r.objective,
                status = string(r.status),
            )
            @printf(
                "%-10s %5d %-9s build %7.3fs solve %8.3fs vars %6d bins %5d cons %6d obj %.4g %s\n",
                row.family,
                n,
                name,
                row.build,
                row.solve,
                row.vars,
                row.bins,
                row.cons,
                row.objective,
                row.status
            )
            push!(rows, row)
            flush(stdout)
        end
    end
    dir = joinpath(@__DIR__, "results")
    mkpath(dir)
    stamp = Dates.format(now(), "yyyymmdd-HHMMSS")
    md = joinpath(dir, "formulations-$stamp.md")
    open(md, "w") do io
        println(io, "# Formulations in a dispatch model\n")
        println(
            io,
            "- date: $(now())\n- julia: $VERSION\n- threads: $(Threads.nthreads())\n",
        )
        println(
            io,
            "Minimize the total cost of n units with piecewise affine cost functions of 5 planes (convex) or 5 + 5 planes (non-convex, difference of convex), subject to a demand. Times in seconds, median of $(opts.repeats) runs; solve time as reported by HiGHS.\n",
        )
        println(
            io,
            "| family | units | formulation | build | solve | variables | binaries | constraints | objective | status |",
        )
        println(io, "|---|---:|---|---:|---:|---:|---:|---:|---:|---|")
        for r ∈ rows
            @printf(
                io,
                "| %s | %d | %s | %.3f | %.3f | %d | %d | %d | %.6g | %s |\n",
                r.family,
                r.units,
                r.formulation,
                r.build,
                r.solve,
                r.vars,
                r.bins,
                r.cons,
                r.objective,
                r.status
            )
        end
    end
    println("\nResults written to\n  $md")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
