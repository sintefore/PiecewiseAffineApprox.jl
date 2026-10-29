# Benchmark of the approximation methods on synthetic and real datasets.
#
# Usage (from the package root):
#
#     julia --project=benchmark -t auto benchmark/run.jl [options]
#
# Options:
#     --full            also run the larger scaling cases
#     --cases=REGEX     only cases whose name matches REGEX
#     --methods=REGEX   only methods whose name matches REGEX
#
# Writes a CSV file and a markdown report to benchmark/results.

using PiecewiseAffineApprox
using HiGHS, JuMP, StableRNGs, Statistics, Printf, Dates, Logging
import Pkg

include("datasets.jl")

const OPTIMIZER = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
const SEED = 2024

"""
A method to benchmark: `make(case)` returns the algorithm for a case, which is
only run if `applies(case)`. Methods with `random = true` are run twice with
the same seed to check that the result is reproducible.
"""
struct Method
    name::String
    applies::Function
    make::Function
    random::Bool
end

Method(name, applies, make) = Method(name, applies, make, false)

is1d(c) = dims(c) == 1
small(n) = c -> npoints(c) ≤ n

const METHODS = [
    # 1D interpolation with a fixed or selected number of planes
    Method("Interpol-l2", is1d,
        c -> Interpol(; planes = c.planes, metric = :l2, optimizer = OPTIMIZER)),
    Method("Interpol-l2-bic", is1d,
        c -> Interpol(; planes = 20, metric = :l2, selection = :bic,
            optimizer = OPTIMIZER)),
    Method("Interpol-l2-adjust", is1d,
        c -> Interpol(; planes = c.planes, metric = :l2,
            adjust = :breakpoints, optimizer = OPTIMIZER)),
    Method("Interpol-max", is1d,
        c -> Interpol(; planes = c.planes, metric = :max,
            optimizer = OPTIMIZER)),
    # Magnani & Boyd
    Method("Cluster-l2", c -> true,
        c -> Cluster(; planes = c.planes, rng = StableRNG(SEED)), true),
    Method("Cluster-l1", small(1000),
        c -> Cluster(; planes = c.planes, metric = :l1, optimizer = OPTIMIZER,
            rng = StableRNG(SEED)), true),
    # Kazda & Li
    # The band LP grows quickly with the number of points in several
    # dimensions (about 100 s for 1000 points in 2D)
    Method("Progressive-band", c -> is1d(c) || npoints(c) ≤ 1000,
        c -> Progressive(; tolerance = c.tolerance, strict = :none,
            optimizer = OPTIMIZER)),
    Method("Progressive-outer", c -> c.convexdata,
        c -> Progressive(; tolerance = c.tolerance, metric = :max,
            optimizer = OPTIMIZER)),
    Method("FullOrder", c -> c.convexdata && npoints(c) ≤ 2000,
        c -> FullOrder(; optimizer = OPTIMIZER)),
    # Mixed integer formulation, only for small instances
    Method("MILP-l1", small(100),
        c -> MILP(; planes = c.planes, metric = :l1, maxtime = 30,
            optimizer = OPTIMIZER)),
]

# Errors relative to the range of the function: root mean square and maximum
function errors(pwa, f::FunctionEvaluations, r)
    e = [evaluate(pwa, x) - v for (x, v) ∈ zip(f.points, f.values)]
    return sqrt(mean(abs2, e)) / r, maximum(abs, e) / r
end

function run_method(m::Method, c::Case)
    alg = m.make(c)
    GC.gc()
    stats = @timed approx(c.data, c.curvature, alg)
    pwa, time = stats.value, stats.time
    # Fast runs are repeated to reduce the timing noise
    for _ ∈ 1:2
        time < 1 || break
        time = min(time, @elapsed approx(c.data, c.curvature, m.make(c)))
    end
    r = valuerange(c)
    data_rms, data_max = errors(pwa, c.data, r)
    test_rms, test_max =
        isnothing(c.test) ? (NaN, NaN) : errors(pwa, c.test, r)
    reproducible = missing
    if m.random
        again = approx(c.data, c.curvature, m.make(c))
        reproducible = again == pwa
    end
    return (; planes = PiecewiseAffineApprox._planes(pwa), data_rms, data_max, test_rms,
        test_max, time, alloc_mib = stats.bytes / 2^20,
        reproducible, status = "ok")
end

const FAILED = (; planes = 0, data_rms = NaN, data_max = NaN, test_rms = NaN,
    test_max = NaN, time = NaN, alloc_mib = NaN, reproducible = missing)

function run_case(m::Method, c::Case)
    try
        # Warnings, e.g. about the selected number of planes, are expected and
        # would only clutter the output
        return with_logger(ConsoleLogger(stderr, Logging.Error)) do
            run_method(m, c)
        end
    catch e
        msg = first(split(sprint(showerror, e), '\n'))
        return (; FAILED..., status = "error: " * first(msg, 80))
    end
end

# Compile each method on a tiny case so that the timings exclude compilation
function warmup(methods)
    tiny = [synthetic("warmup-1d", x -> exp(x[1]), 1, 20; planes = 3),
        synthetic("warmup-2d", logsumexp, 2, 30; planes = 3)]
    for m ∈ methods, c ∈ tiny
        m.applies(c) && run_case(m, c)
    end
end

function parse_args(args)
    full = "--full" ∈ args
    pattern(key) = begin
        a = findfirst(startswith("--$key="), args)
        isnothing(a) ? r"" : Regex(split(args[a], '='; limit = 2)[2])
    end
    return (; full, cases = pattern("cases"), methods = pattern("methods"))
end

fmt(x::AbstractFloat) = isnan(x) ? "" : @sprintf("%.4g", x)
fmt(::Missing) = ""
fmt(x) = string(x)

const COLUMNS = (:planes, :data_rms, :data_max, :test_rms, :test_max, :time,
    :alloc_mib, :reproducible, :status)

function environment()
    commit = try
        readchomp(`git -C $(dirname(@__DIR__)) rev-parse --short HEAD`)
    catch
        "unknown"
    end
    deps = Dict(p.name => p.version for p ∈ values(Pkg.dependencies()))
    return [
        "date" => string(now()),
        "commit" => commit,
        "julia" => string(VERSION),
        "threads" => string(Threads.nthreads()),
        "cpu" => Sys.cpu_info()[1].model,
        "HiGHS.jl" => string(get(deps, "HiGHS", "?")),
        "JuMP" => string(get(deps, "JuMP", "?")),
    ]
end

function write_results(rows, env)
    dir = joinpath(@__DIR__, "results")
    mkpath(dir)
    stamp = Dates.format(now(), "yyyymmdd-HHMMSS")
    header = ["case", "group", "dims", "points", "truek", "method",
        string.(COLUMNS)...]
    csv = joinpath(dir, "results-$stamp.csv")
    open(csv, "w") do io
        println(io, join(header, ','))
        for r ∈ rows
            fields = [r.case, r.group, r.dims, r.points, something(r.truek, ""),
                r.method, (fmt(r.result[k]) for k ∈ COLUMNS)...]
            println(io, join(replace.(string.(fields), ',' => ';'), ','))
        end
    end
    md = joinpath(dir, "results-$stamp.md")
    open(md, "w") do io
        println(io, "# Benchmark results\n")
        for (k, v) ∈ env
            println(io, "- $k: $v")
        end
        println(io, "\nErrors are relative to the range of the function ",
            "(`data_*` against the data, `test_*` against the noise-free ",
            "truth at 2000 test points). Time in seconds, allocations in MiB.")
        for case ∈ unique(r.case for r ∈ rows)
            rs = filter(r -> r.case == case, rows)
            r1 = first(rs)
            k = isnothing(r1.truek) ? "" : ", true planes $(r1.truek)"
            println(io, "\n## $case\n\n$(r1.dims)D, $(r1.points) points$k\n")
            println(io, "| method | ", join(string.(COLUMNS), " | "), " |")
            println(io, "|---", "|---:"^(length(COLUMNS) - 1), "|---|")
            for r ∈ rs
                println(io, "| $(r.method) | ",
                    join((fmt(r.result[c]) for c ∈ COLUMNS), " | "), " |")
            end
        end
    end
    return csv, md
end

function main(args = ARGS)
    opts = parse_args(args)
    cases = [synthetic_cases(); nonconvex_cases(); real_cases()]
    opts.full && append!(cases, scaling_cases())
    filter!(c -> occursin(opts.cases, c.name), cases)
    methods = filter(m -> occursin(opts.methods, m.name), METHODS)
    env = environment()
    @info "Benchmark" cases = length(cases) methods = length(methods) env
    warmup(methods)
    rows = []
    for c ∈ cases, m ∈ methods
        m.applies(c) || continue
        result = run_case(m, c)
        @printf("%-28s %-20s planes %3d  test rms %-9s data rms %-9s %8.2fs  %s\n",
            c.name, m.name, result.planes, fmt(result.test_rms),
            fmt(result.data_rms), result.time, result.status)
        push!(rows, (; case = c.name, group = c.group, dims = dims(c),
            points = npoints(c), truek = c.truek, method = m.name, result))
    end
    csv, md = write_results(rows, env)
    println("\nResults written to\n  $csv\n  $md")
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
