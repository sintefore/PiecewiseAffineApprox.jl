# The test files are run in groups of about equal run time, each group in its
# own Julia process, with the groups in parallel. Much of the time of the tests
# is compilation, so the processes use -O0, which compiles faster.
#
# Arguments select the files whose names contain one of them, e.g.
# `Pkg.test(test_args = ["dc", "milp"])`. With the argument "--serial", the
# selected files run one after the other in this process instead.
using Test

const GROUPS = [
    ["test_dc.jl"],
    [
        "test_milp.jl",
        "test_vertex_enum.jl",
        "test_json.jl",
        "test_types.jl",
        "test_record.jl",
        "test_quadratic.jl",
    ],
    ["test_kazda_li.jl", "test_everywhere.jl", "test_cluster.jl", "test_linopt.jl"],
    ["test_aqua.jl", "test_interpolation.jl", "test_algorithms.jl"],
]
const FILES = reduce(vcat, GROUPS)

# Every test file must belong to a group, so that none is left out
let ungrouped = setdiff(filter(startswith("test_"), readdir(@__DIR__)), FILES)
    isempty(ungrouped) ||
        error("Test files not in a group in runtests.jl: $(join(ungrouped, ", "))")
end

function selected(args)
    names = filter(!=("--serial"), args)
    isempty(names) && return FILES
    files = filter(f -> any(a -> occursin(a, f), names), FILES)
    isempty(files) && error("No test files match $(join(names, ", "))")
    return files
end

if "--serial" ∈ ARGS || get(ENV, "PWA_TEST_WORKER", "") == "true"
    include("setup.jl")
    @testset "PiecewiseAffineApprox" begin
        with_logger(ConsoleLogger(devnull, Logging.Debug)) do
            for file ∈ selected(ARGS)
                include(file)
            end
        end
    end
else
    files = selected(ARGS)
    groups = filter(!isempty, [filter(∈(files), g) for g ∈ GROUPS])
    # The processes inherit the options of this one, e.g. for coverage and
    # deprecation warnings, and use the same environment
    function command(group)
        cmd = `$(Base.julia_cmd()) -O0 --project=$(Base.active_project())
            $(joinpath(@__DIR__, "runtests.jl")) $group`
        return addenv(cmd, "PWA_TEST_WORKER" => "true")
    end
    passed = fill(false, length(groups))
    @sync for (i, group) ∈ enumerate(groups)
        @async begin
            output = IOBuffer()
            t = @elapsed begin
                process = run(
                    pipeline(command(group); stdout = output, stderr = output);
                    wait = false,
                )
                wait(process)
            end
            passed[i] = success(process)
            status = passed[i] ? "passed" : "FAILED"
            println(
                "\n── Group $i ($status in $(round(t; digits = 1)) s): ",
                join(group, ", "),
            )
            print(String(take!(output)))
        end
    end
    @testset "PiecewiseAffineApprox" begin
        for (i, group) ∈ enumerate(groups)
            @testset "Group $i: $(join(group, ", "))" begin
                @test passed[i]
            end
        end
    end
end
