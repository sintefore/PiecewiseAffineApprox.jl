using JuMP
using HiGHS
using Juniper
using Ipopt
using Pajarito
using Hypatia
using Random

using PiecewiseAffineApprox
PWA = PiecewiseAffineApprox
optimizer = HiGHS.Optimizer

f(x) = x[1]^2 + x[2]^2
vals = PWA._sample_uniform(f, [(-1, 1), (-1, 1)], 10)
pwa = approx(
    vals,
    Convex(),
    Cluster(optimizer = optimizer, planes = 10, metric = :l2),
)

struct Instance
    size::Any
    rnd::Any
    pwa::Any
end

function Instance(size, pwa)
    Random.seed!(42)
    rnd = rand(Float64, size)

    return Instance(size, rnd, pwa)
end

function linear_model(instance)
    model = Model()
    I = 1:instance.size
    @variable(model, x[I], Int)
    @variable(model, y[I])
    z = []
    for i ∈ I
        zi = PWA.pwaffine(model, [x[i], y[i]], instance.pwa)
        @constraint(model, zi >= x[i] + y[i] + 2 + instance.rnd[i])
        push!(z, zi)
    end
    @constraint(model, sum(x) + 5 <= 0.5 * sum(y))
    @objective(model, Min, sum(z))

    return model
end

function nonlin_model(instance)
    model_nl = Model()
    I = 1:instance.size
    @variable(model_nl, x[I], Int)
    @variable(model_nl, y[I])
    @variable(model_nl, z[I])
    @constraint(model_nl, [i in I], z[i] >= x[i]^2 + y[i]^2)
    @constraint(model_nl, [i in I], z[i] >= x[i] + y[i] + 2 + instance.rnd[i])
    @constraint(model_nl, sum(x) + 5 <= 0.5 * sum(y))
    @objective(model_nl, Min, sum(z))

    return model_nl
end

# Setup solvers 
highs = HiGHS.Optimizer

ipopt = optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
juniper = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)
oa_solver = optimizer_with_attributes(
    HiGHS.Optimizer,
    MOI.Silent() => true,
    "mip_feasibility_tolerance" => 1e-8,
    "mip_rel_gap" => 1e-6,
)
hypatia = optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => true)
pajarito = optimizer_with_attributes(
    Pajarito.Optimizer,
    "time_limit" => 60,
    "oa_solver" => oa_solver,
    "conic_solver" => hypatia,
)

results = NamedTuple[]

for sz ∈ StepRange(5, 5, 20)
    instance = Instance(sz, pwa)

    model_lin = linear_model(instance)
    model_nl = nonlin_model(instance)

    set_optimizer(model_lin, highs)
    time_opt = @timed optimize!(model_lin)
    obj_val = objective_value(model_lin)
    push!(
        results,
        (
            size = sz,
            formulation = "pwa",
            solver = "highs",
            time = time_opt.time,
            obj_val = obj_val,
        ),
    )

    #set_optimizer(model_nl, juniper)
    #optimize!(model_nl)

    set_optimizer(model_nl, pajarito)
    time_opt = @timed optimize!(model_nl)
    obj_val = objective_value(model_nl)
    push!(
        results,
        (
            size = sz,
            formulation = "nonlin",
            solver = "highs",
            time = time_opt.time,
            obj_val = obj_val,
        ),
    )
end
