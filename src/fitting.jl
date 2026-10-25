# The interface between the approximation algorithms and the rest of the
# package. An algorithm is a subtype of ApproximationAlgorithm that implements _fit for
# convex data; everything else (concave and difference of convex functions,
# sampling of functions, validation) is built on the traits below, which have
# defaults reading the fields of the same name.

"""
    ApproximationAlgorithm

The abstract type of the algorithms for computing approximations, passed to
[`approx`](@ref). An algorithm is a struct with its options as fields, usually
defined with `@kwdef`.

A new algorithm implements [`_fit`](@ref), which approximates convex data by
a convex piecewise affine function. Concave data is approximated by negation,
and non-convex data as a difference of convex functions, see [`DC`](@ref), by
the package. The options of the algorithm are read by the trait functions
[`planes`](@ref), [`concave_planes`](@ref), [`metric`](@ref), [`strict`](@ref)
and [`optimizer`](@ref), which by default return the fields of the same name,
and changed with [`with`](@ref), which by default uses the keyword constructor.
The functions [`supported_dimensions`](@ref) and
[`requires_convex_data`](@ref) describe what the algorithm accepts.
"""
ApproximationAlgorithm

"""
    _fit(a::ApproximationAlgorithm, data::FunctionEvaluations{D}) -> PWAFunc{Convex,D}

Approximate the `data` by a convex piecewise affine function. Implemented by
every algorithm. The data has at least one point and a number of variables
in [`supported_dimensions`](@ref) when called through [`approx`](@ref), and
the `metric` and `strict` options have been checked against
[`supported_metrics`](@ref) and [`supported_strict`](@ref). The algorithm
should validate its other options and throw an `ArgumentError` for
unsupported values.

The method may also return a tuple `(f, details)` with a named tuple of
diagnostics, which [`fit`](@ref) reports in the [`FitResult`](@ref).
"""
function _fit(a::ApproximationAlgorithm, ::FunctionEvaluations{D}) where {D}
    return throw(
        ArgumentError(
            "$(nameof(typeof(a))) does not implement _fit for functions of " *
            "$D variables, see the documentation of ApproximationAlgorithm",
        ),
    )
end

# The function and the diagnostics from the return value of _fit
_unwrap(f::AbstractPWAFunc) = (f, (;))
_unwrap(x::Tuple{AbstractPWAFunc,NamedTuple}) = x

"""
    FitResult

The result of [`fit`](@ref): the approximation `func`, its `error` from the
data in `metric` (see [`approximation_error`](@ref)), the `time` in seconds,
the `algorithm`, and `details`, a named tuple of diagnostics from the
algorithm, e.g. the errors of the trials of a [`Restart`](@ref), the scores
of a [`Select`](@ref) or the status of the [`MILP`](@ref) solve.
"""
struct FitResult{F<:AbstractPWAFunc,A<:ApproximationAlgorithm,N<:NamedTuple}
    func::F
    error::Float64
    metric::Symbol
    time::Float64
    algorithm::A
    details::N
end
function Base.show(io::IO, r::FitResult)
    print(io, "FitResult: ")
    summary(io, r.func)
    print(io, ", error = ", round(r.error; sigdigits = 4), " (", r.metric, ")")
    print(io, ", time = ", round(r.time; sigdigits = 3), " s")
    isempty(r.details) ||
        print(io, ", details: ", join(string.(keys(r.details)), ", "))
    return
end

# The metric of an algorithm for reporting errors: its own, or :l1
function _report_metric(a::ApproximationAlgorithm)
    m = metric(a)
    return isnothing(m) || m ∉ (:l1, :l2, :rms, :max) ? :l1 : m
end

# A FitResult from the fitted function and the diagnostics
function _result(data, a, func, details, time)
    m = _report_metric(a)
    return FitResult(
        func,
        approximation_error(func, data, m),
        m,
        time,
        a,
        details,
    )
end

# The field of the algorithm, or the default if it has no such field
function _option(a::ApproximationAlgorithm, name::Symbol, default)
    return hasfield(typeof(a), name) ? getfield(a, name) : default
end

"""
    planes(a::ApproximationAlgorithm)

The number of planes of the approximation computed by `a`, or `nothing` if
the algorithm determines it itself. Defaults to the field `planes`.
"""
planes(a::ApproximationAlgorithm) = _option(a, :planes, nothing)

"""
    concave_planes(a::ApproximationAlgorithm)

The number of planes of the concave part of a difference of convex
approximation computed by `a`, see [`DC`](@ref). Defaults to the field
`concave_planes`, and otherwise to [`planes`](@ref).
"""
function concave_planes(a::ApproximationAlgorithm)
    return _option(a, :concave_planes, planes(a))
end

"""
    metric(a::ApproximationAlgorithm)

The metric in which `a` measures the deviation from the data, `:l1`, `:l2`,
`:rms` or `:max`, or `nothing` if the algorithm has no such option. Defaults to
the field `metric`.
"""
metric(a::ApproximationAlgorithm) = _option(a, :metric, nothing)

"""
    strict(a::ApproximationAlgorithm)

Whether the approximation computed by `a` lies below the data (`:outer`),
above it (`:inner`) or neither (`:none`). Defaults to the field `strict`, and
to `:none` without one.
"""
strict(a::ApproximationAlgorithm) = _option(a, :strict, :none)

"""
    optimizer(a::ApproximationAlgorithm)

The optimizer used by `a`, or `nothing`. Defaults to the field `optimizer`.
"""
optimizer(a::ApproximationAlgorithm) = _option(a, :optimizer, nothing)

"""
    with(a::ApproximationAlgorithm; options...)

A copy of `a` with the given options changed, e.g. `with(a; planes = 3)`.
Defaults to calling the keyword constructor of the type of `a` with its
fields, so algorithms defined with `@kwdef` need no method. Throws an
`ArgumentError` for an option the algorithm does not have.
"""
function with(a::ApproximationAlgorithm; options...)
    T = typeof(a)
    for name ∈ keys(options)
        hasfield(T, name) ||
            throw(ArgumentError("$(nameof(T)) has no option $name"))
    end
    fields = (name => getfield(a, name) for name ∈ fieldnames(T))
    return T.name.wrapper(; fields..., options...)
end

"""
    supported_metrics(a::ApproximationAlgorithm)

The values of the `metric` option that `a` accepts, checked by [`approx`](@ref)
before the algorithm is called. Defaults to `(:l1, :l2, :rms, :max)`.
"""
supported_metrics(::ApproximationAlgorithm) = (:l1, :l2, :rms, :max)

"""
    supported_strict(a::ApproximationAlgorithm)

The values of the `strict` option that `a` accepts, checked by [`approx`](@ref)
before the algorithm is called. Defaults to `(:none, :inner, :outer)`.
"""
supported_strict(::ApproximationAlgorithm) = (:none, :inner, :outer)

# Check the metric and strict options of an algorithm, if it has them
function _validate_options(a::ApproximationAlgorithm)
    m = metric(a)
    isnothing(m) || _check_option(a, "metric", m, supported_metrics(a))
    _check_option(a, "strict", strict(a), supported_strict(a))
    return a
end

"""
    supported_dimensions(a::ApproximationAlgorithm)

The numbers of variables of the functions `a` can approximate. Defaults to
all, [`Interpol`](@ref) returns `1:1`.
"""
supported_dimensions(::ApproximationAlgorithm) = 1:typemax(Int)

"""
    requires_convex_data(a::ApproximationAlgorithm)

Whether `a` requires data that can be interpolated by a convex function, as
[`Progressive`](@ref) and [`FullOrder`](@ref) do. Such algorithms are not used
to refit the parts of a difference of convex approximation. Defaults to
`false`.
"""
requires_convex_data(::ApproximationAlgorithm) = false

# Throw an ArgumentError unless the option has one of the allowed values
function _check_option(a::ApproximationAlgorithm, name, value, allowed)
    value ∈ allowed && return value
    return throw(
        ArgumentError(
            "Unsupported $name = :$value for $(nameof(typeof(a))), use one of " *
            join(":" .* string.(allowed), ", "),
        ),
    )
end

"""
    approximation_error(f::AbstractPWAFunc, data::FunctionEvaluations, metric = :l1)
    approximation_error(f::AbstractPWAFunc, X::AbstractMatrix, z, metric = :l1)

The deviation of the approximation `f` from the `data` per point: the mean
absolute deviation (`:l1`), the root mean square deviation (`:l2` or `:rms`)
or the maximum absolute deviation (`:max`). The points can also be given as
the columns of the matrix `X` with the values `z`.
"""
function approximation_error(
    f::AbstractPWAFunc,
    data::FunctionEvaluations,
    metric = :l1,
)
    return _deviation_measure(
        [abs(evaluate(f, x) - z) for (x, z) ∈ data],
        metric,
    )
end
function approximation_error(
    f::AbstractPWAFunc,
    X::AbstractMatrix,
    z::AbstractVector,
    metric = :l1,
)
    return _deviation_measure(
        [abs(evaluate(f, x) - zᵢ) for (x, zᵢ) ∈ zip(eachcol(X), z)],
        metric,
    )
end
function _deviation_measure(r, metric)
    if metric == :l1
        return mean(r)
    elseif metric == :l2 || metric == :rms
        return sqrt(mean(abs2, r))
    elseif metric == :max
        return maximum(r)
    end
    return throw(
        ArgumentError(
            "Unsupported metric :$metric for the approximation error, use " *
            ":l1, :l2 or :max",
        ),
    )
end

"""
    Select(inner::ApproximationAlgorithm; criterion = :aic, tolerance = 0.0, folds = 5, metric)

Select the number of planes for the algorithm `inner`, which is fitted with
`1` to `planes(inner)` planes, by `criterion`:

- `:tolerance`: the fewest planes with an [`approximation_error`](@ref) of at
  most `tolerance` in `metric`.
- `:aic`, `:bic`: minimize the Akaike or Bayesian information criterion,
  assuming normally (`:l2`) or Laplace (`:l1`) distributed errors, with
  `D + 1` parameters per plane for functions of `D` variables.
- `:cv`: `folds`-fold cross validation, with the distinct points assigned to
  folds by interleaving them in sorted order (repeated measurements at a
  point are held out together), and the one standard error rule: the fewest
  planes with an error on the held-out points within one standard error of
  the minimum.

The `metric` defaults to that of `inner`, or `:l1`. The `selection` option of
[`Interpol`](@ref) is a shorthand for this wrapper. `Select` works with any
algorithm, e.g. `Select(Cluster(planes = 8))` or
`Select(MILP(planes = 6; optimizer); criterion = :tolerance, tolerance = 0.01)`.
A warning is given when the maximum number of planes is selected.
"""
@kwdef struct Select{A<:ApproximationAlgorithm} <: ApproximationAlgorithm
    inner::A
    criterion::Symbol = :aic
    tolerance::Float64 = 0.0
    folds::Int = 5
    metric::Symbol = something(PiecewiseAffineApprox.metric(inner), :l1)
end
Select(inner::ApproximationAlgorithm; kwargs...) = Select(; inner, kwargs...)

planes(s::Select) = planes(s.inner)
concave_planes(s::Select) = concave_planes(s.inner)
strict(s::Select) = strict(s.inner)
optimizer(s::Select) = optimizer(s.inner)
supported_dimensions(s::Select) = supported_dimensions(s.inner)
requires_convex_data(s::Select) = requires_convex_data(s.inner)
function with(s::Select; options...)
    for name ∈ keys(options)
        hasfield(typeof(s.inner), name) ||
            hasfield(Select, name) ||
            throw(
                ArgumentError(
                    "Neither Select nor $(nameof(typeof(s.inner))) has an option $name",
                ),
            )
    end
    inner = (k => v for (k, v) ∈ options if hasfield(typeof(s.inner), k))
    outer = (k => v for (k, v) ∈ options if !hasfield(typeof(s.inner), k))
    s =
        isempty(inner) ? s :
        Select(; _fields(s)..., inner = with(s.inner; inner...))
    return isempty(outer) ? s : Select(; _fields(s)..., outer...)
end
_fields(a) = (name => getfield(a, name) for name ∈ fieldnames(typeof(a)))

function Base.show(io::IO, s::Select)
    print(io, "Select(", s.inner, "; criterion = :", s.criterion)
    s.criterion == :tolerance && print(io, ", tolerance = ", s.tolerance)
    s.criterion == :cv && print(io, ", folds = ", s.folds)
    return print(io, ", metric = :", s.metric, ")")
end

function _fit(s::Select, data::FunctionEvaluations{D}) where {D}
    _check_option(s, "criterion", s.criterion, (:tolerance, :aic, :bic, :cv))
    K = planes(s.inner)
    isnothing(K) && throw(
        ArgumentError(
            "Select requires an algorithm with a planes option, got " *
            "$(nameof(typeof(s.inner)))",
        ),
    )
    K ≥ 1 || throw(ArgumentError("Select requires planes ≥ 1"))
    s.criterion == :tolerance &&
        s.tolerance ≤ 0 &&
        throw(ArgumentError("criterion = :tolerance requires tolerance > 0"))
    s.criterion == :cv &&
        s.folds < 2 &&
        throw(ArgumentError("criterion = :cv requires folds ≥ 2"))
    s.criterion ∈ (:aic, :bic) &&
        s.metric == :max &&
        throw(
            ArgumentError(
                "criterion = :$(s.criterion) requires metric :l1 or :l2",
            ),
        )
    N = length(data)
    fitp(p, d) = approx(d, Convex(), with(s.inner; planes = p))

    # The score of each number of planes, reported in the details
    scores = Float64[]
    k, fits = if s.criterion == :cv
        E = _cv_errors(s, data, K, fitp)
        F = size(E, 1)
        scores = vec(mean(E, dims = 1))
        se = vec(std(E, dims = 1)) ./ sqrt(F)
        kᵐⁱⁿ = argmin(scores)
        findfirst(≤(scores[kᵐⁱⁿ] + se[kᵐⁱⁿ]), scores), nothing
    else
        fits = [fitp(p, data) for p ∈ 1:K]
        e = [approximation_error(f, data, s.metric) for f ∈ fits]
        if s.criterion == :tolerance
            scores = e
            k = findfirst(≤(s.tolerance), e)
            if isnothing(k)
                @warn "Tolerance $(s.tolerance) not reached with $K planes " *
                      "(error $(minimum(e)))"
                k = K
            end
            k, fits
        else
            # -2 log-likelihood for Gaussian (:l2) or Laplace (:l1) errors
            # from the mean deviation per point, with D + 1 parameters per
            # plane. Errors at the level of rounding errors are treated as
            # exact fits, so that these are compared equal.
            scale = max(maximum(abs, data.values), floatmin())
            e = max.(e, 1e-12 * scale)
            m2ll =
                s.metric == :max ? e :
                s.metric == :l1 ? 2N .* log.(2 .* e) .+ 2N :
                N .* log.(2π .* e .^ 2) .+ N
            nparams = [(D + 1) * _planes(f) for f ∈ fits]
            penalty = s.criterion == :aic ? 2.0 : log(N)
            scores = m2ll .+ penalty .* nparams
            argmin(scores), fits
        end
    end
    if s.criterion != :tolerance && k == K
        @warn "The maximum number of planes ($K) was selected by " *
              ":$(s.criterion), consider increasing planes"
    end
    @debug "Select chose $k planes by :$(s.criterion)"
    func = isnothing(fits) ? fitp(k, data) : fits[k]
    return func, (criterion = s.criterion, planes = k, scores = scores)
end

# Cross validation errors E[f, p] for fold f and p = 1, ..., K planes, with
# the distinct points assigned to the folds in an interleaved manner (sorted
# along the first coordinate, so that every fold covers the range of the
# data). Repeated measurements at a point are held out together.
function _cv_errors(s::Select, data::FunctionEvaluations, K, fit)
    N = length(data)
    distinct = sort!(unique(data.points))
    F = min(s.folds, length(distinct))
    F ≥ 2 || throw(ArgumentError("criterion = :cv requires folds ≥ 2"))
    group = Dict(p => mod(i - 1, F) + 1 for (i, p) ∈ enumerate(distinct))
    fold = [group[p] for p ∈ data.points]
    E = fill(Inf, F, K)
    for f ∈ 1:F
        train = findall(≠(f), fold)
        test = findall(==(f), fold)
        length(train) ≥ 2 ||
            throw(ArgumentError("Too few points for $F-fold cross validation"))
        dtrain = FunctionEvaluations(data.points[train], data.values[train])
        dtest = FunctionEvaluations(data.points[test], data.values[test])
        for p ∈ 1:K
            E[f, p] = approximation_error(fit(p, dtrain), dtest, s.metric)
        end
    end
    return E
end

"""
    Restart(inner::ApproximationAlgorithm; trials = 20, rng = Random.default_rng(), metric)

Run the randomized algorithm `inner` `trials` times from different random
starts, and return the best approximation as measured by `metric` (default:
that of `inner`, or `:l1`). The algorithm must have an `rng` option, which is
set to a seeded generator for each trial, with the seeds drawn from `rng` up
front, so that the result is reproducible and independent of the number of
threads. If julia is started with multiple threads, the trials run in
parallel; consider how many threads are beneficial, particularly with a
commercial solver whose license may restrict the number of simultaneous
instances. Errors in a trial are reported directly.

The `trials` option of [`Cluster`](@ref) is a shorthand for this wrapper.
"""
@kwdef struct Restart{A<:ApproximationAlgorithm,R<:AbstractRNG} <:
              ApproximationAlgorithm
    inner::A
    trials::Int = 20
    rng::R = Random.default_rng()
    metric::Symbol = something(PiecewiseAffineApprox.metric(inner), :l1)
end
Restart(inner::ApproximationAlgorithm; kwargs...) = Restart(; inner, kwargs...)

planes(r::Restart) = planes(r.inner)
concave_planes(r::Restart) = concave_planes(r.inner)
strict(r::Restart) = strict(r.inner)
optimizer(r::Restart) = optimizer(r.inner)
supported_dimensions(r::Restart) = supported_dimensions(r.inner)
requires_convex_data(r::Restart) = requires_convex_data(r.inner)
function with(r::Restart; options...)
    for name ∈ keys(options)
        hasfield(typeof(r.inner), name) ||
            hasfield(Restart, name) ||
            throw(
                ArgumentError(
                    "Neither Restart nor $(nameof(typeof(r.inner))) has an option $name",
                ),
            )
    end
    # Options of both go to the wrapper, as its rng and metric are the ones
    # that matter for the restarts
    inner = (
        k => v for (k, v) ∈ options if
        hasfield(typeof(r.inner), k) && !hasfield(Restart, k)
    )
    outer = (k => v for (k, v) ∈ options if hasfield(Restart, k))
    r =
        isempty(inner) ? r :
        Restart(; _fields(r)..., inner = with(r.inner; inner...))
    return isempty(outer) ? r : Restart(; _fields(r)..., outer...)
end

function Base.show(io::IO, r::Restart)
    return print(
        io,
        "Restart(",
        r.inner,
        "; trials = ",
        r.trials,
        ", metric = :",
        r.metric,
        ")",
    )
end

# Whether the algorithm has an rng option, through wrappers
_has_rng(a::ApproximationAlgorithm) = hasfield(typeof(a), :rng)
_has_rng(s::Select) = _has_rng(s.inner)

function _fit(r::Restart, data::FunctionEvaluations)
    r.trials ≥ 1 || throw(ArgumentError("Restart requires trials ≥ 1"))
    _has_rng(r.inner) || throw(
        ArgumentError(
            "Restart requires an algorithm with an rng option, got " *
            "$(nameof(typeof(r.inner)))",
        ),
    )
    # One seed per trial up front, so that the result does not depend on how
    # the trials are scheduled across threads
    seeds = rand(r.rng, UInt64, r.trials)
    tasks = [
        @spawn _fit(with(r.inner; rng = StableRNG(seeds[i])), data) for
        i ∈ 1:r.trials
    ]
    fits = map(tasks) do task
        try
            _unwrap(fetch(task))
        catch e
            # Report the error of the trial rather than the failure of the task
            e isa TaskFailedException ? throw(e.task.result) : rethrow()
        end
    end
    errors = [approximation_error(f, data, r.metric) for (f, _) ∈ fits]
    best = argmin(errors)
    @debug "Restart: best of $(r.trials) trials has error $(errors[best]) ($(r.metric))"
    func, inner = fits[best]
    return func,
    (trials = r.trials, errors = errors, best = best, inner = inner)
end
