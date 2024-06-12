using PiecewiseAffineApprox
using GLMakie
using HiGHS

x = LinRange(0, 1, 20)
f(x) = first(x)^2
pwa = approx(
    f,
    [(0, 1)],
    Convex(),
    MILP(optimizer = HiGHS.Optimizer, planes = 3),
)
p = plot(x, f.(x), pwa)

using CairoMakie
save(joinpath(@__DIR__, "..", "docs", "approx.svg"), p; backend = CairoMakie)

function pwademo(x, f, Ns = 1:5; opt = HiGHS.Optimizer, C = Convex())
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1])

    scatter!(ax, x, f.(x), color = :red, markersize = 8)

    linerecords = []
    for n in Ns
        add_mplane!(ax, x, f, C, opt, n, linerecords)
        sleep(1)
    end

    return (; fig, ax, linerecords)
end

function add_mplane!(ax, x, f, C, opt, n, linerecords)
    x̄ = LinRange(minimum(x), maximum(x), 100)
    pwa = approx(f, [(0, 1)], C, MILP(optimizer = opt, planes = n))
    for ol in linerecords
        delete!(ax, ol)
    end
    empty!(linerecords)

    en = [-Inf for _ in x̄]
    # Plot each cutting plane
    for plane ∈ pwa.planes
        l = [evaluate(plane, i, C) for i ∈ x̄]
        nl = lines!(ax, x̄, l, linestyle = :dash)
        push!(linerecords, nl)
        en = max.(en, l)
    end
    # Plot effective envelope
    e = lines!(ax, x̄, en, color = :black)
    return push!(linerecords, e)
end

# With CairoMakie to save single frame as svg:

p = pwademo(x, f, 3:3)
save("approx.svg", p.fig; backend = CairoMakie)

# Create animation
record(p.fig, joinpath(@__DIR__, "..", "docs", "approxanim.mp4")) do io
    for n = 1:5
        add_mplane!(p.ax, x, f, Convex(), HiGHS.Optimizer, n, p.linerecords)
        for _ = 1:30
            recordframe!(io)
        end
    end
end

function rotated_plot()
    GLMakie.activate!()
    # 3D plot
    I = 100
    xmat = 2 * rand(2, I) .- 1
    x = [Tuple(xmat[:, i]) for i = 1:size(xmat, 2)]
    z = [p[1]^2 + p[2]^2 for p in x]
    vals = FunctionEvaluations(x, z)
    pwa = approx(
        vals,
        Convex(),
        Cluster(; optimizer = Xpress.Optimizer, planes = 9, strict = :none),
        # Cluster(; optimizer = HiGHS.Optimizer, planes = 9, strict = :none),
    )
    p = plot(vals, pwa)
    save(joinpath(@__DIR__, "..", "docs", "approx_3D.png"), p)

    # 3D animation
    x = xmat'
    xmin = minimum(x[:, 1])
    xmax = maximum(x[:, 1])

    ymin = minimum(x[:, 2])
    ymax = maximum(x[:, 2])
    x̄ = LinRange(xmin, xmax, 20)
    ȳ = LinRange(ymin, ymax, 20)
    fig = Figure(size = (1000, 1000))
    ax1 = Axis3(fig[1:2, 1:2])
    C = Convex()

    for p ∈ pwa.planes
        f = [
            evaluate(p, [x̄[i], ȳ[j]], C) for i ∈ eachindex(x̄),
            j ∈ eachindex(ȳ)
        ]
        surface!(ax1, x̄, ȳ, f)
    end
    scatter!(ax1, x[:, 1], x[:, 2], z, color = :red, markersize = 8)

    record(
        fig,
        joinpath(@__DIR__, "..", "docs", "rotation.mp4"),
        1:120,
    ) do frame
        return ax1.azimuth[] = 1.7pi + 0.3 * sin(2pi * frame / 120)
    end
end
rotated_plot()
