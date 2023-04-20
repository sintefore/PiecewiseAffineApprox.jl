module GLMakieExt

using GLMakie
using PiecewiseAffineApprox

function GLMakie.plot(x, z, pwl::PWLFunc{Convex,2})
    xmin = minimum(x[1, :])
    xmax = maximum(x[1, :])

    ymin = minimum(x[2, :])
    ymax = maximum(x[2, :])

    x̄ = LinRange(xmin, xmax, 20)
    ȳ = LinRange(ymin, ymax, 20)
    fig = Figure(resolution = (1000, 1000))
    ax1 = Axis3(fig[1:2, 1:2])
    ax2 = Axis(fig[1, 3])
    ax3 = Axis(fig[2, 3])

    scatter!(ax1, x[1, :], x[2, :], z, color = :red, markersize = 8)
    for p ∈ pwl.planes
        f = [evaluate(p, [x̄[i], ȳ[j]]) for i ∈ 1:length(x̄), j ∈ 1:length(ȳ)]
        surface!(ax1, x̄, ȳ, f)
    end
    l1 = PiecewiseAffineApprox._approx_error(x, z, pwl, :l1)
    l2 = PiecewiseAffineApprox._approx_error(x, z, pwl, :l2)
    lmax = PiecewiseAffineApprox._approx_error(x, z, pwl, :max)
    ax1.title = "l1 = $(round(l1, digits=2)), l2 = $(round(l2, digits=2)), max = $(round(lmax, digits=2))"
    z̄ = [abs(evaluate(pwl, x[:, i]) - z[i]) for i ∈ 1:length(z)]
    θ = 20 / maximum(z̄)
    ax2.title = "Max error = $(round(maximum(z̄),digits=2))"
    scatter!(ax2, x[1, :], x[2, :]; markersize = θ * z̄)

    𝒫 = PiecewiseAffineApprox._update_partition(x, pwl)
    for j ∈ 1:length(𝒫)
        x̄ = x[:, 𝒫[j]]
        scatter!(ax3, x̄[1, :], x̄[2, :], marker = :xcross)
    end

    return display(fig)
end

function GLMakie.plot(input::FunctionEvaluations{2}, pwl::PWLFunc{Convex,2})
    x = [p[i] for i ∈ 1:2, p ∈ input.points]
    z = input.values
    return GLMakie.plot(x, z, pwl)
end

end