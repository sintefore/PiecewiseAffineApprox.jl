module MakieExt

using Makie
using PiecewiseAffineApprox

function Makie.plot(x, z, pwa::PWAFunc{C,2}) where {C}
    xmin = minimum(x[1, :])
    xmax = maximum(x[1, :])

    ymin = minimum(x[2, :])
    ymax = maximum(x[2, :])

    x̄ = LinRange(xmin, xmax, 20)
    ȳ = LinRange(ymin, ymax, 20)
    fig = Figure(size = (1000, 1000))
    ax1 = Axis3(fig[1:2, 1:2])
    ax2 = Axis(fig[1, 3])
    ax3 = Axis(fig[2, 3])

    for p ∈ pwa.planes
        f = [
            evaluate(p, [x̄[i], ȳ[j]], C) for i ∈ eachindex(x̄),
            j ∈ eachindex(ȳ)
        ]
        surface!(ax1, x̄, ȳ, f)
    end
    scatter!(ax1, x[1, :], x[2, :], z, color = :red, markersize = 8)

    l1 = PiecewiseAffineApprox._approx_error(x, z, pwa, :l1)
    l2 = PiecewiseAffineApprox._approx_error(x, z, pwa, :l2)
    lmax = PiecewiseAffineApprox._approx_error(x, z, pwa, :max)
    ax1.title = "l1 = $(round(l1, digits=2)), l2 = $(round(l2, digits=2)), max = $(round(lmax, digits=2))"
    z̄ = [abs(evaluate(pwa, x[:, i]) - z[i]) for i ∈ eachindex(z)]
    θ = 20 / maximum(z̄)
    ax2.title = "Max error = $(round(maximum(z̄),digits=2))"
    scatter!(ax2, x[1, :], x[2, :]; markersize = θ * z̄)

    𝒫 = PiecewiseAffineApprox._update_partition(x, pwa)
    for j ∈ eachindex(𝒫)
        x̄ = x[:, 𝒫[j]]
        scatter!(ax3, x̄[1, :], x̄[2, :], marker = :xcross)
    end

    return fig
end

function Makie.plot(input::FunctionEvaluations{2}, pwa::PWAFunc{C,2}) where {C}
    x = [p[i] for i ∈ 1:2, p ∈ input.points]
    z = input.values
    return Makie.plot(x, z, pwa)
end

function Makie.plot(pwa::PWAFunc{C,2}) where {C}
    verts = PiecewiseAffineApprox.vertices_by_subsets(pwa)
    x = [v[1][i] for i in 1:2, v in verts]
    z = [v[2] for v in verts]
    return Makie.plot(x, z, pwa)
end

function Makie.plot(x, y, pwa::PWAFunc{C,1}) where {C}
    active = getactive(C)
    inactive = getinactive(C)
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], limits = (nothing, ylims(y)))

    x̄ = LinRange(minimum(x), maximum(x), 100)

    scatter!(ax, x, y, color = :red, markersize = 8)
    en = [inactive(Inf, -Inf) for _ ∈ x̄]
    for plane ∈ pwa.planes
        f = [evaluate(plane, i, C) for i ∈ x̄]
        lines!(ax, x̄, f, linestyle = :dash)
        en = active.(en, f)
    end
    lines!(ax, x̄, en, color = :black)

    return fig
end

getactive(::Type{Convex}) = max
getactive(::Type{Concave}) = min
getinactive(::Type{Convex}) = min
getinactive(::Type{Concave}) = max
function ylims(y)
    L, U = extrema(y)
    diff = (U - L) * 0.1
    return (L - diff, U + diff)
end

end
