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

    l1 = approximation_error(pwa, x, z, :l1)
    l2 = approximation_error(pwa, x, z, :l2)
    lmax = approximation_error(pwa, x, z, :max)
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
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], limits = (nothing, ylims(y)))

    x̄ = LinRange(minimum(x), maximum(x), 100)

    scatter!(ax, x, y, color = :red, markersize = 8)
    _plot_planes!(ax, x̄, pwa)

    return fig
end

# Plot the planes of a one-dimensional pwa function (dashed) and the function
# itself (solid) on the axis, returning the plot of the function
function _plot_planes!(
    ax,
    x̄,
    pwa::PWAFunc{C,1};
    color = nothing,
    planes = true,
    label = nothing,
) where {C}
    active = getactive(C)
    inactive = getinactive(C)
    en = [inactive(Inf, -Inf) for _ ∈ x̄]
    for plane ∈ pwa.planes
        f = [evaluate(plane, i, C) for i ∈ x̄]
        if planes && isnothing(color)
            lines!(ax, x̄, f, linestyle = :dash)
        elseif planes
            lines!(ax, x̄, f, linestyle = :dash, color = (color, 0.5))
        end
        en = active.(en, f)
    end
    return lines!(ax, x̄, en; color = something(color, :black), label)
end

"""
    plot(x, y, f::DCFunc{1}; parts = true)
    plot(x, z, f::DCFunc{2}; parts = true)
    plot(input::FunctionEvaluations, f::DCFunc; parts = true)

Plot a difference of convex approximation together with the data points. With
`parts = true` (default), the convex part z₁ and the concave part z₂ of the
function z = z₁ + z₂ are shown in separate axes, with their planes.
"""
function Makie.plot(x, y, f::DCFunc{1}; parts = true)
    fig = Figure(size = (600, parts ? 700 : 400))
    ax = Axis(fig[1, 1], limits = (nothing, ylims(y)))
    x̄ = LinRange(minimum(x), maximum(x), 200)

    scatter!(ax, x, y, color = :red, markersize = 8)
    lines!(ax, x̄, [evaluate(f, (i,)) for i ∈ x̄], color = :black)

    if parts
        # The planes of the parts extend far beyond the parts themselves, so
        # the axis is limited to the range of the parts
        v = [
            PiecewiseAffineApprox._value_or_zero(part, (i,)) for
            part ∈ (f.convex, f.concave), i ∈ x̄
        ]
        ax2 = Axis(
            fig[2, 1],
            title = "Convex and concave parts",
            limits = (nothing, ylims(vec(v))),
        )
        _plot_planes!(ax2, x̄, f.convex; color = :blue, label = "z₁ (convex)")
        _plot_planes!(
            ax2,
            x̄,
            f.concave;
            color = :orange,
            label = "z₂ (concave)",
        )
        axislegend(ax2, position = :ct)
    end

    return fig
end

function Makie.plot(x, z, f::DCFunc{2}; parts = true)
    xmin, xmax = extrema(x[1, :])
    ymin, ymax = extrema(x[2, :])
    x̄ = LinRange(xmin, xmax, 40)
    ȳ = LinRange(ymin, ymax, 40)
    fig = Figure(size = (1000, 1000))
    ax1 = Axis3(fig[1:(parts ? 3 : 2), 1:2])
    ax2 = Axis(fig[1, 3], title = "Error at the data points")

    v = [evaluate(f, (x̄[i], ȳ[j])) for i ∈ eachindex(x̄), j ∈ eachindex(ȳ)]
    surface!(ax1, x̄, ȳ, v)
    scatter!(ax1, x[1, :], x[2, :], z, color = :red, markersize = 8)

    l1 = approximation_error(f, x, z, :l1)
    l2 = approximation_error(f, x, z, :l2)
    lmax = approximation_error(f, x, z, :max)
    ax1.title = "l1 = $(round(l1, digits=2)), l2 = $(round(l2, digits=2)), max = $(round(lmax, digits=2))"
    z̄ = [abs(evaluate(f, x[:, i]) - z[i]) for i ∈ eachindex(z)]
    θ = 20 / max(maximum(z̄), eps())
    ax2.title = "Max error = $(round(maximum(z̄),digits=2))"
    scatter!(ax2, x[1, :], x[2, :]; markersize = θ * z̄)

    if parts
        x̄ = LinRange(xmin, xmax, 20)
        ȳ = LinRange(ymin, ymax, 20)
        for (k, (part, C, name)) ∈ enumerate((
            (f.convex, Convex, "z₁ (convex)"),
            (f.concave, Concave, "z₂ (concave)"),
        ))
            ax = Axis3(fig[1+k, 3], title = name)
            for p ∈ part.planes
                v = [
                    evaluate(p, [x̄[i], ȳ[j]], C) for i ∈ eachindex(x̄),
                    j ∈ eachindex(ȳ)
                ]
                surface!(ax, x̄, ȳ, v)
            end
        end
    end

    return fig
end

function Makie.plot(input::FunctionEvaluations{2}, f::DCFunc{2}; kwargs...)
    x = [p[i] for i ∈ 1:2, p ∈ input.points]
    return Makie.plot(x, input.values, f; kwargs...)
end
function Makie.plot(input::FunctionEvaluations{1}, f::DCFunc{1}; kwargs...)
    return Makie.plot(first.(input.points), input.values, f; kwargs...)
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
