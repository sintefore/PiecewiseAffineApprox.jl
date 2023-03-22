function plot!(p::Plots.Plot, pwl::PWLFunc{C,1}, xlims) where {C}
    x̄ = LinRange(xlims[1], xlims[2], 100)

    for plane ∈ pwl.planes
        f = [evaluate(plane, x̄[i]) for i ∈ eachindex(x̄)]
        Plots.plot!(p, x̄, f, legend = :none, linestyle = :dash)
    end
    f = [evaluate(pwl, x̄[i]) for i ∈ eachindex(x̄)]
    Plots.plot!(p, x̄, f, legend = :none, linecolor = :black)

    return p
end

function plot(pwl::PWLFunc{C,1}, xlims) where {C}
    return plot!(Plots.plot(), pwl, xlims)
end
