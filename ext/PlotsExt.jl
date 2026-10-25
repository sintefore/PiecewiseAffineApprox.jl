module PlotsExt

using PiecewiseAffineApprox
using Plots

function Plots.plot!(p::Plots.Plot, pwa::PWAFunc{C,1}, xlims) where {C}
    x̄ = LinRange(xlims[1], xlims[2], 100)

    for plane ∈ pwa.planes
        f = [evaluate(plane, x̄[i], C) for i ∈ eachindex(x̄)]
        Plots.plot!(p, x̄, f, legend = :none, linestyle = :dash)
    end
    f = [evaluate(pwa, x̄[i]) for i ∈ eachindex(x̄)]
    Plots.plot!(p, x̄, f, legend = :none, linecolor = :black)

    return p
end

function Plots.plot(pwa::PWAFunc{C,1}, xlims) where {C}
    return plot!(Plots.plot(), pwa, xlims)
end

# A difference of convex function, with its convex and concave parts in a
# second panel when parts = true
function Plots.plot(f::DCFunc{1}, xlims; parts = true)
    x̄ = LinRange(xlims[1], xlims[2], 200)
    v = [evaluate(f, (x̄[i],)) for i ∈ eachindex(x̄)]
    p₁ = Plots.plot(x̄, v, legend = :none, linecolor = :black)
    parts || return p₁
    p₂ = Plots.plot(title = "Convex and concave parts")
    for (part, color, label) ∈
        ((f.convex, :blue, "z₁ (convex)"), (f.concave, :orange, "z₂ (concave)"))
        isempty(part.planes) && continue
        for plane ∈ part.planes
            g = [
                evaluate(plane, x̄[i], typeof(part).parameters[1]) for
                i ∈ eachindex(x̄)
            ]
            Plots.plot!(
                p₂,
                x̄,
                g,
                linestyle = :dash,
                linecolor = color,
                linealpha = 0.5,
                label = "",
            )
        end
        g = [evaluate(part, (x̄[i],)) for i ∈ eachindex(x̄)]
        Plots.plot!(p₂, x̄, g, linecolor = color, label = label)
    end
    return Plots.plot(p₁, p₂, layout = (2, 1), size = (600, 700))
end

end
