module PlotsExt

using PiecewiseAffineApprox
using Plots

function Plots.plot!(p::Plots.Plot, pwa::PWAFunc{C,1}, xlims) where {C}
    x̄ = LinRange(xlims[1], xlims[2], 100)

    for plane ∈ pwa.planes
        f = [evaluate(plane, x̄[i]) for i ∈ eachindex(x̄)]
        Plots.plot!(p, x̄, f, legend = :none, linestyle = :dash)
    end
    f = [evaluate(pwa, x̄[i]) for i ∈ eachindex(x̄)]
    Plots.plot!(p, x̄, f, legend = :none, linecolor = :black)

    return p
end

function Plots.plot(pwa::PWAFunc{C,1}, xlims) where {C}
    return plot!(Plots.plot(), pwa, xlims)
end

end
