# PiecewiseAffineApprox

Add convex (or concave) piecewise linear approximations of functions or a set of points to optimization models modelled in [JuMP](https://jump.dev/). 

This package provides two main methods to fit a set of points: 

1. creates and solves a MILP to fit a set of points, and adds the resulting linear constraints to the optimization model. This method is partially based on [Toriello & Vielma, 2012](https://doi.org/10.1016/j.ejor.2011.12.030). 
2. uses a herusitic to fit the set of points. This method is based on [Magnani & Boyd, 2009](https://doi.org/10.1007/s11081-008-9045-3).

For non-convex functions, consider using [PiecewiseLinearOpt.jl](https://github.com/joehuchette/PiecewiseLinearOpt.jl).

## Usage

```julia
using JuMP, PiecewiseAffineApprox, HiGHS

m = Model(HiGHS.Optimizer)
@variable(m, x)
# Create a piecewise linear approximation to x^2 on the interval [-1, 1]
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), Optimized(optimizer = HiGHS.Optimizer, planes=5))
# Add the pwa function to the model
z = pwaffine(m, x, pwa)
# Minimize
@objective(m, Min, z)
# Check approximation/solution at x = 0.5
@constraint(m, x >= 0.5)
optimize!(m)
value(z) # 0.2653
```

To keep dependencies light, PiecewiseAffineApprox does not include plotting by default. If the Plots package is loaded
before using the module, some simple plotting routines will be available

The following demonstrates how this can be achieved:

```julia
using Plots, PiecewiseAffineApprox, HiGHS

function plotquademo(N = 3, opt = HiGHS.Optimizer)
    pwa = approx(x -> x[1]^2, [(0,1)], Convex(), Optimized( optimizer = opt; planes=N))
    x = LinRange(0, 1, 20)
    p = plot(x, x.^2, seriestype=:scatter, markershape=:x, ylims=(-0.5,1))
    PiecewiseAffineApprox.plot!(p, pwa, (0,1))
    return p
end

p = plotquademo()
savefig(p, "approx.svg")
```
![](docs/approx.svg)

Animation showing the accuracy when adding more cuts:

```julia
function gifdemo()
    anim = Animation()
    for i in 1:5
        frame(anim,plotquademo(i))
    end
    gif(anim,"approxanim.gif";fps=1)
end
gifdemo()
```
![](docs/approxanim.gif)
