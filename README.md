# PiecewiseLinearApprox

Helper functions to add approximations of convex functions to optimization models modelled in [JuMP](https://jump.dev/). Currently supporting univariate monotonically increasing functions. This package creates and solves a MILP to fit a set of points, and adds the resulting linear constraints to the optimization model. This method is based on [Toriello & Vielma, 2012](https://doi.org/10.1016/j.ejor.2011.12.030), following the formulation/notation used in [Kazda & Li, 2018](https://doi.org/10.3390/pr6100198)

For non-convex functions, consider using [PiecewiseLinearOpt.jl](https://github.com/joehuchette/PiecewiseLinearOpt.jl).

## Usage

```julia
using JuMP, PiecewiseLinearApprox, Xpress
d = range(0,1,step=0.1)
fd = [i^2 for i in d]

m = Model()
@variable(m, x)
@variable(m, test_y)
# Compute and add constraints approximating f(d)
y = convexlinearization(m, x, d, fd; opt=Xpress.Optimizer ,K=5, y=test_y)
# Minimize
@objective(m, Min, y)
set_optimizer(m, Xpress.Optimizer)
# Check approximation/solution at x = 0.5
@constraint(m, x >= 0.5)
optimize!(m)
value(m[:test_y]) # 0.25153374233128833
```

To keep dependencies light, PiecewiseLinearApprox does not include plotting of the resulting approximation. Plotting can be useful, however, the following demonstrates how this can be achieved:

```julia
using PiecewiseLinearApprox,Xpress,Plots

function plotquademo(N=3,opt=Xpress.Optimizer)
    points = [i=>i^2 for i in 0:0.1:1]
    m = bestlinearization(points,N)
    set_optimizer(m, opt)
    optimize!(m)
    c = JuMP.value.(m[:c])
    d = JuMP.value.(m[:d])
    x = [p.first for p in points]
    y = [p.second for p in points]
    p = plot(x,y,markershape=:x,label="Measurements",legend=:topleft,ylims=(-0.5,1))
    for i = 1:length(c)
        plot!(p,[0,1],[d[i],d[i]+c[i]],label="Constraint $i")
    end
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
    gif(anim,"approxanim.gif";fps=0.3)
end
```
![](docs/approxanim.gif)



[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
