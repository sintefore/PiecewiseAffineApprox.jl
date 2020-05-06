""" 
```
convexlinearization(m, x, d, fd; K, opt, y)
````
Compute and add constraints best fitting the given input points.

## Input

* m: JuMP Model
* x: JuMP variable, first axis of function to be approximated
* d: array of points to be approximated
* fd: function values of points d
* K: number of constraints to be calculated - higher improves accuracy and computation time
* opt: optimizer (e.g. Xpress.Optimizer or Cbc.Optimizer)
* y: JuMP variable, second axis, if not given, a variable y will be added to model m

## Output

* y: reference to y-variable in model m

"""
function convexlinearization(m, x, d, fd; K=3, opt=Cbc.Optimizer,y=nothing)
    @assert length(d) == length(fd)
    pts = [d[i] => fd[i] for i = 1:length(d)]

    prs = optimize_points(pts,opt;numpoints=K)
    # TODO: Check variable name not already used in m (will error)
    if isnothing(y)
        @variable(m, y) 
    end
    for p ∈ prs
        @constraint(m, y >= p.first + p.second * x )
    end
    return y
end

"""
```
bestlinearization(points, K)
```
Construct model to be used to compute best piecewise fit, using K segments

## input

* points: array of points to be approximated
* K: number of segments to use (higher gives better accuracy and longer computation)

## output

* m: JuMP model
"""
function bestlinearization(points,K=3)
    
    I = length(points) 
    M_big = maximum((i.second for i in points))

    m = Model()

    @variable(m, t)
    @variable(m, c[1:K])
    @variable(m, d[1:K])
    @variable(m, z[1:I, 1:K],Bin)
    @variable(m, f[1:I])

    @objective(m, Min, t)

    # s.t. 
    for (i,p) in enumerate(points)
        x_i = p.first
        y_i = p.second

        @constraint(m, abs(y_i)*t >= f[i] - y_i)
        @constraint(m, abs(y_i)*t >= y_i - f[i])
        for k ∈ 1:K
            @constraint(m, f[i] >= c[k]*x_i + d[k])
            @constraint(m, f[i] <= c[k]*x_i + d[k] + M_big*(1-z[i,k]))
        end
        @constraint(m, sum(z[i,k] for k = 1:K) == 1)
    end
    return m
end

@memoize function optimize_points(points,optimizer=Cbc.Optimizer;numpoints=3)
    
    m = bestlinearization(points,numpoints)
    set_optimizer(m,optimizer)
    optimize!(m)

    c = JuMP.value.(m[:c])
    d = JuMP.value.(m[:d])

    return [d[i]=>c[i] for i in 1:length(c)]
end
