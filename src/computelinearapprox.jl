""" 
Compute and add constraints best fitting the given input points.

input:
    * m: JuMP Model
    * x: JuMP variable, first axis of function to be approximated
    * d: array of points to be approximated
    * fd: function values of points d
    * K: number of constraints to be calculated - higher improves accuracy and computation time
    * opt: optimizer (e.g. Xpress.Optimizer or Cbc.Optimizer)
    * y: JuMP variable, second axis, if not given, a variable y will be added to model m
output:
    * y: reference to y-variable in model m

    convexlinearization(m, x, d, fd; K=3, opt, y=nothing)
"""
function convexlinearization(m, x, d, fd; K=3, opt,y=nothing)
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
Construct model to be used to compute best piecewise fit, using N segments

input: 
    * points: array of points to be approximated
    * K: number of segments to use (higher gives better accuracy and longer computation)
output:
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

function optimize_points(points,optimizer;numpoints=3)
    
    m = bestlinearization(points,numpoints)
    set_optimizer(m,optimizer)
    optimize!(m)

    c = JuMP.value.(m[:c])
    d = JuMP.value.(m[:d])

    return [d[i]=>c[i] for i in 1:length(c)]
end

function convex_pwlinear_interpolate(m, x, d, fd, opt; nseg=3, y=nothing)
    @assert length(d) == length(fd)
    pts = [d[i] => fd[i] for i = 1:length(d)]
    prs = interpolatepw_convex(pts, opt; nseg=nseg)
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
Finds a convex piecewise linear approximation to the given points
input:
    * points: array of points to be approximated
    * optimizer: the optimizer to be used
    * nseg: the nuber of segments in the interpolant
    * penalty: the measure for deviations

output: 
    * array of line segments given with  interception and slope
"""
function interpolatepw_convex(points, optimizer; nseg=5, penalty="l1")
    ipts = interpolatepw(points, optimizer; nseg=nseg, penalty=penalty)
    cpts = convexify(ipts, optimizer)
    return lines(cpts)
end


"""
Find a piecewise linear interpolant with a given number of breakpoints
input:
    * points: array of points to be approximated
    * optimizer: the optimizer to be used
    * nseg: the nuber of segments in the interpolant
    * penalty: the measure for deviations

output: 
    * array of interpolating points
"""
function interpolatepw(points, optimizer; nseg=5, penalty="l1")
    
    N = length(points)
    𝒩 = 1:N 

    X = first.(points)
    Y = last.(points)

    # Find gradients
    C = [(i < j ? (Y[j] -Y[i]) / (X[j]-X[i]) : 0)  for i ∈ 𝒩, j ∈ 𝒩]

    # Calculate penalties 
    if penalty == "l1"
        P = [(i < j ? sum(abs(C[i,j] * (X[k] - X[i]) + Y[i] - Y[k]) for k ∈ i:j) : 0) for i ∈ 𝒩, j ∈ 𝒩]
    elseif penalty == "l2"
        P = [(i < j ? sum((C[i,j] * (X[k] - X[i]) + Y[i] - Y[k])^2 for k ∈ i:j) : 0) for i ∈ 𝒩, j ∈ 𝒩]
    end
    
    # Setup optimization model
    m = Model()

    @variable(m, z[𝒩,𝒩], Bin)

    # Minimize total penalty
    @objective(m, Min, sum(P[i,j]*z[i,j] for i ∈ 𝒩, j ∈ 𝒩))

    # Number of line segments in interpolant
    @constraint(m, sum(z[i,j] for i ∈ 𝒩, j ∈ 𝒩) == nseg )

    # Only forward segments allowed
    for i ∈ 𝒩, j ∈ 𝒩 
        if i >= j 
            @constraint(m, z[i,j] == 0)
        end
    end 

    # Path modelling
    for j ∈ 𝒩
        @constraint(m, (j > 1 ? sum(z[i,j] for i ∈ 𝒩) : 1) == (j < N ? sum(z[j,i] for i ∈ 𝒩 ) : 1))
    end

    set_optimizer(m,optimizer)
    optimize!(m)
    
    zᴼᵖᵗ = JuMP.value.(m[:z])
    ipts = collect(X[i] => Y[i] for i ∈ 𝒩, j ∈ 𝒩 if zᴼᵖᵗ[i,j] == 1)
    push!(ipts, X[N] => Y[N])

    return ipts

end

"""
Convexify a one dimensional piecewise linear function by adjusting the y-values.

input:
    * points: array of breakpoints for the pwl
    * optimizer: the optimizer to be used

output: 
    * array of adjusted breakpoints
"""
function convexify(points, optimizer)
    
    N = length(points)
    X = first.(points)
    Y = last.(points)

    m = Model()

    @variable(m, d_pos[1:N] >= 0)
    @variable(m, d_neg[1:N] >= 0)

    @objective(m,Min,sum(d_pos[i] + d_neg[i] for i=1:N))
    for i=2:N-1
        @constraint(m, (Y[i] + d_pos[i] - d_neg[i] - Y[i-1] - d_pos[i-1] + d_neg[i-1]) / 
            (X[i] -X[i-1]) <= (Y[i+1] + d_pos[i+1] - d_neg[i+1] - Y[i] - d_pos[i] + d_neg[i]) / (X[i+1] -X[i]))
    end

    set_optimizer(m,optimizer)
    optimize!(m)

    dp = JuMP.value.(m[:d_pos])
    dn = JuMP.value.(m[:d_neg])

    return [X[i] => Y[i] + dp[i] - dn[i] for i=1:N]

end

function lines(points)
    N = length(points)
    X = first.(points)
    Y = last.(points)

    c = [(Y[i+1] - Y[i]) / (X[i+1] - X[i]) for i ∈ 1:(N-1)]
    d = [Y[i] - c[i]*X[i] for i ∈ 1:(N-1)]
    return [d[i] => c[i] for i ∈ 1:(N-1)]
end

function isconvex(points; tol=1e-5)
    c = last.(lines(points))
    for i=1:length(c)-1
        if c[i] > c[i+1] + tol 
            return false
        end
    end
    return true
end