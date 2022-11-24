using Plots
using LinearAlgebra

function my_forward_euler(f,y0,T,dt)
"""
    my_forward_euler(f,y0,tf,dt)
Solve the ODE y'(t) = f(t,y) with initial condition y(0) = y0
on the time interval [0,T]
"""

    N = Integer(ceil(T/dt))
    t = 0:dt:N*dt
    y = zeros(N+1)  
    y[1] = y0

    for i = 1:N
        y[i+1] = y[i] .+ f(t[i],y[i]) .* dt
    end
    return t,y 
end

function my_backward_euler(A,y0,T,dt)
"
    my_backward_euler(A,y0,T,dt)
Solve the linear ODE y' = Ay via the backward Euler
method with initial condition y(0) = y0 over the
time interval [0,T]
"
    N = Integer(ceil(T/dt))
    t = 0:dt:N*dt
    y = zeros(N+1)  
    y[1] = y0
    for i = 1:N
        y[i+1] = (1/(1-A*dt))*y[i]
    end
    return t,y
end

function forward_euler(f,y0,T,dt)
    """
        forward_euler(f,y0,T,dt)
    solve the system y' = f(t,y) with the forward Euler method
    on the interval [0,T]. Assumes y is a vector. Generalization
    of my_forward_euler to multidimensional systems
    """
    N = size(y0)[1]
    n = Int(ceil(T/dt))
    y = zeros(N,n+1)
    t = dt.*(0:n)
    y[:,1] = y0 
    for i = 1:n
        y[:,i+1] = y[:,i] + dt.*f(t[i],y[:,i])
    end
    return t,y
end

#=
A = [0 1; -1 0]
f(t,y) = A*y
y0 = [1;0]
T = 6*pi
dt = 0.01
t,y = forward_euler(f,y0,T,dt)
Plots.plot(t,y')
=#

"""
    trapezoidal(f,J,y0,T,dt)
Solve y' = f(y,t) with the linearized trapezoidal 
method on the interval [0,T]

y0 = initial condition 
dt = timestep 
J = J(y,t) = Jacobian of f
"""

function trapezoidal(f,J,y0,T,dt)
    N = size(y0)[1]
    n = Int(ceil(T/dt))
    t = dt.*(0:n)
    y = zeros(N,n+1)
    y[:,1] = y0
    for i = 1:n
        A = I-(dt/2).*J(y[:,n],t[n+1])
        b = f(y[:,i],t[i+1])+f(y[:,i],t[i])
        y[:,i+1] = y[:,i] + (dt/2).*(A\b) 
    end
    return t,y
end