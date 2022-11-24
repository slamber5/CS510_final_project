include("euler_methods.jl")
include("my_solvers.jl")
using LinearAlgebra

"""
    burgers_euler_periodic(u0,T,dt)
Solve the Burgers equation on the spatial inverval [0,1]
with periodic boundary conditions using forward Euler.
u0 = initial condition with u0[1] = u0(0) and 
u0[N] = u0(1)
"""
function burgers_euler_periodic(u0,T,dt)
    Nx = size(u0)[1]-1
    dx = 1/Nx 
    x = dx.*(0:Nx)
    D = symmetric_D_periodic(Nx+1,dx)
    f(t,u) = -(1/2).*D*(u.*u)
    t, u = forward_euler(f,u0,T,dt)
    return x, t, u
end


"""
    burgers_trapezoidal_periodic(u0,T,dt)
Solve the Burgers equation on the spatial inverval [0,1]
with periodic boundary conditions using the linearized
trapezoidal method. 
u0 = initial condition with u0[1] = u0(0) and 
u0[N] = u0(1)
"""
function burgers_trapezoidal_periodic(u0,T,dt)
    Nx = size(u0)[1]-1
    dx = 1/Nx 
    x = dx.*(0:Nx)
    D = symmetric_D_periodic(Nx+1,dx)
    f(u,t) = -(1/2).*D*(u.*u)
    J(u,t) = -D*Diagonal(u)
    t,u = trapezoidal(f,J,u0,T,dt)
    return x, t, u 
end


"""
    burgers_upwind_conservative(u0,T,dt)
Solve the Burgers equation with the upwind conservative scheme
and Dirichlet boundary conditions. Equivalent to forward Euler with
an asymmetric first derivative matrix. Equivalent to forward Euler 
timestepping with an asymmetric difference derivative operator. ul = 
upwind boundary (upwind direction assumed to be at negative x 
direction)
"""
function burgers_upwind_conservative(u0,T,dt)
    Nx = size(u0)[1]-1
    dx = 1/Nx 
    x = dx.*(0:Nx)
    D = upwind_D_Dirichlet(Nx,dx)
    f(t,u) = [0;(-1/2).*D*(u.*u)]
    t,u = forward_euler(f,u0,T,dt)
    return x,t,u 
end


"""
    burgers_moc(u0,T)
Solve the Burgers equation with the method of characteristics, given
an initial condition vector u0 whose components are values of 
the initial condition function at evenly spaced points between
0 and 1 (inclusive). Returns a vector of x coordinates and 
corresponding u values (the x and u values at time T after the initial
time).
"""

function burgers_moc(u0,T,dt)
    Nx = size(u0)[1]-1
    Nt = Int(ceil(T/dt))
    t = dt.*(0:Nt)
    dx = 1/Nx 
    x0 = dx.*(0:Nx)
    x = x0 .+ u0*t'
    return x,t,u0
end

ic_tanh(x;a=10,b=1/4) = 1/2 - (1/2)*tanh(a*(x - b))


