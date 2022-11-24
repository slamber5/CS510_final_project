include("my_solvers.jl")
include("euler_methods.jl")

function heat1d(f,F,k,ul,ur,X,T,dt,N)
    """
        heat1d(F,k,N,dt)
    Solve the 1d heat equation 
    
    u_t(x,t) = k*u_xx(x,t) + F(x,t)

    on the spatial interval 
    [0,X] for the time interval [0,T] with fixed
    boundary conditions u(0,t) = ul and u(X,t) = ur
    and initial condition f(x)
    
    -Initial and boundary conditions are assumed to be 
    compatible

    N = number of subdivisions of spatial interval 
    N+1 = number of solution points
    dt = timestep

    """
    # space array
    dx = X/N
    x = dx.*(0:N)
    # Number of time intervals 
    n = Int(ceil(T/dt))
    # Solution array 
    u = zeros(N+1,n+1)
    # Initial condition 
    u[:,1] = f.(x)
    # Boundary conditions
    u[1,:] .= ul 
    u[end,:] .= ur  

    # Centered difference matrix 
    dd = ddMatrix(N)

    # Solve the equation with forward Euler time integration
    # Bracketing y between BCs enables arbitrary constants
    # to be chosen as boundary values 
    fstep(t,y) = k.*dd*[ul;y;ur]./dx^2 .+ F(x[2:end-1],t)
    y0 = u[2:end-1,1]
    t,y = forward_euler(fstep,y0,T,dt)
    u[2:end-1,:] = y  
    return x, t, u 
end


