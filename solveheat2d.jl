include("my_solvers.jl")

"""
    heat2d(k,u0,ut,ub,ul,ur,F,lx,ly,nx,ny,dt,T)
Solve the 2d heat equation
k = thermal conductivity
u0 = array of initial conditions 
ut,ub,ul,ur: vectors of boundary conditions 
    (t <-> top <-> +y, r <-> right <-> +x)
    - Note that these boundary conditions must be compatible with each
    other and the initial condition 
    - x indexed by row (first index), y by column (second index)
F = source term (F(x,y,t))
lx,ly: side lengths of solution region
nx,ny: number of subintervals in each direction 
dt = timestep 
T = solution final time (initial time 0)
"""

function heat2d(k,u0,ut,ub,ul,ur,F,lx,ly,nx,ny,dt,T)
    # Allocate array for solution and set initial conditions 
    nt = Int(ceil(T/dt))
    t = 0:dt:nt*dt
    u = zeros(nx+1,ny+1,nt+1)
    u[:,:,1] = u0 
    u[1,:,:] .= ul
    u[end,:,:] .= ur 
    u[:,1,:] .= ub 
    u[:,end,:] .= ut  
    
    # Derivative matrix
    ddx = ddMatrix(nx)
    ddy = ddMatrix(ny)

    # Evaluate solution 
    for i = 1:nt
        u[2:end-1,2:end-1,i+1] = u[2:end-1,2:end-1,i] + 
        k * dt .* ((1/dx*2) .* )  ##### You were here #####
    end
    

    
end