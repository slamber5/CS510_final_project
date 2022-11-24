include("../euler_methods.jl")
include("../my_solvers.jl")
include("../heat_MOL.jl")

using Plots 

# Solve the heat equation with simple initial 
# and boundary conditions and a simple source term


F(x,t) .= 0
function f(x)
    return sin.(x)
end

function plotHeatSol(f,F,k,X,T,N)
    ul = 0
    ur = 0
    dx = X/N
    dt  = 0.5*dx^2/(2*k) # Half the maximum stable timestep
    n = Int(ceil(T/dt))
    x,t,u = heat1d(f,F,k,ul,ur,X,T,dt,N)

    p = scatter(x,u[:,1])
    scatter!(p,x,u[:,Int(ceil(0.5*n))])
    scatter!(p,x,u[:,n])

    ndense = 1000
    tdense = T.*(0:1/ndense:1)
    xdense = X.*(0:1/ndense:1)
    uexact = sin.(xdense)*exp.(-1 .* tdense)'
    plot!(p,xdense,uexact[:,1])
    plot!(p,xdense,uexact[:,Int(ceil(0.5*ndense))])
    plot!(p,xdense,uexact[:,ndense])
    display(p)
    # Is fstep returning all zeros? no
    # Is forward Euler on a single point returning a constant?
    plot(t,u[50,:])
    # No, but u is decreasing very slowly
    # Is the timestep set incorrectly somewhere?
    # Would appear that no, it isn't 
    # Unclear what to do so I'll try another function
end


function g(x)
    return x.*(1 .-x)
end

# That didn't work either
# Still isn't evolving quickly enough
# If I increase the final time does it appear 
# to be doing the right thing?
# Well this is strange-- it isn't complaining,
# but T isn't set anywhere in this file
# Let's fix that and set it to 2 explicitly

#plotHeatSol(f,F,1,pi,2,100)

# Didn't change anything-- predictably
# Set it to 10 and see if anything happens:

#plotHeatSol(f,F,1,pi,10,100)
#plotHeatSol(f,F,1,pi,1000,100)

# Ok, so it seems to be doing the right thing, but just
# far too slowly
# Aha
# Literally Brittany talked about this in class 
# You have to divide by dx^2 in the derivative matrix 

# plotHeatSol(f,F,1,pi,2,10)

# Lovely 



