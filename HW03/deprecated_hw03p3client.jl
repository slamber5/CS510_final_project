include("../euler_methods.jl")
include("../my_solvers.jl")
include("../heat_MOL.jl")

using Plots 

# Solve the heat equation with simple initial 
# and boundary conditions and a simple source term

k = 2
F(x,t) = (k*pi^2 - 2)*exp(-2*t).*sin.(pi.*x)
f(x) = sin.(pi.*x)
ul = 0
ur = 0
X = 1
T = 2
N = 100
dx = X/N
dt  = 0.5*dx^2/(2*k) # Half the maximum stable timestep
n = Int(ceil(T/dt))
x,t,u = heat1d(f,F,k,ul,ur,X,T,dt,N)

p = scatter(x,u[:,1])
scatter!(p,x,u[:,Int(ceil(0.5*n))])
scatter!(p,x,u[:,n])

ndense = 1000
tdense = 0:1/ndense:T
xdense = 0:1/ndense:X
uexact = sin.(pi.*xdense)*exp.(-2 .* tdense)'
plot!(p,xdense,uexact[:,1])
plot!(p,xdense,uexact[:,Int(ceil(0.5*ndense))])
plot!(p,xdense,uexact[:,ndense])


