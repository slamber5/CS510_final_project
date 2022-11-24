using Plots 

include("../euler_methods.jl")
pwd()
######### Problem 1 ##########
A = -1
y0 = 1
f(t,y) = A.*y
T = 3
dt = 0.1 # max stable timestep is -2/A = 2

# Forward Euler solution 
tf,yf = my_forward_euler(f,y0,T,dt)

# Backward Euler solution 
tb,yb = my_backward_euler(A,y0,T,dt)

# Exact solution 
tex = 0:0.01:T
yex = y0.*exp.(A.*tex)

# Plots 
p = plot(tex,yex,label="Exact")
scatter!(tf,yf,label = "Forward Euler",marker=:circle)
title!("Forward Euler and Exact Solution")
display(p)

p = plot(tex,yex,label="Exact")
scatter!(tb,yb,label = "Backward Euler",marker=:circle)
title!("Backward Euler and Exact Solution")
display(p)
