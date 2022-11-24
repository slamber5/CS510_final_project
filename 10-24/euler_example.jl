using Plots 

# Write forward Euler to solve an affine IVP 
# y' = Ay + b on 0 ≤ t ≤ Tf
# with y(0) = y0 

# Example system: t0 = 0, A = [-3 13; -5 -1]
# y0 = [3; -10], b = [0; 0]

t0 = 0
y0 = [3; 10]
Tf = 3;
dt = 0.1;
N = Integer(Tf/dt)
b = [0;0]
A = [-3 13; -5 -1]

y = zeros(2,N+1)
t = zeros(N+1)
y[:,1] = y0
for n = 1:N
    y[:,n+1] = y[:,n] + dt.*(A*y[:,n]+b)
    t[n+1] = t[n]+dt
end
# plot initial condition 
plot([y[1,1]],[y[2,1]],marker=(:circle,5), color = :blue, legend = false, xlims=(-100,100),ylims=(-60,60))
# plot solution 
plot!(y[1,:],y[2,:])
# animation
for n = 2:N+1
    p = plot!([y[1,n]],[y[2,n]],marker=(:circle,5), color = :blue, legend = false, xlims=(-60,60),ylims=(-60,60))
    display(p)
    sleep(0.05)
end