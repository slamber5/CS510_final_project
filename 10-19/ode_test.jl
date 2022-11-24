using Plots 

# solve the IVP y' = -3y with initial condition 
# y(0) = 17 on the time interval 0 ≤ t ≤ Tf with exact
# solution y(t) = y0*exp(-3t)

#=
Tf = 4
dt = 0.01
N = Integer(Tf/dt)
t = 0:dt:Tf 
t = Vector{Float64}(undef,N+1)
# Could also do t = zeros(N+1)
y = Vector{Float64}(undef,N+1)

t[1] = 0
y[1] = 17

for n = 1:N
    y[n+1] = y[n] + dt*(-3).*y[n]
    t[n+1] = t[n] + dt
end

scatter(t,y,label="approximation",
    shape=:circle,color=:green)

tfine = 0:dt/20:Tf
yexact = 17 .*exp.(-3 .*tfine)
plot!(tfine,yexact,label="exact")
=#

function my_forward_euler(t0,Tf,dt, y0, f)
    N = Integer(Tf/dt)
    t = 0:dt:Tf 
    t = Vector{Float64}(undef,N+1)
    # Could also do t = zeros(N+1)
    y = Vector{Float64}(undef,N+1)

    t[1] = t0
    y[1] = y0

    for n = 1:N
        y[n+1] = y[n] + dt.*f(t[n],y[n])
        t[n+1] = t[n] + dt
    end
    return (t,y)
end

function test_equationRHS(t,y)
    return λ*y
end


