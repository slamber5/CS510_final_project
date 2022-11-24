using Plots
include("../euler_methods.jl")

function test_nodrive()
    B(t) = [0 1;
            -1 0]
    f(y,t) = B(t)*y

    y0 = [1;0]
    T = 2*2*pi 
    dt = 0.5
    J(y,t) = B(t) 

    #=
    f(y,t) = -y
    y0 = [1]
    T = 4*pi 
    dt = 0.01
    J(y,t) = -1 
    =#
    t,y = trapezoidal(f,J,y0,T,dt)
    p = plot(title="Linearized trapezoidal method,dt=0.5",
    legend=:bottomright)
    xlabel!(p,"t")
    plot!(p,t,y[1,:],marker=:circle,label="y1")
    plot!(p,t,y[2,:],marker=:circle,label="y2")
    plot!(p,t,cos.(t),label="cos(t)")
    plot!(p,t,-sin.(t),label="-sin(t)")
    display(p)
    savefig("linearized_trapezoidal_dt_0p5.pdf")
end

function test_drive()
    B = [0 1;-1 0]
    f(y,t) = B*y + [0;sin(t)]

    y0 = [1;0]
    T = 50*2*pi 
    dt = 0.5
    J(y,t) = B 
    t,y = trapezoidal(f,J,y0,T,dt)
    p = plot(title="Linearized trapezoidal method,dt=0.5\n(problem with resonant drive)",
    legend=:bottomright)
    xlabel!(p,"t")
    plot!(p,t,y[1,:],marker=:circle,label="y1")
    plot!(p,t,y[2,:],marker=:circle,label="y2")
    plot!(p,t,cos.(t),label="cos(t)")
    plot!(p,t,-sin.(t),label="-sin(t)")
    display(p)
    savefig("linearized_trapezoidal_drive_dt_0p5.pdf")
end

function compare_euler_trapezoidal()
    A = [0 1; -1 0]
    f(t,y) = A*y
    g(y,t) = f(t,y)
    J(y,t) = A 
    y0 = [1;0]
    T = 4*pi
    dt = 0.1
    tt,yt = trapezoidal(g,J,y0,T,dt)
    te,ye = forward_euler(f,y0,T,dt)
    t = range(0,T,1000)
    y = [cos.(t) -sin.(t)]'
    p = plot(title = "Forward Euler versus Trapezoidal")
    plot!(p,tt,yt',label="Trapezoidal",marker=:circle)
    plot!(p,te,ye',label="Forward Euler",marker=:square)
    plot!(p,t,y',label="Exact")
    display(p)
end

#test_drive()
compare_euler_trapezoidal()