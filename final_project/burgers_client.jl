using Plots 
using Printf
include("../burgers_solvers.jl")

function plot_burgers_euler_periodic(T,nplot;save=false)
    x = 0:0.01:1
    u0 = sin.(pi.*x).^2
    #u0 = zeros(length(x))
    #u0[1:div(end,2)] .= 1
    dt = 0.01
    Nt = Int(ceil(T/dt))
    x,t,u = burgers_euler_periodic(u0,T,dt)
    p = plot(x,u[:,1],
    label=@sprintf("T = %.2f",t[1]),
    legend =:topleft,title = 
    @sprintf("Forward Euler timestepping,dt = %.2f",dt))
    xlabel!(p,"x")
    ylabel!(p,"u")
    n = 7:7:35
    for i in nplot
        plot!(p,x,u[:,i],label=@sprintf("T = %.2f",t[i]))
    end
    display(p)
    if save
        savefig("burgers_euler_periodic.pdf")
    end
end

function plot_burgers_trapezoidal_periodic(T,nplot;save=false)
    x = 0:0.01:1
    u0 = sin.(pi.*x).^2
    dt = 0.01
    Nt = Int(ceil(T/dt))
    x,t,u = burgers_trapezoidal_periodic(u0,T,dt)
    p = plot(x,u[:,1],
    label=@sprintf("T = %.2f",t[1]),
    legend =:topleft,title = 
    @sprintf("Trapezoidal timestepping,dt = %.2f",dt))
    xlabel!(p,"x")
    ylabel!(p,"u")
    #=
    for j = 1:ncurves
        plot!(p,x,u[:,j*div(Nt,ncurves)],label=@sprintf("T = %.2f",t[j*div(Nt,ncurves)]))
    end
    =#
    for i in nplot
        plot!(p,x,u[:,i],label=@sprintf("T = %.2f",t[i]))
    end
    display(p)
    if save
        savefig("burgers_trapezoidal_periodic.pdf")
    end
end

function plot_burgers(;T=1,nplot=10:10:50,dt=0.01,dx=0.01,
    ic=ic_tanh,save=false,method="upwind_conservative",
    legpos=:topright)
    x = 0:0.01:1
    Nt = Int(ceil(T/dt))
    t = dt.*(0:Nt)
    u0 = ic.(x)
    if method == "upwind_conservative"
        x,t,u = burgers_upwind_conservative(u0,T,dt)
        tstr = "Upwind conservative, Dirichlet left BC"
    elseif method == "trapezoidal_periodic"
        x,t,u = burgers_trapezoidal_periodic(u0,T,dt)
        tstr = "Trapezoidal symmetric, periodic BC"
    elseif method == "euler_periodic"
        x,t,u = burgers_euler_periodic(u0,T,dt)
        tstr = "Forward Euler symmetric, periodic BC"
    elseif method == "moc"
        x,t,u = burgers_moc(u0,T,dt)
        tstr = "Method of characteristics"
    else
        @printf("\nWarning: invalid solver code")
    end
    p = plot(legend =legpos,title = 
        @sprintf("%s, dt = %.2f",tstr,dt))
    xlabel!(p,"x")
    ylabel!(p,"u")
    if method != "moc"
        plot!(p,x,u[:,1],label=@sprintf("t = %.2f",t[1]))
        for i in nplot
            plot!(p,x,u[:,i],label=@sprintf("t = %.2f",t[i]))
        end
    else
        plot!(p,x[:,1],u,label=@sprintf("t = %.2f",t[1]))
        for i in nplot
            plot!(p,x[:,i],u,label=@sprintf("t = %.2f",t[i]))
        end
    end
    display(p)
    if save
        savefig(@sprintf("%s.pdf",tstr))
    end
end


#plot_burgers_euler_periodic(50,10:10:50,false)
#plot_burgers_trapezoidal_periodic(50,10:10:50,false)

#plot_burgers(50,10:10:50;save=true)
#plot_burgers(50,10:10:50;save=true,method="moc")
ic = x->exp(-(10*(x-0.25))^2)
plot_burgers(save=false,nplot=(11:10:101),method="moc",ic = ic)
plot_burgers(save=false,nplot=(11:10:101),method="upwind_conservative",ic = ic)


#### YOU WERE HERE #####
# Overlay method of characteristics solution on numerical solution(s)