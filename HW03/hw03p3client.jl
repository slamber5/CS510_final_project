include("../euler_methods.jl")
include("../my_solvers.jl")
include("../heat_MOL.jl")

using Plots 
using Printf

# Solve the heat equation with simple initial 
# and boundary conditions and a simple source term
function F(x,t)
    return 2*(pi^2-1)*exp(-2*t).*sin.(pi.*x)
end

function f(x)
    return sin.(pi.*x)
end


function plotHeatSol(f,F,k,X,T,N)
    """
        plotHeatSol(f,F,k,X,T,N)
    Plot numerical and analytical solutions to the heat equation.
    """
    ul = 0
    ur = 0
    dx = X/N
    dt  = 0.5*dx^2/(2*k) # Half the maximum stable timestep
    n = Int(ceil(T/dt))
    x,t,u = heat1d(f,F,k,ul,ur,X,T,dt,N)

    p = scatter(x,u[:,1],label=@sprintf("N, t = %.1f",t[1]))
    midind = Int(ceil(0.5*n))
    scatter!(p,x,u[:,midind],label = @sprintf("N, t = %.1f",t[midind]))
    scatter!(p,x,u[:,n],label = @sprintf("N, t = %.1f",t[end]))

    ndense = 1000
    tdense = T.*(0:1/ndense:1)
    xdense = X.*(0:1/ndense:1)
    uexact = sin.(pi.*xdense)*exp.(-2 .* tdense)'
    plot!(p,xdense,uexact[:,1],label=@sprintf("A, t = %.1f",tdense[1]))
    midind_dense = Int(ceil(0.5*ndense))
    plot!(p,xdense,uexact[:,midind_dense],
        label=@sprintf("A, t = %.1f",tdense[midind_dense]))
    plot!(p,xdense,uexact[:,ndense],
        label=@sprintf("A, t = %.1f",tdense[end]))
    display(p)
    title!(p,"Numerical and analytical solutions\n" * 
    "(N: numerical; A: analytical)")
    ylabel!(p,"u(x,t)")
    xlabel!(p,"x")
    savefig(p,"heat1dsol.pdf")
end

function plotTimesteps(f,F,k,X,T,N)
    """
        plotTimesteps(f,F,k,X,T,N)
    Generate plots of solutions of the heat equation for different
    values of the timestep.
    """
    ul = 0
    ur = 0
    dx = X/N
    tfacts = [0.5;1;1.083]
    for (i,tf) in enumerate(tfacts)
        dt  = tf*dx^2/(2*k) # Half the maximum stable timestep
        n = Int(ceil(T/dt))
        x,t,u = heat1d(f,F,k,ul,ur,X,T,dt,N)

        p = plot(x,u[:,1],label=@sprintf("N, t = %.1f",t[1]),marker=:circle)
        midind = Int(ceil(0.5*n))
        plot!(p,x,u[:,midind],label = @sprintf("N, t = %.1f",t[midind]),
            marker=:circle)
        plot!(p,x,u[:,n],label = @sprintf("N, t = %.1f",t[end]),
            marker=:circle)

        ndense = 1000
        tdense = T.*(0:1/ndense:1)
        xdense = X.*(0:1/ndense:1)
        uexact = sin.(pi.*xdense)*exp.(-2 .* tdense)'
        plot!(p,xdense,uexact[:,1],label=@sprintf("A, t = %.1f",tdense[1]))
        midind_dense = Int(ceil(0.5*ndense))
        plot!(p,xdense,uexact[:,midind_dense],
            label=@sprintf("A, t = %.1f",tdense[midind_dense]))
        plot!(p,xdense,uexact[:,ndense],
            label=@sprintf("A, t = %.1f",tdense[end]))
        display(p)
        title!(p,"Numerical and analytical solutions\n" * 
        @sprintf("dt = %.3f x dt_max",tf))
        ylabel!(p,"u(x,t)")
        xlabel!(p,"x")
        savefig(p,@sprintf("tfact%d.pdf",i))
    end
end

function L2norm(v,dx)
    """
        Compute the discrete L2 norm with point spacing dx.
    """
    return sqrt(dx)*norm(v)
end

function errsdx(f,F,k)
    """
        errsdx(f,F,k,X,T,N)
    Compute the L2 error for different values of dx.
    """
    ul = 0
    ur = 0
    lambda = 0.1
    X = 1
    T = 0.1
    errs = zeros(4)
    Ns = [10,20,40,80]
    for (i,N) in enumerate(Ns)
        dx = X/N 
        dt = lambda*dx^2
        x,t,u = heat1d(f,F,k,ul,ur,X,T,dt,N)
        uexact = exp(-2*T).*sin.(pi.*x)
        errs[i] = L2norm(u[:,end]-uexact,dx)
    end
    ratios = errs[1:end-1]./errs[2:end]
    rates = log2.(ratios)
    return errs, ratios, rates  
end

# Function call for part b:
plotHeatSol(f,F,2,1,1,10)

# Function call for part c:
plotTimesteps(f,F,2,1,1,10)

# Function calls for part d:
errs,ratios,rates = errsdx(f,F,2)
display(map(x -> @sprintf("%.3e",x), errs))
display(map(x -> @sprintf("%.3e",x), ratios))
display(map(x -> @sprintf("%.3e",x), rates))








