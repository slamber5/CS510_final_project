using Plots 
using Printf 

include("../my_solvers.jl")


function plotCGSolTime(Ns,eps=1e-6,maxiter=1000;nruns=1)
    """
        plotCGSolTime(Ns)
    plot the time required to solve a linear system
    with the conjugate gradient method. 
    Ns = list of matrix sizes
    eps = tolerance  
    nruns: number of trial runs 
    """
    p = plot(legend=false)
    title!("CG performance at various matrix sizes")
    xlabel!("Matrix size N")
    ylabel!("Solution time (s)")
    for i = 1:nruns
        @printf("run number %d\n",i)
        t = zeros(size(Ns))
        @printf("Solving: N = %d\n",N)
        for (i,N) in enumerate(Ns)
            A = spdrmat(N)
            b = rand(N)
            x0 = zeros(N)
            t[i] = @elapsed conj_grad(A,x0,b,eps,maxiter)
            @printf("Solved with N = %d in %.4e seconds\n",
                N,t[i])
        end
        Plots.scatter!(p,Ns,t)
    end
    display(p)
end

function testCGAccuracy(Ns,eps=1e-6,maxiter=1000;nruns=3)
    """
        testCGAccuracy(Ns)
    Test the accuracy of the CG solver for various 
    matrix sizes N. Record the relative error and plot 
    it as a function of N.
    Ns: list of matrix sizes
    eps: tolerance
    nruns: number of repetitions
    """
    p = plot(legend=false)
    title!(p,"Relative error of CG solver")
    xlabel!("N (matrix size)")
    ylabel!("Relative error")
    @printf("\n")
    for i = 1:nruns
        @printf("Run number %d\n",i)
        errs = zeros(size(Ns))
        for (i,N) in enumerate(Ns)
            x0 = zeros(N)
            A = spdrmat(N)
            xtrue = rand(N)
            b = A*xtrue
            x = conj_grad(A,x0,b,eps,maxiter)
            errs[i] = norm(x-xtrue)/norm(xtrue)
            @printf("Solved with N = %d\n",N)
        end
        Plots.scatter!(p,Ns[2:end],errs[2:end])
    end
    display(p)
end

# Solve with LU factorization and CG and display results
N = 5
B = rand(N,N)
A = transpose(B)*B + identity(N)
b = rand(N)

xLU = solvePLU(A,b)

x0 = zeros(N)
eps = 1e-6
maxiter = 1000
xCG = conj_grad(A,x0,b,eps,maxiter)

@printf("LU factorization solution: \n")
display(xLU)
@printf("CG solution: \n")
display(xCG)


# Time CG factorization and plot performance
Ns = 1:100:2001
#plotCGSolTime(Ns,nruns=3)
testCGAccuracy(Ns)

