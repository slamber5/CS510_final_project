using Plots
using CurveFit
using Printf
using SparseArrays

function backsub(U,b)
"""
    backsub(U,b)
Solve Ux = b via back-substitution and return the solution x.

Notes:
U must be an upper-triangular matrix of size N with nonzero elements 
along the diagonal; b must be a column vector of length N 
"""
    
    N = size(b)[1]
    x = zeros(N,1)
    for i = 0:N-1
        x[N-i] = (1/U[N-i,N-i])*(b[N-i]-sum(U[N-i,N-i+1:N].*x[N-i+1:N]))
    end
    return x
end

function forwardsub(L,b)
"""
    forwardsub(L,b)
Solve Lx = b via forward-substitution. 
L must be a lower-triangular matrix of size NxN with nonzero elements 
along the diagonal; b must be a column vector of length N
"""
    N = size(b)[1]
    x = zeros(N,1)
    for i = 1:N
        x[i] = (1/L[i,i])*(b[i]-sum(L[i,1:i-1].*x[1:i-1]))
    end
    return x
end

function pivot_matrix(A,k)
"""
    pivot_matrix(A,k)
Compute the permutation matrix P such that left-multiplying A by 
P puts the largest-magnitude element in the kth column at or below the
diagonal on the diagonal

returns: permutation matrix P 
"""
    N = size(A)[1]
    # Get elements of kth column on and below diagonal
    col = A[:,k]
    # Zero the first k-1 elements of col so that they won't be picked
    # for the pivot
    col[1:k-1] .= 0 
    # Find index j of maximum element
    j = findmax(abs.(col))[2][1]
    # Create permutation matrix swapping rows j and k 
    P = identity(N)
    P[k,:],P[j,:] =  P[j,:],P[k,:]
    return P 
end

#=
function identity(N)
"""
DEPRECATED; use imat(N) instead.
    identity(N)
return the NxN identity matrix.
"""
    I = zeros(N,N)
    for i = 1:N
        I[i,i] = 1
    end
    return I 
end
=#

function imat(N)
    """
        imat(N)
    return the NxN identity matrix.
    """
        I = zeros(N,N)
        for i = 1:N
            I[i,i] = 1
        end
        return I 
    end

function Lstep(U, k)
"""
    Lstep(U, k)
Compute the matrix Ls required to eliminate the below-diagonal elements of 
column k from the matrix U via the multiplication Ls*U
Assumes U[k,k] != 0
"""
    N = size(U)[1]
    # Vector to be inserted into I to form Lk:
    col = zeros(N,1) 
    col[k] = 1
    col[k+1:N] = -U[k+1:N,k]./U[k,k]
    Ls = identity(N)
    Ls[:,k] = col
    return Ls
end

function LstepInv(Ls,k)
"""
    LstepInv(Lstep, k)
Compute the inverse of the row-reduction matrix Lstep given
that the column containing the elimination coefficients is k.
Lstep must be of the form returned by Lstep(), ie an identity
matrix but containing a single column with possibly-nonzero 
coefficients below the diagonal.
"""
    N = size(Ls)[1]
    LsInv = identity(N)
    LsInv[k+1:N,k] = -Ls[k+1:N,k]
    return LsInv
end

function computeLU(A)
"""
    computeLU(A)
compute the LU factorization of a matrix A.
A must be a square matrix that admits an LU factorization (ie, 
must be a matrix such zero-valued pivot elements don't appear
during Gaussian elimination).
"""
    N = size(A)[1]
    L = identity(N)
    U = A
    for k = 1:N
        Ls = Lstep(U,k)
        U = Ls*U
        L = L*LstepInv(Ls,k)
    end
    return L,U
end

function computePLU(A)
"""
    computePLU(A)
Determine P, L, and U such that PA = LU where P is a permutation matrix,
L is lower-triangular, and U is upper-triangular. A must be a square
matrix. 
"""
    N = size(A)[1]
    U = A
    L = identity(N)
    P = identity(N)
    for k = 1:N-1
        Ps = pivot_matrix(U,k)
        U = Ps*U
        Ls = Lstep(U,k)
        U = Ls*U 
        P = Ps*P 
        L = Ps*L*Ps*LstepInv(Ls,k)
    end
    return P, L, U 
end

function testPLU(N)
"""
    testPLU(N)
Test the accuracy of the PLU factorization algorithm. 
Generates a random NxN matrix B, computes the PLU factorization
of transpose(B)*B + I using computePLU, and finds the 
square of the matrix norm of the difference of PA and LU  
"""
    B = rand(N,N)
    A = transpose(B)*B+identity(N) 
    P,L,U = computePLU(A)
    return sum((P*A - L*U).^2)
end

function solvePLU(A,b)
"""
    solvePLU(A,b)
Solves Ax = b via PLU factorization. A is NxN, b is Nx1
"""
    P,L,U = computePLU(A)
    y = forwardsub(L,P*b)
    x = backsub(U,y)
    return x    
end

function testsolvePLU(N)
    """
        testPLU(N)
    Test the accuracy of solvePLU. 
    Solves Ax = b with solvePLU and outputs the squared error |Ax - b|^2. 
    """
        B  = rand(N,N)
        A = transpose(B)*B + identity(N)
        b = rand(N,1)
        x = solvePLU(A,b)
        return sum((A*x - b).^2)
    end

function showPLUExample(N)
"""
    showPLUExample()
Generate a matrix A and show P, L, and U as calculated
by computePLU(A). Also show the products PA and LU and 
show the squared error |PA - LU|^2
"""
    B = rand(N,N)
    A = transpose(B)*B+identity(N)
    P,L,U = computePLU(A)
    println("PA = ")
    display(P*A)
    println("LU = ")
    display(L*U)
    @printf("Error = %.4e",sum((P*A-L*U).^2)) 
end

function plotPLUTestResults(N,data,title,ylabel)
"""
    plotPLUTestResults()
Plot results of PLU algorithm testing (convenience function).
N = array of matrix sizes (scalars, not tuples; matrices are square)
data = vector of errors or execution times corresponding to entries of N 
"""
    xplot = log.(N)./log(10)
    yplot = log.(data)./log(10)
    p = scatter(xplot,yplot,
        xlabel = "log_10(N)", ylabel = ylabel,legend=false)    
    a0,a1 = linear_fit(xplot,yplot)
    plot!(p,xplot,a0 .+ a1 .* xplot)
    title!(@sprintf("%s\nfit slope = %.2f",title,a1))
    display(p)
end

function plotsolvePLUAccuracy(ntrials,maxpow10)
"""
    plotsolvePLUAccuracy()
Solve Ax = b, where A is an NxN matrix, using solvePLU for different
values of N and plot the error |Ax-b|^2 as a function of N 
"""
    # N-values are logarithmically spaced between 10 and 10^maxpow10
    N = Int.(floor.(10 .^range(1,maxpow10,ntrials)))
    errors  = zeros(ntrials)
    for (i,n) in enumerate(N)
        @printf("\n i = %d, n = %d",i,n)
        errors[i] = testsolvePLU(n)
    end
    plotPLUTestResults(N, errors, "PLU Solver Error","log_10(|Ax - b|^2)")
    #=
    p = scatter(log.(N),log.(errors),title = "PLU solver error",
        xlabel = "log(N)", ylabel = "log(|Ax - b|^2)")
    a0,a1 = linear_fit(log.(N),log.(errors))
    plot!(p,log.(N),a0 .+ a1 .* log.(N))
    title!(@sprintf("PLU Accuracy\nfit slope = %.2f",a1))
    display(p)
    =#
end

function plotPLUAccuracy(ntrials,maxpow10)
    """
        plotPLUAccuracy()
    Compute the squared matrix norm of 
    the difference PA - LU for the LU factorization of A at various
    N and plot the data 
    """
        # N-values are logarithmically spaced between 10 and 10^maxpow10
        N = Int.(floor.(10 .^range(1,maxpow10,ntrials)))
        errors  = zeros(ntrials)
        for (i,n) in enumerate(N)
            @printf("\n i = %d, n = %d",i,n)
            errors[i] = testPLU(n)
        end
        plotPLUTestResults(N, errors, "PLU Factorization Error","log_10(|PA-LU|^2)")
        #=
        p = scatter(log.(N),log.(errors),title = "PLU factorization error",
            xlabel = "log(N)", ylabel = "log(|PA-LU|^2)")
        a0,a1 = linear_fit(log.(N),log.(errors))
        plot!(p,log.(N),a0 .+ a1 .* log.(N))
        title!(@sprintf("PLU Accuracy\nfit slope = %.2f",a1))
        display(p)
        =#
    end

function plotPLUPerformance(ntrials,maxpow10)
"""
    plotPLUPerformance()
Plot the time required to solve Ax = b, where A is an NxN matrix,
at different values of N. 
"""  
    # N-values are logarithmically spaced between 10 and 10^maxpow10
    N = Int.(floor.(10 .^range(1,maxpow10,ntrials)))
    times  = zeros(ntrials)
    for (i,n) in enumerate(N)
        @printf("\n i = %d, n = %d",i,n)
        B = rand(n,n)
        A = transpose(B)*B + identity(n)
        b = rand(n,1)
        times[i] = @elapsed solvePLU(A,b)
    end
    # Plot the performance data 
    plotPLUTestResults(N, times, "PLU Performance","log_10(elapsed time (s))")
    #=
    p = scatter(log.(N),log.(times),
    xlabel="log(N)", ylabel="log(Elapsed time (s))")
    # Fit a line to the data
    a0,a1 = linear_fit(log.(N),log.(times))
    plot!(p,log.(N),a0 .+ a1 .* log.(N))
    title!(@sprintf("PLU Performance\nfit slope = %.2f",a1))
    display(p)
    =#
end

function norm(v)
"""
    norm(v)
Compute the norm of the vector (or array) v.
"""
    return sqrt(sum(v.^2))

end

function normsq(v)
"""
    compute the squared norm of v
"""
    return sum(v.^2)
end

function matconj(v,A)
"""
    compute the conjugated inner product
transpose(v)*A*v
"""
    return transpose(v)*A*v
end

function spdrmat(N)
    """
        pdtestmat(N)
    Return a symmetric positive-definite random NxN 
    matrix
    """
    B = rand(N,N)
    return B'*B + identity(N)
end

function conj_grad(A,x0,b,eps,maxiter)
"""
    conj_grad(A,x0,b,eps,maxiter)
Solve Ax = b using the conjugate gradient method with 
initial guess x0. 
eps: tolerance
maxiter: maximum number of iterations

Algorithm from Schewchuk, "An Introduction to the 
Conjugate Gradient Method Without the Agonizing
Pain, Edition 1 1/4," 1994
"""

    N = size(A)[1]
    iter = 0
    x = x0
    r = b - A*x0 # initial residual
    err =  norm(r)/norm(x0)
    d = r # initial step direction
    while err >= eps && iter <= maxiter
        alpha = normsq(r)/matconj(d,A)
        x = x + alpha .* d
        rnext = r - alpha.*A*d 
        beta = normsq(rnext)/normsq(r)
        d = rnext + beta.*d
        r = rnext
        err = norm(r)/norm(x)
        iter += 1
        #@printf("\niter = %d\n",iter)
        #@printf("err = %f\n",err)
    end
    return x 
end

# Create second derivative matrix
function ddMatrix(N;make_sparse=true)
    """
        ddMatrix(N)
    compute the discrete second derivative matrix as either a dense
    array or a sparse array
    N+1 = number of nodes in spatial discretization
    """
    if !make_sparse 
        # Create A in dense format
        A = zeros(N-1,N+1)
        for i = 1:N-1
            A[i,i:i+2] = [1 -2 1]
        end
    else   
        # Create A in sparse format 
        rowInds = [1:N-1;1:N-1;1:N-1]
        colInds = [1:N-1;2:N;3:N+1]
        vals = [ones(N-1); -2 .* ones(N-1); ones(N-1)]
        A = sparse(rowInds,colInds,vals) 
    end
    return A
end


"""
    symmetric_D_periodic(N,dx)
Compute the symmetric NxN first derivative matrix with 
periodic boundary conditions and node spacing dx 
"""
function symmetric_D_periodic(N,dx)
    row = zeros(N)
    row[1] = 1
    row[N-1] = -1
    D = zeros(N,N)
    for i = 1:N
        D[i,:] = circshift(row,i)
    end
    D = (1/(2*dx))*D
    return D 
end

"""
    upwind_D_Dirichlet(N,dx)
Compute the upwind first derivative for a problem with clamped
boundary conditions. N = (number of nodes)-1 = number of spatial 
intervals
"""
function upwind_D_Dirichlet(N,dx;makesparse=true)
    if makesparse
        D = sparse([1:N;1:N],[1:N;2:N+1],[-1 .*ones(N);1 .*ones(N)])
    else
        D = zeros(N,N+1)
        row = zeros(N+1)
        row[1] = -1
        row[2] = 1
        for i = 1:N
            D[i,:] = circshift(row,i-1)
        end
    end
    D = D./dx
    return D
end