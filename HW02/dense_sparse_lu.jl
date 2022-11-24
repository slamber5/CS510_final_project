using LinearAlgebra
using SparseArrays
using Plots

include("../my_solvers.jl")

# Create A in dense or sparse format
function ddMatrix(N;make_sparse=true)
    """
        ddMatrix(N)
    compute the discrete second derivative matrix as either a dense
    array or a sparse array.
    """
    if !make_sparse 
        # Create A in dense format
        A = zeros(N,N)
        A[1,1:2] = [-2 1]
        for i = 2:N-1
            A[i,i-1:i+1] = [1 -2 1]
        end
        A[N,N-1:N] = [1 -2]
    else   
        # Create A in sparse format 
        rowInds = [1:N-1;1:N;2:N]
        colInds = [2:N;1:N;1:N-1]
        vals = [ones(N-1); -2 .* ones(N); ones(N-1)]
        A = sparse(rowInds,colInds,vals) 
    end
    return A
end

# Compute and time LU factorization with LinearAlgebra
function timeLU(Ns; make_sparse=true,nruns=1)
    """
        timeLU(Ns)
    Time how long it takes Julia to compute the LU factorization of a 
    second derivative matrix as defined in ddMatrix. 
    Ns = vector of matrix sizes 
    nruns = number of trials for each N value
    make_sparse: whether or not to use a sparse array
    solve: whether or not to time 
    returns: t, (length(Ns) x nruns) array of factorization times 
    """
    
    t = zeros(length(Ns),nruns)
    for i in 1:nruns 
        for (j,N) in enumerate(Ns)
            A = ddMatrix(N,make_sparse=make_sparse)
            b = ones(N)
            t[j,i] = @elapsed lu(A)
        end
    end
    return Ns,t
end

# Solve Ax = b and time 
function timesolve(Ns; make_sparse=true,nruns=1)
    """
        timesolve(Ns)
    Time how long it takes Julia to solve Ax = b with LU factorization
    when A is a second derivative matrix and b is a vector of ones. 
    Ns = vector of matrix sizes 
    nruns = number of trials for each N value
    make_sparse: whether or not to use a sparse array
    solve: whether or not to time 
    returns: t, (length(Ns) x nruns) array of factorization times 
    """
    
    t = zeros(length(Ns),nruns)
    for i in 1:nruns 
        for (j,N) in enumerate(Ns)
            A = ddMatrix(N,make_sparse=make_sparse)
            b = ones(N)
            F = lu(A)
            t[j,i] = @elapsed F\b
        end
    end
    return Ns,t
end
