"""
    upper_triangulate(A,b)

Put an augmented system in upper triangular form using elementary 
row operations (ie, Gaussian elimination)
Method currently doesn't support partial pivoting

"""
function upper_triangulate(A,b)
    N = size(A)[1]
    M = [A b]
    for j = 1:N-1
        eliminate!(M,j)
    end
    return M
end

function eliminate!(M,j)
    pivot = M[j,j] # Pivot element
    pivot_row = M[j,:]

    for k = j+1:N
        fac = M[k,j]/pivot
        M[k,:] = M[k,:]-fac*pivot_row 
    end
end


A = Matrix{Float64}(undef,4,4)
A[1,:] = [1 2 -1 1]
A[2,:] = [-1 1 2 -1]
A[3,:] = [2 -1 2 2]
A[4,:] = [1 1 -1 2]

b = [6; 3; 18; 8]

upper_triangulate(A,b)

#=
N = size(A)[1]
# Form augmented matrix 
M = [A b]

# Iterate over columns 1:N-1
for j = 1:N-1
    # Now working on column j
    # Eliminate all entries below a_jj using row j 
    pivot = M[j,j] # Pivot element
    pivot_row = M[j,:]

    for k = j+1:N
        fac = M[k,j]/pivot
        M[k,:] = M[k,:]-fac*pivot_row 
    end
end

display(M)
=#