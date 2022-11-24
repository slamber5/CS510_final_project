using Plots 
include("dense_sparse_lu.jl")

# Generate plots 
functs = [timeLU,timesolve] 
titles = ["LU factorization performance","Backslash solve performance"]
for (i,funct) in enumerate(functs)
    Ns = [2;Integer.(ceil.(10 .^collect(1:0.2:3.6)))]
    display(Ns)
    Ns_sparse,tsparse = funct(Ns,nruns=3,make_sparse=true)
    Ns_dense,tdense = funct(Ns,nruns=3,make_sparse=false)
    p=plot(legend=:bottomright,title=titles[i],
    xlabel="Matrix size (N)", ylabel="Computation time (s)")
    plot!(p, Ns_sparse[2:end],tsparse[2:end,:],xaxis=:log,yaxis=:log,marker=:circle,
        color = :1,labels=["sparse" "" ""])
    plot!(p, Ns_dense[2:end],tdense[2:end,:],xaxis=:log,yaxis=:log,marker=:circle,
        color = :2,labels=["dense" "" ""])
    display(p)
end

