include("../my_solvers.jl")

N = 4
ntrials = 10
maxpow10 = 3

showPLUExample(N)
plotPLUAccuracy(ntrials,maxpow10)
plotsolvePLUAccuracy(ntrials,maxpow10)
plotPLUPerformance(ntrials,maxpow10)