using Distributions
using StatsBase
using PythonPlot
using DataFrames
using Optim
using NLsolve
using Setfield
using CSV
using LinearAlgebra
using ProgressBars
include("estimators.jl")
include("rem_theory.jl")
include("models.jl")

folder = "./output/fig2"
mkpath(folder)

FIG_PATH = "/Users/elevien/Dropbox (Dartmouth College)/Apps/Overleaf/Finite Lineages Plos Journal/paper/figures"



c_range = collect(-0.1:0.01:0.1)



m_reps =5
L = 1000
M = 5000

for k in 1:length(c_range)
    c = c_range[k]
    println("running c = "*string(c))
    ar1s = []
    θ = (b = [1-c; c;;], v = [0.1*(1-c^2);;]);
    for i in ProgressBar(1:m_reps)
        times = runar([0],θ,M*L)
        n = vcat([collect(1:L) for k in 1:M]...) .-1;
        lineages = vcat([ones(L)*k for k in 1:M]...);
        ar1 = DataFrame((lineage=lineages,gt = times[:,1],t=cumsum( times[:,1]),n=n));
        ar1[:,:rep] = i*ones(length(times))
        ar1.t = vcat([ar1[ar1.lineage .==l,:t] .- ar1[ar1.lineage .==l,:t][1] for l in unique(ar1.lineage)]...);
        push!(ar1s,ar1)
    end
    CSV.write(folder*"/ar1"*string(k)*".csv",vcat(ar1s...))
end