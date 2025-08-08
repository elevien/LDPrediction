# -------------------------------------------------------------
# setup 

using Distributions
using StatsBase
using PythonPlot
using DataFrames
using Optim
using NLsolve
using Setfield
using CSV

include("estimators.jl")


# save to output folder
cd(dirname(@__FILE__))
folder = "./output"
mkpath(folder)

# -------------------------------------------------------------
# First do FTE

# setup 
αrng = exp.(collect(log.(1):0.4:log.(700.0)));
Mrng = collect(100:30:1000);
n_samples = 100

# compute FDE estimates
println("running FDE ...")
Lhats_FDE = []
for m in Mrng
    θ = (n = 5,M = m,t = 10, τ0=1,σ=0.1,α=0.03)
    @time Y = hcat([[Lhat_FDE(@set θ.α = a) for k in 1:n_samples] for a in αrng]...); 
    push!(Lhats_FDE,Y)
end


df_FDE = vcat([DataFrame(hcat(Lhats_FDE[k],Mrng[k]*ones(n_samples)),:auto) for k in 1:length(Lhats_FDE)]...);
CSV.write(folder*"/df_FDE",df_FDE)
CSV.write(folder*"/alpha_FDE",DataFrame(alpha = αrng))


# -------------------------------------------------------------
# now do FTE

αrng_FTE = exp.(collect(log.(1):0.6:log.(400.0)));
Mrng_FTE = collect(100:50:1000);
Lhats_FTE = []
for m in Mrng_FTE
    θ = (n = 5,M = m,t = 10, τ0=1,σ=0.1,α=0.03)
    @time Y = hcat([[Lhat_FTE(@set θ.α = a) for k in 1:n_samples] for a in αrng_FTE]...); 
    push!(Lhats_FTE,Y)
end

df_FTE = vcat([DataFrame(hcat(Lhats_FTE[k],Mrng_FTE[k]*ones(n_samples)),:auto) for k in 1:length(Lhats_FTE)]...);
CSV.write(folder*"/df_FTE",df_FTE)
CSV.write(folder*"/alpha_FTE",DataFrame(alpha = αrng_FTE))