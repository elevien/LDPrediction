
# ----------------------------------------------------------------------
# estimators

function FDE(T::Vector{Float64},n::Int64)
    μ = mean(T)
    dT = T .- μ
    f(x) = log(mean(exp.(-dT .* x)))/n .- x .* μ ./n .+log(2)
    sol = nlsolve(f, [0.69])
    return sol.zero[1]
end


function FTE(N::Vector{Int64},t::Float64)
    μ = mean(N)
    dN = N .- μ
    return log(mean(exp.(dN .* log(2))))/t .+ μ*log(2)/t
end

# ----------------------------------------------------------------------
# # estimators

function scalingdata_FDE(df::DataFrame,sizes::Vector{Int64})

    αs = []
    FDEs = []
    Ms = []
    for g in sizes
        dfg = grouplineages(df,g)
        T = combine(groupby(dfg,:group),:t => last => :tf).tf
        push!(FDEs,FDE(T,g))
        a = sqrt(g ./ length(unique(dfg.group)))
        push!(αs,a)
        push!(Ms,length(unique(dfg.group)))
    end
    DataFrame(:α => Vector{Float64}(αs),:M => Ms,:n => sizes,:fde => Vector{Float64}(FDEs))
end



# ------------------------------------------------------------------------------------
# auxillary functions

function grouplineages(X::DataFrame,group_size::Int64)
    dfs = []
    g = 0
    for l in unique(X.lineage)

        Xl = X[X.lineage .== l,:]
        k = 1
        while k+group_size<length(Xl.t)
            df = Xl[k:(k+group_size),:]
            df[:,:group] = ones(length(df.t))*g
            df[:,:t] = df[:,:t] .- df[:,:t][1]
            df[:,:n] = df[:,:n] .- df[:,:n][1]
            push!(dfs,df)
            k+=group_size+1
            g+=1
        end
    end
    df=  vcat(dfs...)
end


# ----------------------------------------------------------------------
# these are used to make data for fig3, but eventually hope to replace 
# them with functions above in data analysis pipeline


function Lhat_FDE(θ)
    n = trunc(Int,θ.α^2*log(θ.M)/(2*log(2)))
    T = rand(Normal(0,sqrt(n)*θ.σ),θ.M)
    f(x) = log(mean(exp.(-T .* x)))/n .- x .*θ.τ0 .+log(2)
    sol = nlsolve(f, [0.69])
    return sol.zero[1]
end


function Lhat_FTE(θ)
    t = trunc(Int,θ.α^2*log(θ.M)/(2*log(2)))
    N = zeros(θ.M)
    for k in 1:θ.M
        s = 0
        n = 0
        while s<t
            s = s+ rand(Normal(θ.τ0,θ.σ))
            n = n+1
        end
        N[k] = n
    end     
    dN =  N .- mean(N)
    return log(mean(exp.(dN .* log(2))))/t .+ mean(N)*log(2)/t
end
