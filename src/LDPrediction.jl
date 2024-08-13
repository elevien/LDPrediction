module LDPrediction

    include("estimators.jl")
    include("rem_theory.jl")
    export Lhat_FDE,Lhat_FTE,fed_rem_ht,fed_rem_lt,βc
end