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