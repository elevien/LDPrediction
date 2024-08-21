
βc = 2sqrt(log(2))

function fed_rem_ht(β)
    -β/4 - log(2)/β
end


function fed_rem_lt(β)
    -sqrt(log(2))
end

function fed_rem(β)
   β<βc ? fed_rem_ht(β) : fed_rem_lt(β)
end



function f_ht(z,θ)
    β = z* θ.σ * θ.α
    (2/θ.α^2)*(-β*fed_rem_ht(β)-log(2)) - z*θ.τ0 + log(2)
end

function f_lt(z,θ)
    β = z*θ.σ*θ.α
    (2/θ.α^2)*(-β*fed_rem_lt(β)-log(2)) - z*θ.τ0 + log(2)
end

function Lhat_rem_theory(θ)
    zc = βc/(θ.σ*θ.α)
    f = z -> z > zc ? f_lt(z,θ)  : f_ht(z,θ)
    F(z) = f.(z)
    sol = nlsolve(F, [1.0],ftol=10e-8,xtol=10e-20,iterations=10^8)
    max(sol.zero...)
end

function Lhat_rem_theory2(θ)
    σ2y = θ.σ^2/(1-θ.c)^2
    b = log(2)*θ.α*sqrt(2*σ2y)
    return log(2)/θ.τ0 + 1/θ.α^2*(-b*fed_rem(b) - log(2))
end