
# first order autoregressive model
function AR1(prev::Vector{Float64},θ)
    tau_prev = prev[2]
    birth_prev = prev[1]
    tau = max(0.00001,tau_prev*θ.a + (1-θ.a) + sqrt(1-θ.a^2)rand(Normal(0,θ.σ)))
    time = tau_prev + birth_prev
    return [time,tau]
end


# cell size control model
function CSC(prev::Vector{Float64},θ)
    birth_prev = prev[1]
    size_prev = prev[2]
    growth_prev = prev[3]
    tau_prev = prev[4]
    growth = max(0.00001,growth_prev*θ.ag + (1-θ.ag) + sqrt(1-θ.ag^2)rand(Normal(0,θ.σg)))
    size = max(0.00001,size_prev*θ.aM + (1-θ.aM) + sqrt(1-θ.aM^2)rand(Normal(0,θ.σM)))
    tau = 1/growth * log(2*size/size_prev)
    time = tau_prev + birth_prev
    return [time,size,growth,tau]
end



function fitar(y)
    L = length(y[:,1])
    X = hcat(ones(L-1),y[1:L-1,:])
    b = X\y[2:end,:]
    v = mean((y[2:end,:] .- X*b).^2,dims=1)'
    return θ = (b = Matrix(b),v = Matrix(v))
end

function runar(y0,θ,n)
    n_vars = length(y0)
    yp = zeros((n,n_vars))
    yp[1,:] = y0
    for k in 2:n
        yp[k,:] = θ.b' * vcat(1,yp[k-1,:]) .+ rand(MvNormal(zeros(n_vars),diagm(θ.v[:,1])))
    end
    return yp
end