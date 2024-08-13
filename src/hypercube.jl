function generate_hypercube_points(n::Int)
    if n <= 0
        throw(ArgumentError("Dimension n must be a positive integer"))
    end
    
    paths = []
    
    for i in 0:2^n-1
        bin_str = string(i, base=2)
        bin_str = "0"^(n - length(bin_str)) * bin_str
        point = [parse(Int, c) for c in bin_str]
        push!(paths, point)
    end
    
    return paths
end

function rprojection(point)
    1-sum([(1-point[i])*2.0^(-i) for i in eachindex(point)])
end

function overlap(σ1,σ2)
    n = length(σ1)
    ov = n
    k =0
    for k in 1:n
        if σ1[k] == σ2[k]
            ov = ov - 1
        end
    end
    return ov/n
end
