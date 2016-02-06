function compute_level(y::Vector{Float64}, β::Vector{Float64}, N::Int)
    
    for j in 1:(N-1)
        for i in (j+1):N
            y[i] = y[i] + β[j]
        end
    end
    
    return y
end

function compute_slope(y::Vector{Float64}, β::Vector{Float64}, N::Int)
    
    for j in 1:(N-1)
        for i in (j+1):N
            y[i] = y[i] + (i-j)*β[j]
        end
    end
    
    return y
end

function compute_spike(y::Vector{Float64}, β::Vector{Float64}, N::Int)

    for i in 1:N
        y[i] = y[i] + β[i]
    end

    return y
end

function compute_sin()

    for i in 1:N
        for f in 1:nfreqs
            y[i] = y[i]+sin(freq[f]*i)*β[f]
        end
    end
    return y
end

function compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β::Vector{Float64}, N::Int, ɛ = 1e-5::Float64)
    
    BIC = 0.0
    err = zeros(N)
    
    N2 = size(β,1)
    
    for i in 1:N
        err[i] = y[i] - y_hat[i]
    end
    
    k = 0 ::Int
    
    for i in 1:N2
        if abs(β[i]) > ɛ
            k = k + 1
        end
    end
    
    BIC = N2*var(err) + k * ln(N2)
    
    return BIC
    
end
    
    
    
    