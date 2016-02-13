
function compute_estimate()

end

function compute_step(y, β, μ, σ, IT)
    
    for j in 1:(N-1)
        for i in (j+1):N
            y[i] = y[i] + β[j]
        end
    end
    
    return y
end

function compute_slope(y::Vector{Float64}, β::Vector{Float64}, IT)
    
    for j in 1:(N-1)
        for i in (j+1):N
            y[i] = y[i] + (i-j)*β[j]
        end
    end
    
    return y
end

function compute_spike(y::Vector{Float64}, β::Vector{Float64}, IT)

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


    
    
    
    