
# STEP
function xdy_step(y::Vector{Float64},L::Vector{Float64},U::Vector{Float64},N::Int)
    xdy = Vector{Float64}(N-1) 
    sumY = 0.0
    
    for i in 1:(N-1)
        sumY = sumY + y[i]
        xdy[i] = sumY*(U[i]-L[i])
    end
    
    return xdy
end

#SLOPE
function xdy_slope(y::Vector{Float64},N::Int,...)
    xdy = zeros(N-1) 

    y_s = sum(y)

    for i in 1:N

        for j in i:N
            xdy[i] += y[j] * (1+j-i)
        end
        
        xdy[i] =( xdy[i] -mu[i]*y_s)/sigma[i]

    end

    return xdy
end

#SPIKE
function xdy_spike(y::Vector{Float64},N::Int,...)
    xdy = Vector{Float64}(N-1) 

    y_s = sum(y)

    for i in 1:N
        
        xdy[i] =( y[i] -mu[i]*y_s)/sigma[i]

    end

    return xdy33
end

#sin
function xdy_sin(y::Vector{Float64},N::Int,...)
    xdy = zeros(N-1) 

    y_s = sum(y)

    for i in 1:N#nfreqs

        for j in 1:N
            xdy[i] += y[j] * sin(f[i]*j)
        end
        
        xdy[i] =( xdy[i] -mu[i]*y_s)/sigma[i]

    end

    return xdy
end

#cos
function xdy_sin(y::Vector{Float64},N::Int,...)
    xdy = zeros(N-1) 

    y_s = sum(y)

    for i in 1:N

        for j in 1:N
            xdy[i] += y[j] * cos(f[i]*j)
        end
        
        xdy[i] =( xdy[i] -mu[i]*y_s)/sigma[i]

    end

    return xdy
end


