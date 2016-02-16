


# STEP
#correct here to remove U an L
function xdy_step(IT,y::Vector{Float64},d)
    N = IT.obs

    xdy = Vector{Float64}(N-1) 
    sumY = 0.0
    
    for i in 1:(N-1)
        sumY = sumY + y[i]
        xdy[i] = sumY*(-1.0/d.σ[STEP][i])
    end
    
    return xdy
end

#SLOPE
function xdy_slope(IT,y::Vector{Float64},d)
    N = IT.obs

    xdy = zeros(N) 

    y_s = sum(y)

    for i in 1:N

        for j in i:N
            xdy[i] += y[j] * (1.0+j-i)
        end

        xdy[i] =( xdy[i] -d.μ[SLOPE][i]*y_s )/d.σ[SLOPE][i]

    end

    return xdy
end

#SPIKE
function xdy_spike(IT,y::Vector{Float64},d)
    N = IT.obs

    xdy = zeros(N) 

    y_s = sum(y)

    for i in 1:N
        
        xdy[i] =( y[i] -d.μ[SPIKE][i]*y_s)/d.σ[SPIKE][i]

    end

    return xdy
end

#sin
function xdy_sin(IT,y::Vector{Float64},d)
    N = IT.obs
    
    nf = IT.nelements[4] 
    
    xdy0 = zeros(nf) 

    y_s = sum(y)

    for i in 1:nf

        for j in 1:N

            xdy0[i] += y[j] * sin(d.fs[i]*j)
        end
        
        xdy0[i] =( xdy0[i] -d.μ[SIN][i]*y_s)/d.σ[SIN][i]

    end

    return xdy0
end

#cos
function xdy_cos(IT,y::Vector{Float64},d)
    N = IT.obs

    nf = IT.nelements[5]

    xdy = zeros(N-1) 

    y_s = sum(y)

    for i in 1:nf

        for j in 1:N
            xdy[i] += y[j] * cos(d.fc[i]*j)
        end
        
        xdy[i] =( xdy[i] -d.μ[COS][i]*y_s)/d.σ[COS][i]

    end

    return xdy
end


xdy = Vector{Function}(5)

xdy[STEP] = xdy_step
xdy[SLOPE] = xdy_slope
xdy[SPIKE] = xdy_spike
xdy[SIN] = xdy_sin
xdy[COS] = xdy_cos
