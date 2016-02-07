const STEP = 1
const SLOPE = 2
const SPIKE = 3
const SIN = 4
const COS = 5



# STEP
#correct here to remove U an L
function xdy_step(IT,y::Vector{Float64},d)
    N = IT.obs

    xdy = Vector{Float64}(N-1) 
    sumY = 0.0
    
    for i in 1:(N-1)
        sumY = sumY + y[i]
        xdy[i] = sumY*(-1.0/d.σt[i])
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

        xdy[i] =( xdy[i] -d.μl[i]*y_s )/d.σl[i]

    end

    return xdy
end

#SPIKE
function xdy_spike(IT,y::Vector{Float64},d)
    N = IT.obs

    xdy = zeros(N) 

    y_s = sum(y)

    for i in 1:N
        
        xdy[i] =( y[i] -d.μp[i]*y_s)/d.σp[i]

    end

    return xdy
end

#sin
function xdy_sin(IT,y::Vector{Float64},d)
    N = IT.obs
    
    nf = IT.nelements[4] 
    
    xdy = zeros(N) 

    y_s = sum(y)

    for i in 1:nf

        for j in 1:N
            xdy[i] += y[j] * sin(d.fs[i]*j)
        end
        
        xdy[i] =( xdy[i] -d.μs[i]*y_s)/d.σs[i]

    end

    return xdy
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
        
        xdy[i] =( xdy[i] -d.μc[i]*y_s)/d.σc[i]

    end

    return xdy
end


