const STEP = 1
const SLOPE = 2
const SPIKE = 3
const SIN = 4
const COS = 5



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
function xdy_slope(IT,y::Vector{Float64},d)
    N = IT.obs

    xdy = zeros(N) 

    y_s = sum(y)

    for i in 1:N

        for j in i:N
            xdy[i] += y[j] * (1+j-i)
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
            xdy[i] += y[j] * sin(d.f[i]*j)
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
            xdy[i] += y[j] * cos(d.f[i]*j)
        end
        
        xdy[i] =( xdy[i] -d.μc[i]*y_s)/d.σc[i]

    end

    return xdy
end


