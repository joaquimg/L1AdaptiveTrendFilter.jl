# Computes the inner products between each component and the given response

# STEP component
#correct here to remove U an L
function xdy_step(IT::iterator, y::Vector{Float64}, d::dataCD)
    N = IT.obs
    Nel = IT.nelements[STEP]
    xdy = zeros(Nel)
    y_s = sum(y)

    @inbounds for i in 1:Nel
        for j in i+1:N
            xdy[i] += y[j]
        end
        xdy[i]=(xdy[i] -d.μ[STEP][i]*y_s)
    end

    return xdy
end

# SLOPE component
function xdy_slope(IT::iterator, y::Vector{Float64}, d::dataCD)
    N = IT.obs
    Nel = IT.nelements[SLOPE]
    xdy = zeros(Nel)
    y_s = sum(y)

    @inbounds for i in 1:Nel
        for j in i:N
            xdy[i] += y[j] * (j-i)
        end
        xdy[i] =(xdy[i] -d.μ[SLOPE][i]*y_s)
    end
    return xdy
end

# SPIKE component
function xdy_spike(IT::iterator, y::Vector{Float64}, d::dataCD)
    N = IT.obs
    Nel = IT.nelements[SPIKE]
    xdy = zeros(Nel)
    y_s = sum(y)

    @inbounds for i in 1:Nel
        xdy[i] =(y[i] -d.μ[SPIKE][i]*y_s)
    end
    return xdy
end

# SINE component
function xdy_sin(IT::iterator, y::Vector{Float64}, d::dataCD)
    N = IT.obs
    nf = IT.nelements[SIN]
    xdy0 = zeros(nf)
    y_s = sum(y)

    @inbounds for i in 1:nf
        for j in 1:N
            xdy0[i] += y[j] * sin(d.fs[i]*j)
        end
        xdy0[i] =(xdy0[i] -d.μ[SIN][i]*y_s)
    end
    return xdy0
end

# COSINE component
function xdy_cos(IT::iterator, y::Vector{Float64}, d::dataCD)
    N = IT.obs
    nf = IT.nelements[COS]
    xdy = zeros(nf)
    y_s = sum(y)

    @inbounds for i in 1:nf
        for j in 1:N
            xdy[i] += y[j] * cos(d.fc[i]*j)
        end
        xdy[i] = (xdy[i]-d.μ[COS][i]*y_s)
    end
    return xdy
end

xdy = Vector{Function}(TOTALCOMPONENTS)
xdy[STEP] = xdy_step
xdy[SLOPE] = xdy_slope
xdy[SPIKE] = xdy_spike
xdy[SIN] = xdy_sin
xdy[COS] = xdy_cos
