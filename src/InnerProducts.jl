

function load_xdy11(y::Vector{Float64},L::Vector{Float64},U::Vector{Float64},N::Int)
    xdy11 = Vector{Float64}(N-1) 
    sumY = 0.0
    
    for i in 1:(N-1)
        sumY = sumY + y[i]
        xdy11[i] = sumY*(U[i]-L[i])
    end
    
    return xdy11
end

function load_xdy22(y::Vector{Float64},N::Int,...)
    xdy22 = Vector{Float64}(N-1) 

    return xdy22
end
function load_xdy33(y::Vector{Float64},N::Int,...)
    xdy33 = Vector{Float64}(N-1) 

    return xdy33
end
function load_xdy44(y::Vector{Float64},N::Int,...)
    xdy44 = Vector{Float64}(N-1) 

    return xdy44
end