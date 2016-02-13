

#lambda max
function build_xdy0(N::Int,y::Vector{Float64}) 
    xdy0 = 0.0
    for i in 1.0:N
            xdy0 = xdy0 + i*y[i]
    end
    return xdy0
end


function findλmax(N::Int,λmax=0.0::Int,sigma0)
    
    #λmax   = 0.0
    temp  = 0.0
    
    maxPos = 0
    
    for i in 1:(N-1)
        
        temp = abs(sigma[i]*xdy11[i])
        
        if temp > λmax
            
            λmax = temp
            maxPos = i
            
        end
        
    end
    if maxPos > 0
        λmax = λmax/sigma[maxPos]
    else
        λmax = λmax/sigma0
    end
    
    return λmax
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

#achar um mx por coponente depois juntar tudo...

