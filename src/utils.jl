

#lambda max
function build_xdy0(N::Int,y::Vector{Float64}) 
    xdy0 = 0.0
    for i in 1.0:N
            xdy0 = xdy0 + i*y[i]
    end
    return xdy0
end


function findλmax(IT,xdy)
    
    λmax   = 0.0
    temp  = 0.0
    
    maxPos = 0
    
    for i in IT.components
        
        temp = maxabs(xdy[i])
        
        if temp > λmax
            
            λmax = temp

        end
        
    end
    
    return λmax
end

function computeλvec(IT,xdy,numλ)

    λmax = findλmax(IT,xdy)

    logarit = 1
    if logarit == 1

        vec =  λmax*(logspace(1,0,numλ)-1.0)/9.0 

    end
    return vec
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

