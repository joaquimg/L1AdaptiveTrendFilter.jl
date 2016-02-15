



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

function computeλvec(IT,xdy,numλ; logarit = true )

    λmax = findλmax(IT,xdy)


    if logarit 

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

