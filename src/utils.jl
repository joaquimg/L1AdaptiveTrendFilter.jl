
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
    else

        vec = λmax*(linspace(1,0,numλ))

    end
    return vec
end



function compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β, IT; ɛ = 1e-5::Float64)
    
    BIC = 0.0
    err = zeros(y)
    
    N2 = IT.ttelements
    
    for i in 1:IT.obs
        err[i] = y[i] - y_hat[i]
    end
    
    k = 0.0::Float64
    
    for i in IT.components
        for j in IT.elements[i]
            if abs(β[i][j]) > ɛ
                k += 1.0
            end
        end
    end
    
    BIC = N2 * var(err) + k * log(N2)
    
    return BIC
    
end
function compute_BIC(y::Vector{Float64}, β, IT, d; ɛ = 1e-5::Float64, std = 1)

    y_hat, β_new = compute_estimate(zeros(y), IT, β, d)
    
    if std == 1
        BIC = compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β_new, IT)
    else
        BIC = compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β, IT)
    end

    return BIC
end

