

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

#achar um mx por coponente depois juntar tudo...

