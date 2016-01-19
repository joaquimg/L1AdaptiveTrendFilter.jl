#has to be moved inside a function



A=[]
λ_set = 1:10
flagConv = false
β_CD = zeros(numλ,N)
β_tilde = zeros(N)
NN=convert(Float64,N)
pos = 1



for k in 1:numλ
    
    
    for i in 1:I #iter num
    #if flagConv == false
        
        for j in 1:NN#ℵ  #(N)#-1)
            
            temp = 0.0
            
            for l in A
                temp = temp + GM(l,j, L, U, N)*β_tilde[l]
            end
                    
            β_ols = β_tilde[j] + 1.0/NN * (xdy[j] - temp)
            
            α = 0.0
            
            if abs(β_ols) <= α
                
                β_tilde[j] = 0.0
                
                setdiff!(A,j)
                
                else
                
                β_tilde[j] = sign(b_ols)*(abs(β_ols) - α)
                
                union!(A,j)
                
            end
            
        end
        
    end
    
    for j in 1:(N-1)
        β_CD[j, pos] = β_tilde[j]
    end
    pos = pos + 1
end
                



