#has to be moved inside a function



A=[]
λ_set = 1:10
flagConv = false
β_CD = zeros(numλ,N)
β_tilde = zeros(N)
NN=convert(Float64,N)
pos = 1



for k in 1:1#numλ #Lambda LOOP
    
    change = true
    for i in 1:I #Main outer LOOP: fix maximum
    #if flagConv == false
        
        if !change
            break
        end

        for c1 in IT.components
            for j in 1:IT.nelements[c1]

                temp = 0.0

                #GramMatrix LOOP (double loop: component and element)
                for c2 in IT.components
                    for l in activeSet[c2]
                        if l
                            temp = temp + GM[c1,c2](l,j, data)*β_tilde[c2][l]
                        end
                    end
                end

                β_ols = β_tilde[c1][j] + 1.0/NN * (xdy[c1][j] - temp)

                α = 0.0
                
                if abs(β_ols) <= α
                    
                    if activeSet[c1][j]
                        β_tilde[c1][j] = 0.0 #talvez nem precise
                        activeSet[c1][j] = false
                        change = true
                    end
                    
                else
                    
                    β_tilde[c1][j] = sign(β_ols)*(abs(β_ols) - α)
                    
                    if !activeSet[c1][j]
                        activeSet[c1][j] = true
                        change = true
                    end 
                    
                end
            end
        end
        change = true
    end

    
    for j in 1:(N-1)
        β_CD[j, pos] = β_tilde[j]
    end
    pos = pos + 1

end
                



