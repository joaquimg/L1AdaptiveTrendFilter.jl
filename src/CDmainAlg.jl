#has to be moved inside a function


function CD(y,components; f=Float64[])

    #prepare check for dimension sizes
    N = size(y)[1]

    IT = initIT_range(N,components,f)

    d = initData(IT,f, f)

    xdy = initXDY(IT,y,d)

    λ = computeλvec(IT,xdy,10)

    BCD = CoordinateDescent(IT,d,xdy,λ)

    return BCD
end


function CoordinateDescent(IT,d,xdy,λ; sparse = 1)

    if sparse == 1
        BCD,β_tilde,β,activeSet = initSparse(IT)
    else
        BCD,β_tilde,β,activeSet = initDense(IT)
    end

    BIC = Inf::Float64

    β_ols = 0.0 :: Float64
    temp = 0.0 :: Float64

    for k in λ#numλ #Lambda LOOP

        change = true
        for i in 1:IT.maxIter #Main outer LOOP: fix maximum
        #if flagConv == false

            if !change
                break
            end

            for c1 in IT.components
                for j in IT.elements[c1]

                    temp = 0.0

                    #GramMatrix LOOP (double loop: component and element)
                    for c2 in IT.components
                        for l in activeSet[c2]
                            if l
                                temp = temp + GM[c1,c2](l,j, data)*β_tilde[c2][l]
                            end
                        end
                    end

                    β_ols = β_tilde[c1][j] + 1.0/IT.ttelements * (xdy[c1][j] - temp)

                    #α = 0.0

                    if abs(β_ols) <= k#α

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

        push!(BCD,β_tilde)

        BIC_new = compute_BIC(y, β_tilde, IT)
        if BIC_new < BIC
            BIC = BIC_new
            β = β_tilde
        end

    end

    return BCD,β
end




