#has to be moved inside a function


function CD(y,components; f = Vector{Float64}(0))

    ym = mean(y)

    y = y - ym

    #prepare check for dimension sizes
    N = size(y)[1]

    IT = initIT_range(N,components,f)

    d = initData(IT,f, f)

    xdy = initXDY(IT,y,d)

    λ = computeλvec(IT,xdy,10)

    BCD = CoordinateDescent(IT,d,xdy,λ)

    return BCD
end

# coordinate descent for the whole regularization path Λ
function CoordinateDescent(IT, d, xdy, Λ; sparse=0)

  if sparse == 1
    BCD, β_tilde, β, activeSet = initSparse(IT)
  else
    BCD, β_tilde, β, activeSet = initDense(IT)
  end

  BIC = Inf::Float64
  β_ols = 0.0 :: Float64
  partial_fit = 0.0 :: Float64

  # regularization path
  for λ in Λ

    change = true
    # loop until active set converges
    for iter in 1:IT.maxIter
      if !change
        break
      end

      # cycle through every component
      for c1 in IT.components
        for j in IT.elements[c1]

          # compute the partial fit with the components in the active set
          partial_fit = 0.0
          for c2 in IT.components
            for l in 1:size(activeSet[c2])[1]
              if activeSet[c2][l] && (c1, j) != (c2, l)
                partial_fit = partial_fit + GM[c1,c2](j,l,d,IT) * β_tilde[c2][l]
              end
            end
          end

          # univariate ordinary leasts squares coefficient
          β_ols = β_tilde[c1][j] + (1.0/IT.obs) * (xdy[c1][j] - partial_fit)

          # soft thresholding operator
          if abs(β_ols) <= λ
            if activeSet[c1][j]
              β_tilde[c1][j] = 0.0 #talvez nem precise
              activeSet[c1][j] = false
              change = true
            end
          else
            β_tilde[c1][j] = sign(β_ols) * (abs(β_ols) - λ)

            if !activeSet[c1][j]
              activeSet[c1][j] = true
              change = true
            end
          end
        end
      end
    end

    # bayesian information criterion
    βtemp = copy(β_tilde)
    println(βtemp)
    #push!(BCD,copy(β_tilde))
    #println(BCD)
    BIC_new = compute_BIC(y, βtemp, IT, d)
    if BIC_new < BIC
      BIC = BIC_new
      β = βtemp
    end
  end

  return BCD,β
end
