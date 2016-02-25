#has to be moved inside a function


function CD(y,components; f = Vector{Float64}(0), numλ = 40)

    ym = mean(y)

    y = y - ym

    #prepare check for dimension sizes
    const N = size(y)[1]

    const IT = initIT_range(N,components,f,MAXITER=20)

    const d = initData(IT,f, f)

    const xdy = initXDY(IT,y,d)

    const Λ = computeλvec(IT,xdy,numλ)

    @time @fastmath BCD,β1,β2, y_best = CoordinateDescent(IT,d,xdy,Λ,y)
#CoordinateDescent(IT,d,xdy,Λ,y)
    return BCD,β1,β2, y_best
end

# coordinate descent for the whole regularization path Λ
function CoordinateDescent(IT, d, xdy, Λ, y; sparse=0)

  if sparse == 1
    BCD, β_tilde, β, activeSet = initSparse(IT)
  else
    @time BCD, β_tilde, β, activeSet = initDense(IT)
  end
  BIC = Inf::Float64
  β_ols = 0.0::Float64
  partial_fit = 0.0::Float64

  β1 = 0
  β2 = 0
  y_best = 0
#@bp
  # regularization path
  @inbounds for λ in Λ

    change = true
    # loop until active set converges
    for iter in 1:IT.maxIter
      if !change
        break
      end

      change = false
      # cycle through every component
      for c1 in IT.components
        for j in IT.elements[c1]

          # compute the partial fit with the components in the active set
          partial_fit = 0.0
          for c2 in IT.components
            @inbounds for l in 1:size(activeSet[c2])[1]
              @inbounds activeSet[c2][l] ? partial_fit += GM2(c1,c2,j,l,d,IT) * β_tilde[c2][l] *activeSet[c2][l] : true
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
     #println(βtemp)
     push!(BCD,deepcopy(β_tilde))
     #println(BCD)
     β_unbiased = compute_OLS(β_tilde,λ,activeSet,IT,xdy,d)
     BIC_new, y_hat= compute_BIC(y, β_unbiased, IT, d)
     if BIC_new < BIC

       BIC = BIC_new
       y_best = copy(y_hat)
       β1 = deepcopy(β_unbiased)
       β2 = deepcopy(β_tilde)
     end

  end
  return BCD,β1,β2,y_best
end
