
function l1_adaptive_trend_filter(
    y, components; f = Vector{Float64}(0), numλ=40, MAXITER=100,
    lower_bounds=-10e+7*ones(TOTALCOMPONENTS), upper_bounds=10e+7*ones(TOTALCOMPONENTS)
    )

  # subtracting the mean
  y_mean = mean(y)
  y = y - y_mean

  #prepare check for dimension sizes
  const N = size(y)[1]
  # iterator
  const IT = initIT_range(N,components,f,MAXITER=MAXITER)
  # ?
  const d = initData(IT,f, f)
  # inner products between response y and components
  const xdy = initXDY(IT,y,d)
  # build regularization path
  const Λ = computeλvec(IT,xdy,numλ)

  @time @fastmath BCD,β1,β2, y_best = coordinate_descent(
    IT, d, xdy, Λ, y, lower_bounds, upper_bounds
    )
  y_best = y_best + y_mean
  return BCD,β1,β2, y_best
end

# coordinate descent algorithm for the regularization path Λ
function coordinate_descent(IT::iterator, d::dataCD, xdy, Λ::Vector{Float64}, y, lower_bounds, upper_bounds; sparse=0)

  # initializations
  if sparse == 1
    BCD, β_tilde, β, activeSet = initSparse(IT)
  else
    @time BCD, β_tilde, β, activeSet = initDense(IT)
  end

  # memory allocation
  BIC = Inf::Float64
  β_ols = 0.0::Float64
  partial_fit = 0.0::Float64
  β_best_unbiased = 0
  β_best_biased = 0
  y_best = 0

  # regularization path
  c = 0
  @inbounds for λ in Λ
  c+=1
  println(c)
    change = true
    # loop until active set converges
    @inbounds for iter in 1:IT.maxIter
    println(iter)
      if !change
        break
      end

      change = false
      # cycle through every component
      @inbounds for c1 in IT.components
        @inbounds for j in IT.elements[c1]

          # compute the partial fit with the components in the active set
          partial_fit = 0.0
          @inbounds for c2 in IT.components
            @inbounds for l in 1:size(activeSet[c2])[1]
              if activeSet[c2][l] 
               @inbounds partial_fit += GM2(c1,c2,j,l,d,IT) * β_tilde[c2][l] 
              end
            end
          end

          # univariate ordinary leasts squares coefficient
          β_ols =  β_tilde[c1][j] + (1.0/IT.obs) *(xdy[c1][j] - partial_fit)

          # soft thresholding operator
          if abs(β_ols) <= λ
            if activeSet[c1][j]
              β_tilde[c1][j] = 0.0
              activeSet[c1][j] = false
              change = true
            end
          else
            β_tilde[c1][j] = sign(β_ols) * (abs(β_ols) - λ)

            # projection onto the box constraints [lower_bound, upper_bound]
            #β_tilde[c1][j] = max(β_tilde[c1][j], lower_bounds[c1])
            #β_tilde[c1][j] = min(β_tilde[c1][j], upper_bounds[c1])

            if !activeSet[c1][j] && β_tilde[c1][j] != 0.0
              activeSet[c1][j] = true
              change = true
            end
          end

        end
      end
    end

#    # bayesian information criterion
#    push!(BCD,deepcopy(β_tilde))
#    β_unbiased = compute_OLS(β_tilde, λ, activeSet, IT, xdy, d, lower_bounds, upper_bounds)
#    BIC_new, y_hat= compute_BIC(y, β_unbiased, IT, d)
#
#    # save the best fit so far
#    if BIC_new < BIC
#       BIC = BIC_new
#       y_best = copy(y_hat)
#       β_best_unbiased = deepcopy(β_unbiased)
#       β_best_biased = deepcopy(β_tilde)
#    end

  end
  return BCD, β_best_unbiased, β_best_biased, y_best
end
