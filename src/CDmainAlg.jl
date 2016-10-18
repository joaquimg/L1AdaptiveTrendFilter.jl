
function l1_adaptive_trend_filter{T<:Real, S<:Integer}(
    y::Vector{T}, components::Vector{S}; f = Vector{T}(0), numλ::T=T(40), numγ::T=T(10), MAXITER::Integer=500, verbose::Bool=false,
    lower_bounds=-10e+7*ones(T, TOTALCOMPONENTS), upper_bounds=10e+7*ones(T, TOTALCOMPONENTS),
    )

  # subtracting the mean
  y_mean = mean(y)
  y = y - y_mean

  # prepare check for dimension sizes
  const N = size(y)[1]
  # iterator
  IT = initIT_range(N, components, f, MAXITER=MAXITER)
  # compute moments
  const d = initData(IT, f, f)
  # inner products between response y and components
  const xdy = initXDY(IT, y, d)
  # build regularization path for λ
  path, λ_max = compute_λ_path(IT, xdy, numλ, d)
  # build regularization path for γ
  const Γ = compute_γ_path(IT, xdy, numγ, d)

  # initialize weights
  w = Vector{Float64}[]
  for i in 1:TOTALCOMPONENTS
    if in(i,IT.components)
      push!(w, ones(IT.nelements[i]))
    else
      push!(w, ones(0))
    end
  end

  # lasso pass
  @fastmath β_path, y_path, β_best, y_best, λ_best, γ_best = coordinate_descent(
    IT, d, xdy, path, λ_max, [1.0,] , y, lower_bounds, upper_bounds, w ,verbose
    )

  # exclude the components the lasso has set to zero
  update_components!(IT, w, β_best, d, xdy)

  @fastmath β_path, y_path, β_best, y_best, λ_best, γ_best = coordinate_descent(
    IT, d, xdy, path, λ_max, Γ, y, lower_bounds, upper_bounds, w, verbose
   )

  # adding back the mean
  y_best = y_best + y_mean

  if verbose
    print(string(
            "best regularization according to BIC was (λ=", round(λ_best, 3),
            ", γ=", round(γ_best, 3), ") \n"
            ))
  end

  return β_path, y_path, β_best, y_best, λ_best, γ_best
end

# coordinate descent algorithm for the regularization path (Λ x Γ)
function coordinate_descent(
    IT::iterator, d::dataCD, xdy, path, λ_max, Γ::Vector{Float64},
    y, lower_bounds, upper_bounds, w, verbose; sparse=0
    )

  # initializations
  #if sparse == 1
  #  β_path, β_tilde, β, activeSet = initSparse(IT)
  #else
  β_path, β_tilde, β, activeSet = initDense(IT)
  #end

  # memory allocation
  BIC = Inf::Float64
  β_ols = 0.0::Float64
  partial_fit = 0.0
  β_best = 0.0
  y_best = 0
  λ_best = 0.0
  γ_best = 0.0
  path_iteration = 0

  y_path = Vector{Float64}[]
  #for i in 1:TOTALCOMPONENTS
  #  push!(w, zeros(IT.nelements[i]))
  #end

  # regularization path
  @inbounds for γ in Γ

    # clear warm-start
    for i in IT.components
      for j in IT.elements[i]
         β_tilde[i][j] = 0.0
       end
    end

    @inbounds for path_iter in path

      λ = path_iter * λ_max

      change_flag = true

      if verbose
        path_iteration += 1
        print(string(
                "regularizers = (", round(λ, 3), ", ", round(γ, 3),
                ") path iteration = ", path_iteration, ".\n"
                ))
      end

      # loop until active set converges
      @inbounds for iter in 1:IT.maxIter

        println(iter)
        if !change_flag
          break
        end
        change_flag = false

        # cycle through every component
        @inbounds for c1 in IT.components
          @inbounds for j in IT.elements[c1]

            # compute the partial fit with the components in the active set
            partial_fit = 0.0
            @inbounds for c2 in IT.components
              @inbounds for l in IT.elements[c2]
                if activeSet[c2][l] && (c1, j) != (c2, l)
                  @inbounds partial_fit += GM2(c1, c2, j, l, d, IT) * β_tilde[c2][l]
                end
              end
            end

            # univariate ordinary leasts squares coefficient
            inner_prod_partial_residual = (xdy[c1][j] - partial_fit) / IT.obs

            # weighted penalty
            #w[c1][j] = 1.0 / (abs(β_ols)^γ)
#             w[c1][j] = 1.0 / abs(inner_prod_partial_residual)

            # soft thresholding operator
            if abs(inner_prod_partial_residual) <= λ[c1] * w[c1][j]^γ
              if activeSet[c1][j]
                β_tilde[c1][j] = 0.0
                activeSet[c1][j] = false
                change_flag = true
                #println("$(c1)  , $(j)")
              end
            else
              β_tilde[c1][j] = sign(inner_prod_partial_residual) * (abs(inner_prod_partial_residual) - λ[c1] * w[c1][j]^γ) / (d.σ[c1][j]^2)

              # projection onto the box constraints [lower_bound, upper_bound]
              β_tilde[c1][j] = max(β_tilde[c1][j], lower_bounds[c1])
              β_tilde[c1][j] = min(β_tilde[c1][j], upper_bounds[c1])

              if !activeSet[c1][j] #&& β_tilde[c1][j] != 0.0
                activeSet[c1][j] = true
                change_flag = true
                #println("$(c1)  , $(j)")
              end
            end

          end
        end
      end


      push!(β_path, deepcopy(β_tilde))
      # bayesian information criterion
      BIC_new, y_hat = compute_BIC(y, β_tilde, activeSet, IT, d, xdy)

      push!(y_path, copy(y_hat))
    
      if verbose
        print(string(" BIC = ", BIC_new))
      end
      
      # save the best fit so far
      if BIC_new < BIC
        BIC = BIC_new
        y_best = copy(y_hat)
        β_best = deepcopy(β_tilde)
        λ_best = λ
        γ_best = γ
      end

    end
  end

  return β_path, y_path, β_best, y_best, λ_best, γ_best
end


function update_components!(IT, w, β, d, xdy)

  # clear given iterator
  for c in IT.components
    IT.elements[c] = Int[]
  end

  for c in IT.components
    for j in 1:IT.nelements[c]

      if β[c][j] != 0.0
        # pre-compute weight
        w[c][j] = abs(1.0 / β[c][j]) #/ (d.σ[c][j]^2)
        β[c][j] = 0.0
        # update iterator only with nonzero elements
        push!(IT.elements[c], j)
      end

    end
  end

  return nothing
end
