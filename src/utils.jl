
function findλmax(IT,xdy,d)

  λmax   = 0.0
  temp  = 0.0

  maxPos = 0

  for i in IT.components

    #temp = maxabs(xdy[i].*d.σ[i])
    #k=indmax(abs(xdy[i].*d.σ[i]))
    temp,k=findmax(abs(xdy[i].*d.σ[i]))
    temp=temp/d.σ[i][k]
    if temp > λmax

      λmax = temp

    end

  end

  return λmax/IT.obs
end

function compute_λ_path(IT, xdy, numλ, d; logarit=true)

    λmax = findλmax(IT,xdy,d)

    if logarit
        vec = λmax*(logspace(1,0.001,numλ)-1.0)/9.0
    else
        vec = λmax*(linspace(1,0.001,numλ))
    end
    return vec
end

function compute_γ_path(IT, xdy, numγ, d; logarit=true)

    γmax = 1

    if logarit
        vec = γmax * (logspace(1, 0.001, numγ) - 1.0) / 9.0
    else
        vec = γmax*(linspace(1, 0.001, numγ))
    end
    return vec
end

function compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β, IT; ɛ = 1e-7::Float64)

    BIC = 0.0
    err = zeros(y)

    N = IT.obs

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

    BIC = N * log(var(err)) + 2 * k * log(N)

    return BIC

end
function compute_BIC(y::Vector{Float64}, β, IT, d; ɛ = 1e-5::Float64, std = 1)

  y_hat, β_new = compute_estimate(zeros(y), IT, β, d)

  if std == 1
    BIC = compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β_new, IT)
  else
    BIC = compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β, IT)
  end

  return BIC, y_hat
end

function compute_OLS(β_tilde, λ, activeSet,IT, xdy, d, lower_bounds, upper_bounds)

    β_ols = deepcopy(β_tilde)
    for c1 in IT.components
        for j in IT.elements[c1]
            # compute the partial fit with the components in the active set
            partial_fit = 0.0
            for c2 in IT.components
                for l in 1:size(activeSet[c2])[1]
                    if activeSet[c2][l]
                        partial_fit = partial_fit + GM2(c1,c2,j,l,d,IT) * β_tilde[c2][l]
                    end
                end
            end
            # univariate ordinary least squares coefficient
            β_ols[c1][j] = β_tilde[c1][j] + (1.0/IT.obs) * (xdy[c1][j] - partial_fit)

            # projection onto the box constraints [lower_bound, upper_bound]
            #β_ols[c1][j] = max(β_ols[c1][j], lower_bounds[c1])
            #β_ols[c1][j] = min(β_ols[c1][j], upper_bounds[c1])

            # hard thresholding operator
            if abs(β_ols[c1][j]) <= λ #* d.σ[c1][j]
                #β_ols[c1][j] = sign(β_ols[c1][j]) * (abs(β_ols[c1][j]) - λ*d.σ[c1][j])
                β_ols[c1][j] = 0.0
            else
                #β_ols[c1][j] = sign(β_ols[c1][j]) * (abs(β_ols[c1][j]) - λ)
            end
        end
    end
    return β_ols
end
