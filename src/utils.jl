
function findλmax(IT,xdy,d)

  λmax   = 0.0
  temp  = 0.0
  maxPos = 0

  for i in IT.components

    temp, k = findmax(abs(xdy[i].*d.σ[i]))
    temp = temp / d.σ[i][k]

    if temp > λmax
      λmax = temp
    end
  end

  return λmax/IT.obs
end

function compute_λ_path(IT, xdy, numλ, d; logarit=true)

    λ_max = findλmax(IT,xdy,d)

    if logarit
        vec = λ_max*(logspace(1,0.001,numλ)-1.0)/9.0
    else
        vec = λ_max*(linspace(1,0.001,numλ))
    end
    return vec
end

function compute_γ_path(IT, xdy, numγ, d; logarit=true)

    γ_max = 2

    if logarit
        vec = γ_max * (logspace(1, 0.001, numγ) - 1.0) / 9.0
    else
        vec = γ_max * (linspace(1, 0.001, numγ))
    end
    return vec
end

function compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β, IT; ɛ = 1e-5::Float64)

  BIC = 0.0
  err = zeros(y)

  N = IT.obs

  for i in 1:IT.obs
    err[i] = y[i] - y_hat[i]
  end

  k = 0::Int64
  p = 0::Int64

  for i in IT.components
    for j in IT.elements[i]
      p += 1
      if abs(β[i][j]) != 0.0
        k += 1
      end
    end
  end

  BIC = N * log(var(err)) + k * log(N) + 2 * log(binomial(BigInt(p), BigInt(k)))

  return BIC
end
function compute_BIC(y::Vector{Float64}, β, activeSet, IT, d, xdy; ɛ = 1e-5::Float64, std = 1)

  β_ols = compute_OLS(β, activeSet, IT, xdy, d)
  y_hat, β_new = compute_estimate(zeros(y), IT, β_ols, d)

  if std == 1
    BIC = compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β_ols, IT)
  else
    BIC = compute_BIC(y_hat::Vector{Float64}, y::Vector{Float64}, β_ols, IT)
  end

  return BIC, y_hat
end

function compute_OLS(β, activeSet, IT, xdy, d)

  β_ols = deepcopy(β)

  @inbounds for c1 in IT.components
    @inbounds for j in IT.elements[c1]
      if β[c1][j] != 0.0
        partial_fit = 0.0
        @inbounds for c2 in IT.components
          @inbounds for l in IT.elements[c2]
            if β[c2][l] != 0.0 && (c1, j) != (c2, l)
              @inbounds partial_fit += GM2(c1, c2, j, l, d, IT) * β[c2][l]
            end
          end
        end

        # univariate ordinary leasts squares coefficient
        β_ols[c1][j] = (xdy[c1][j] - partial_fit) / (IT.obs * (d.σ[c1][j]^2))
      else
        β_ols[c1][j] = 0.0
      end

    end
  end

  return β_ols
end
