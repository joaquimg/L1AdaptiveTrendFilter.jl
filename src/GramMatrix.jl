#ATENCAO
# considerar mu e sigma como input de cada função
# as outras entradas são i, j (numeros dos elementos)
# pode usar o N como numero total de elementos.
# eu estou usando uma struct chamada iterator (ver em CDtype.jl)
# posso corrigir isso sem problemas
# o principla é escrever as fómulas mesmo

#to complete

GM = Matrix{Function}(2,2)
GM[1, 1] = GM11
GM[1, 2] = GM12
GM[2, 1] = GM12
GM[1, 3] = GM13
GM[3, 1] = GM13
GM[1, 4] = GM14
GM[4, 1] = GM14
GM[1, 5] = GM15
GM[5, 1] = GM15

function GM11(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x step
  GM = Float64[]
  m = min(i, j)
  M = max(i, j)

  GM = (
    m*μ[m]*μ[M]-(M-m)*(1-μ[m])*μ[M]+(N-M)*(1-μ[m])*(1-μ[M])
    )/(σ[m]*σ[M])

  return GM
end

function GM12(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x slope
  GM = Float64[]

  if i > j
    GM = (
      j*μ[j]*μ[i] - (i-j)*(1-μ[j])*μ[i] - μ[i]*(N+1) + μ[i]*μ[j]*(N+1)
      + 0.5*(N+1)^2 - 0.5*N - 0.5*μ[j]*(N+1)^2 + 0.5*μ[j]*(N+1)
      + μ[i]*(i+1) - μ[i]*μ[j]*(i+1) - 0.5*(i+1)^2 + 0.5*i
      + 0.5*μ[j]*(i+1)^2 - 0.5*μ[j]*(i+1)
      )/(σ[i]*σ[j])
  else
    GM = (
      i*μ[i]*μ[j] - (i+1)*μ[i]*μ[j] + 0.5*μ[j]*(i+1)^2 - 0.5*μ[j]*(i+1)
      - μ[i]*(N+1) + μ[i]*μ[j]*(N+1) + 0.5*(N+1)^2 - 0.5*N-0.5*μ[j]*(N+1)^2
      + 0.5*μ[j]*(N+1) + μ[i]*(j+1) - 0.5*(j+1)^2 + 0.5*j
      )/(σ[i]*σ[j])
  end

  return GM
end

function GM13(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x spike
  GM = Float64[]

  if i > j
    GM = (
      j*μ[j]*μ[i] - (N-j)*(1-μ[j])*μ[i] + (1-μ[j])*μ[i] + (1-μ[j])*(1-μ[i])
      )/(σ[i]*σ[j])
  else
    GM = (
      j*μ[j]*μ[i] - μ[i]*μ[j] - μ[j]*(1-μ[i]) - (N-j)*(1-μ[j])*μ[i]
      )/(σ[i]*σ[j])
  end

  return GM
end

function GM14(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x sine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end

function GM15(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x cosine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end

function GM22(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # slope x slope
  GM = Float64[]
  m = min(i, j)
  M = max(i, j)

  GM = (
    (1/6)*N + m*μ[m]*μ[M] - (1/6)*M-0.5*(N+1)^2 + (1/3)*(N+1)^3 + 0.5*μ[M]*(m+1)^2
    - 0.5*μ[M]*(m+1) - μ[M]*(m+1)*m - μ[M]*μ[m]*(m+1) + 0.5*(M+1)^2 - (1/3)*(M+1)^3
    + 0.5*m*(M+1)^2 - 0.5*m(M+1) + 0.5*μ[m]*(M+1)^2 - 0.5*μ[m]*(M+1) + 0.5*M*(M+1)^2
    - 0.5*M*(M+1) - 0.5*m*(N+1)^2 + 0.5*m*(N+1) - 0.5*μ[m]*(N+1)^2 + 0.5*μ[m]*(N+1)
    -0.5*M*(N+1)^2 + 0.5*M*(N+1) - 0.5*μ[M]*(N+1)^2 + 0.5*μ[M]*(N+1) - (M+1)*n*M
    - (M+1)*μ[m]*M + (N+1)*m*M + μ[M]*(N+1)*m
  )/(σ[i]*σ[j])

  return GM
end


