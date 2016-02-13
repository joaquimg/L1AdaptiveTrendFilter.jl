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
GM[2, 2] = GM22
GM[2, 3] = GM23
GM[3, 2] = GM23
GM[2, 4] = GM24
GM[4, 2] = GM24
GM[2, 5] = GM25
GM[5, 2] = GM25
GM[3, 3] = GM33
GM[3, 4] = GM34
GM[4, 3] = GM34
GM[3, 5] = GM35
GM[5, 3] = GM35
GM[4, 4] = GM44
GM[4, 5] = GM45
GM[5, 4] = GM45
GM[5, 5] = GM55

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
    μ[i]*μ[j]*(i+1) - 0.5*(μ[i]*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) + 0.5*μ[i]*sin((t+1)*ω)
    - μ[i]*μ[j] + 0.5*(μ[i]*sin(ω)*cos(ω))/(cos(ω)-1) - 0.5*μ[i]*sin(ω)
    + (μ[i]*μ[j]-μ[j])*(N+1) - 0.5*((μ[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*(μ[i]-1)*sin((N+1*ω)) - (μ[i]μ[j]-μ[j])*(i+1)
    + 0.5*((μ[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1) - 0.5*(μ[i]-1)*sin((i+1)*ω)
  )/(σ[i]*σ[j])

  return GM
end

function GM15(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x cosine
  GM = Float64[]

  GM = (
    μ[i]*μ[j]*(i+1) + 0.5*μ[i]*cos((i+1)*ω) + 0.5*(μ[i]*sin(ω)*sin((i+1)*ω))/(cos(ω)-1)
    - μ[i]*μ[j] - 0.5*μ[i]*cos(ω) - 0.5*(μ[i]*sin(ω)^2)/(cos(ω)-1)
    + (μ[i]*μ[j]-μ[j])*(N+1) - 0.5*((μ[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*(μ[i]-1)*sin((N+1)*ω) - (μ[i]*μ[j]-μ[j])*(i+1)
    + 0.5*((μ[i]-1)*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) - 0.5*(μ[i]-1)*sin((i+1)*ω)
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

function GM23(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # slope x spike
  GM = Float64[]

  if i > j
    GM = (
      j*μ[j]*μ[i] + μ[i]*(N+1)*j + μ[i]*μ[j]*(N+1) - 0.5*μ[i]*(N+1)^2 + 0.5*μ[i]*(N+1)
      - μ[i]*(j+1)*j - μ[i]*μ[j]*(j+1) + 0.5*μ[i]*(j+1)^2 - 0.5*μ[i]*(j+1)
      + (i-j-μ[j])*μ[i] + (i-j-μ[j])*(1-μ[i])
    )/(σ[i]*σ[j])
  else
    GM = (
      j*μ[i]*μ[j] + μ[i]*(N+1)*j + μ[i]*μ[j]*(N+1) - 0.5*μ[i]*(N+1)^2 + 0.5*μ[i]*(N+1)
      - μ[i]*(j+1)*j - μ[i]*μ[j]*(j+1) + 0.5*μ[i]*(j+1)^2 - 0.5*μ[i]*(j+1) - μ[i]*μ[j]
      - μ[j]*(1-μ[i])
    )/(σ[i]*σ[j])
  end

  return GM
end

function GM24(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # slope x sine
  GM = Float64[]

  GM = (
    μ[i]*μ[j]*(i+1) - 0.5*(μ[i]*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) + 0.5*μ[i]*sin((i+1)*ω)
    -μ[i]*μ[j] + 0.5*(μ[i]*sin(ω)*cos(ω))/(cos(ω)-1) - 0.5*μ[i]*sin(ω)
    + (i*μ[j]+μ[i]*μ[j]+0.5*μ[j])*(N+1) - 0.5*μ[j]*(N+1)^2
    + 0.5*(sin(ω)*(N+1)*cos((N+1)*ω))/(cos(ω)-1) - 0.5*(N+1)*sin((N+1)*ω)
    +0.5*((-i-μ[i])*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*((i*cos(ω)+μ[i]*cos(ω)-μ[i]-i-1)*sin((N+1)*ω))/(cos(ω)-1)
    - (i*μ[j]+μ[i]*μ[j]+0.5*μ[j])*(i+1) + 0.5*μ[j]*(i+1)^2
    - 0.5*(sin(ω)*(i+1)*cos((i+1)*ω))/(cos(ω)-1) + 0.5*(i+1)*sin((i+1)*ω)
    - 0.5*((-i-μ[i])*sin(ω)*cos((i+1)*ω))/(cos(ω)-1)
    - 0.5*((i*cos(ω)+μ[i]*cos(ω)-i-μ[i]-1)*sin((i+1)*ω)/(cos(ω)-1)
  )/(σ[i]*σ[j])

  return GM
end

function GM25(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # slope x cosine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end

function GM33(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # spike x spike
  GM = Float64[]

  if i != j
    GM = (
      (N-2)*μ[i]^2 - 2*μ[i]*(1-μ[i])
    )/(σ[i]*σ[j])
  else
    GM = 1
  end

  return GM
end

function GM34(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # spike x sine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end

function GM35(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # spike x cosine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end

function GM44(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # sine x sine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end

function GM45(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # sine x cosine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end

function GM55(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # cosine x cosine
  GM = Float64[]

  GM = (
  )/(σ[i]*σ[j])

  return GM
end
