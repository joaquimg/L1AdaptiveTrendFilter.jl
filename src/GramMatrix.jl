#ATENCAO
# considerar mu e sigma como input de cada função
# as outras entradas são i, j (numeros dos elementos)
# pode usar o N como numero total de elementos.
# eu estou usando uma struct chamada iterator (ver em CDtype.jl)
# posso corrigir isso sem problemas
# o principla é escrever as fómulas mesmo

#to complete



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
function GM11(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
    return GM11(i, j, size(d.μ[STEP]), d.μ[STEP], d.σ[STEP])
end
function GM12(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x slope
  GM = Float64[]

  if i > j
    GM = (
      j*μSLOPE[j]*μSTEP[i] - (i-j)*(1-μSLOPE[j])*μSTEP[i] - μSTEP[i]*(N+1) + μSTEP[i]*μSLOPE[j]*(N+1)
      + 0.5*(N+1)^2 - 0.5*N - 0.5*μSLOPE[j]*(N+1)^2 + 0.5*μSLOPE[j]*(N+1)
      + μSTEP[i]*(i+1) - μSTEP[i]*μSLOPE[j]*(i+1) - 0.5*(i+1)^2 + 0.5*i
      + 0.5*μSLOPE[j]*(i+1)^2 - 0.5*μSLOPE[j]*(i+1)
      )/(σ[i]*σ[j])
  else
    GM = (
      i*μSTEP[i]*μSLOPE[j] - (i+1)*μSTEP[i]*μSLOPE[j] + 0.5*μSLOPE[j]*(i+1)^2 - 0.5*μSLOPE[j]*(i+1)
      - μSTEP[i]*(N+1) + μSTEP[i]*μSLOPE[j]*(N+1) + 0.5*(N+1)^2 - 0.5*N-0.5*μSLOPE[j]*(N+1)^2
      + 0.5*μSLOPE[j]*(N+1) + μSTEP[i]*(j+1) - 0.5*(j+1)^2 + 0.5*j
      )/(σ[i]*σ[j])
  end

  return GM
end
function GM11(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
    return GM11(i, j, size(d.μ[STEP]), d.μ[STEP], d.σ[STEP])
end

function GM13(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x spike
  GM = Float64[]

  if i > j
    GM = (
      j*μSPIKE[j]*μSTEP[i] - (N-j)*(1-μSPIKE[j])*μSTEP[i] + (1-μSPIKE[j])*μSTEP[i] + (1-μSPIKE[j])*(1-μSTEP[i])
      )/(σ[i]*σ[j])
  else
    GM = (
      j*μSPIKE[j]*μSTEP[i] - μSTEP[i]*μSPIKE[j] - μSPIKE[j]*(1-μSTEP[i]) - (N-j)*(1-μSPIKE[j])*μSTEP[i]
      )/(σ[i]*σ[j])
  end

  return GM
end

function GM14(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x sine
  GM = Float64[]

  GM = (
    μSTEP[i]*μSIN[j]*(i+1) - 0.5*(μSTEP[i]*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) + 0.5*μSTEP[i]*sin((t+1)*ω)
    - μSTEP[i]*μSIN[j] + 0.5*(μSTEP[i]*sin(ω)*cos(ω))/(cos(ω)-1) - 0.5*μSTEP[i]*sin(ω)
    + (μSTEP[i]*μSIN[j]-μSIN[j])*(N+1) - 0.5*((μSTEP[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*(μSTEP[i]-1)*sin((N+1*ω)) - (μSTEP[i]μSIN[j]-μSIN[j])*(i+1)
    + 0.5*((μSTEP[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1) - 0.5*(μSTEP[i]-1)*sin((i+1)*ω)
  )/(σ[i]*σ[j])

  return GM
end

function GM15(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # step x cosine
  GM = Float64[]

  GM = (
    μSTEP[i]*μCOS[j]*(i+1) + 0.5*μSTEP[i]*cos((i+1)*ω) + 0.5*(μSTEP[i]*sin(ω)*sin((i+1)*ω))/(cos(ω)-1)
    - μSTEP[i]*μCOS[j] - 0.5*μSTEP[i]*cos(ω) - 0.5*(μSTEP[i]*sin(ω)^2)/(cos(ω)-1)
    + (μSTEP[i]*μCOS[j]-μCOS[j])*(N+1) - 0.5*((μSTEP[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*(μSTEP[i]-1)*sin((N+1)*ω) - (μSTEP[i]*μCOS[j]-μCOS[j])*(i+1)
    + 0.5*((μSTEP[i]-1)*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) - 0.5*(μSTEP[i]-1)*sin((i+1)*ω)
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
      j*μSPIKE[j]*μSLOPE[i] + μSLOPE[i]*(N+1)*j + μSLOPE[i]*μSPIKE[j]*(N+1) - 0.5*μSLOPE[i]*(N+1)^2 + 0.5*μSLOPE[i]*(N+1)
      - μSLOPE[i]*(j+1)*j - μSLOPE[i]*μSPIKE[j]*(j+1) + 0.5*μSLOPE[i]*(j+1)^2 - 0.5*μSLOPE[i]*(j+1)
      + (i-j-μSPIKE[j])*μSLOPE[i] + (i-j-μSPIKE[j])*(1-μSLOPE[i])
    )/(σ[i]*σ[j])
  else
    GM = (
      j*μSLOPE[i]*μSPIKE[j] + μSLOPE[i]*(N+1)*j + μSLOPE[i]*μSPIKE[j]*(N+1) - 0.5*μSLOPE[i]*(N+1)^2 + 0.5*μSLOPE[i]*(N+1)
      - μSLOPE[i]*(j+1)*j - μSLOPE[i]*μSPIKE[j]*(j+1) + 0.5*μSLOPE[i]*(j+1)^2 - 0.5*μSLOPE[i]*(j+1) - μSLOPE[i]*μSPIKE[j]
      - μSPIKE[j]*(1-μSLOPE[i])
    )/(σ[i]*σ[j])
  end

  return GM
end

function GM24(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # slope x sine
  GM = Float64[]

  GM = (
    μSLOPE[i]*μSIN[j]*(i+1) - 0.5*(μSLOPE[i]*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) + 0.5*μSLOPE[i]*sin((i+1)*ω)
    -μSLOPE[i]*μSIN[j] + 0.5*(μSLOPE[i]*sin(ω)*cos(ω))/(cos(ω)-1) - 0.5*μSLOPE[i]*sin(ω)
    + (i*μSIN[j]+μSLOPE[i]*μSIN[j]+0.5*μSIN[j])*(N+1) - 0.5*μSIN[j]*(N+1)^2
    + 0.5*(sin(ω)*(N+1)*cos((N+1)*ω))/(cos(ω)-1) - 0.5*(N+1)*sin((N+1)*ω)
    +0.5*((-i-μSLOPE[i])*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*((i*cos(ω)+μSLOPE[i]*cos(ω)-μSLOPE[i]-i-1)*sin((N+1)*ω))/(cos(ω)-1)
    - (i*μSIN[j]+μSLOPE[i]*μSIN[j]+0.5*μSIN[j])*(i+1) + 0.5*μSIN[j]*(i+1)^2
    - 0.5*(sin(ω)*(i+1)*cos((i+1)*ω))/(cos(ω)-1) + 0.5*(i+1)*sin((i+1)*ω)
    - 0.5*((-i-μSLOPE[i])*sin(ω)*cos((i+1)*ω))/(cos(ω)-1)
    - 0.5*((i*cos(ω)+μSLOPE[i]*cos(ω)-i-μSLOPE[i]-1)*sin((i+1)*ω)/(cos(ω)-1)
  ))/(σ[i]*σ[j])

  return GM
end

function GM25(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # slope x cosine
  GM = Float64[]

  GM = (
      μSLOPE[i]*μCOS[j]*(i+1) + 0.5*μSLOPE[i]*cos((i+1)*ω) + 0.5*(μSLOPE[i]*sin(ω)*sin((i+1)*ω))/(cos(ω)-1)
      -μSLOPE[i]*μCOS[j] - 0.5*μSLOPE[i]*cos(ω) -0.5*(μSLOPE[i]*sin(ω)^2)/(cos(ω)-1)#sera que eh mais
      + (i*μCOS[j]+μSLOPE[i]*μCOS[j]+0.5*μCOS[j])*(N+1)
      - 0.5*μCOS[j]*(N+1)^2 + 0.5*(sin(ω)*(N+1)*cos((N+1)*ω))/(cos(ω)-1)
      - 0.5*(N+1)sin((N+1)*ω) + 0.5*((-i-μSLOPE[i])*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
      + 0.5*(i*cos(ω)+μSLOPE[i]*cos(ω)-i-μSLOPE[i]-1)/(cos(ω)-1)-(i*μCOS[j]+μSLOPE[i]*μCOS[j]+0.5*μCOS[j])*(i+1)
      + 0.5*μCOS[j]*(i+1)^2 - 0.5*(sin(ω)*(i+1)*cos((i+1)*ω))/(cos(ω)-1)
      + 0.5*(i+1)*sin((i+1)*ω) - 0.5*((-i-μSLOPE[i])*sin(ω)*cos((i+1)*ω))/(cos(ω)-1)
      - 0.5*((i*cos(ω)+μSLOPE[i]*cos(ω)-i-μSLOPE[i]-1)*sin((i+1)*ω))/(cos(ω)-1)
  )/(σ[i]*σ[j])

  return GM
end

function GM33(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # spike x spike
  GM = Float64[]

  if i != j
    GM = (
      (N-2)*μSPIKE[i]^2 - 2*μSPIKE[i]*(1-μSPIKE[i])
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
      μSPIKE[i]*μSIN[j]*(N+1) - 0.5*(μSPIKE[i]*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
      + 0.5*μSPIKE[i]*sin((N+1)*ω) - μSPIKE[i]*μSIN[j] + 0.5*(μSPIKE[i]*sin(ω)cos(ω))/(cos(ω)-1)
      - 0.5*μSPIKE[i]*sin(ω) + μSPIKE[i]*(sin(i*ω)-μSIN[j]) + (1-μSPIKE[i])*(sin(i*ω)-μSIN[j])
  )/(σ[i]*σ[j])

  return GM
end

function GM35(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # spike x cosine
  GM = Float64[]

  GM = (
      μSPIKE[i]*μCOS[j]*(N+1) + 0.5*μSPIKE[i]*cos((N+1)*ω)
      + 0.5*(μSPIKE[i]*sin(ω)*sin((N+1)*ω))/(cos(ω)-1) -μSPIKE[i]*μCOS[j]
      - 0.5*μSPIKE[i]*cos(ω) - 0.5*(μSPIKE[i]*sin(ω)^2)/(cos(ω)-1)
      + μSPIKE[i]*(cos(i*ω)-μCOS[j]) + (1-μSPIKE[i])*(cos(i*ω)-μCOS[j])
  )/(σ[i]*σ[j])

  return GM
end

function GM44(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # sine x sine
  GM = Float64[]

  GM = (
      μSIN[i]*μSIN[j]*(N+1) - 0.5*(sin(j)*cos((N+1)*j)*sin(N+1)*i)/(cos(i)-cos(j))
      + 0.5*(sin(i)*sin((N+1)*j)*cos((N+1)*i))/(cos(i)-cos(j))
      - 0.5*sin((N+1)*j)*sin((N+1)*i) - 0.5*(μSIN[j]*sin(i)*cos((N+1)*i))/(cos(i)-1)
      - 0.5*(μSIN[i]*sin(j)*cos((N+1)*j))/(cos(j)-1) + 0.5*μSIN[j]*sin((N+1)*i)
      + 0.5*μSIN[i]*sin((N+1)*j) - μSIN[i]*μSIN[j] + 0.5*(sin(i)*sin(j)*cos(j))/(cos(i)-cos(j))
      - 0.5*(sin(j)*sin(i)*cos(i))/(cos(i)-cos(j)) + 0.5*sin(i)*sin(j)
      + 0.5*(μSIN[j]*sin(i)*cos(i))/(cos(i)-1) + 0.5*(μSIN[i]*sin(j)*cos(j))/(cos(j)-1)
      - 0.5*μSIN[j]*sin(i) - 0.5*μSIN[i]*sin(j)
  )/(σ[i]*σ[j])

  return GM
end

function GM45(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # sine x cosine
  GM = Float64[]

  GM = (
      μSINE[i]*μCOS[j]*(N+1) + 0.5*(sin(i)*cos((N+1)*j)*cos((N+1)*i))/(cos(i)-cos(j))
      - 0.5*cos((N+1)*j)*sin((N+1)*i)
      + 0.5*(sin(j)*sin((N+1)*j)*sin((N+1)*i))/(cos(i)-cos(j))
      - 0.5*(μCOS[j]*sin(i)*cos((N+1)*i))/(cos(i)-1) + 0.5*μSINE[i]cos((N+1)*j)
      + 0.5*μCOS[j]*sin((N+1)*i) + 0.5*(μSINE[i]*sin(j)*sin((N+1)*j))/(cos(j)-1)
      - μSINE[i]*μCOS[j] - 0.5*(sin(i)*cos(j)*cos(i))/(cos(i)-cos(j))
      + 0.5*cos(j)*sin(i) - 0.5*(sin(j)^2*sin(i))/(cos(i)-cos(j))
      + 0.5*(μCOS[j]*sin(i)*cos(i))/(cos(i)-1) - 0.5*μSINE[i]*cos(j)
      - 0.5*sin(i)*μCOS[j] - 0.5*(μSINE[i]*sin(j)^2)/(cos(j)-1)
  )/(σ[i]*σ[j])

  return GM
end

function GM55(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # cosine x cosine
  GM = Float64[]

  GM = (
      μCOS[i]*μCOS[j]*(N+1) - 0.5*cos((N+1)*j)*cos((N+1)*i)
      - 0.5*(sin(i)*cos((N+1)*j)*sin((N+1)*i))/(cos(i)-cos(j))
      + 0.5*(sin(j)*sin((N+1)*j)*cos((N+1)*i))/(cos(i)-cos(j))
      + 0.5*μCOS[j]*cos((N+1)*i) + 0.5*cos((N+1)*j)*μCOS[i]
      + 0.5*(μCOS[j]*sin(i)*sin((N+1)*i))/(cos(i)-1)
      + 0.5*(μCOS[i]*sin(j)*sin((N+1)*j))/(cos(j)-1) - μCOS[i]*μCOS[j]
      + 0.5*cos(i)*cos(j) + 0.5*(sin(i)^2*cos(j))/(cos(i)-cos(j))
      - 0.5*(sin(j)^2*cos(i))/(cos(i)-cos(j)) - 0.5*μCOS[j]*cos(i)
      - 0.5*μCOS[i]*cos(j) - 0.5*(μCOS[j]*sin(i)^2)/(cos(i)-1)
      - 0.5*(μCOS[i]*sin(j)^2)/(cos(j)-1)
  )/(σ[i]*σ[j])

  return GM
end

GM = Matrix{Function}(5,5)
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