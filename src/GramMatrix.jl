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
function GM11(i::Int, j::Int, d, IT)
    return GM11(i, j, IT.obs, d.μ[STEP], d.σ[STEP])
end

function GM12(t::Int, l::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64}, μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64})
  # step x slope
  GM = Float64[]

  if l > t

    GM = (t*μSTEP[t]*μSLOPE[l]-(l-t)*(1-μSTEP[t])*μSLOPE[l]-(N+1)*l-μSLOPE[l]*(N+1)+(N+1)*l*μSTEP[t]+μSLOPE[l]*(N+1)*μSTEP[t]+(1/2)*(N+1)^2-(1/2)*N-(1/2)*μSTEP[t]*(N+1)^2+(1/2)*μSTEP[t]*(N+1)+(l+1)*l+μSLOPE[l]*(l+1)-(l+1)*l*μSTEP[t]-μSLOPE[l]*(l+1)*μSTEP[t]-(1/2)*(l+1)^2+(1/2)*l+(1/2)*μSTEP[t]*(l+1)^2-(1/2)*μSTEP[t]*(l+1))/(σSTEP[t]*σSLOPE[l])
    #GM = (
    #  j*μSLOPE[j]*μSTEP[i] - (i-j)*(1-μSLOPE[j])*μSTEP[i] - μSTEP[i]*(N+1) + μSTEP[i]*μSLOPE[j]*(N+1)
    #  + 0.5*(N+1)^2 - 0.5*N - 0.5*μSLOPE[j]*(N+1)^2 + 0.5*μSLOPE[j]*(N+1)
    #  + μSTEP[i]*(i+1) - μSTEP[i]*μSLOPE[j]*(i+1) - 0.5*(i+1)^2 + 0.5*i
    #  + 0.5*μSLOPE[j]*(i+1)^2 - 0.5*μSLOPE[j]*(i+1)
    #  )/(σSTEP[i]*σSLOPE[j])
  else

    GM = (-(l+1)*l*μSTEP[t]+(1/2)*t-(1/2)*N+(1/2)*(N+1)^2+μSLOPE[l]*(N+1)*μSTEP[t]-μSLOPE[l]*(l+1)*μSTEP[t]-μSLOPE[l]*(N+1)-(1/2)*μSTEP[t]*(N+1)^2+(1/2)*μSTEP[t]*(N+1)+(1/2)*μSTEP[t]*(l+1)^2-(1/2)*μSTEP[t]*(l+1)+μSLOPE[l]*(t+1)+l*μSTEP[t]*μSLOPE[l]-(1/2)*(t+1)^2+(N+1)*l*μSTEP[t]-(N+1)*l+(t+1)*l)/(σSTEP[t]*σSLOPE[l])
    #GM = (
    #  i*μSTEP[i]*μSLOPE[j] - (i+1)*μSTEP[i]*μSLOPE[j] + 0.5*μSLOPE[j]*(i+1)^2 - 0.5*μSLOPE[j]*(i+1)
    #  - μSTEP[i]*(N+1) + μSTEP[i]*μSLOPE[j]*(N+1) + 0.5*(N+1)^2 - 0.5*N-0.5*μSLOPE[j]*(N+1)^2
    #  + 0.5*μSLOPE[j]*(N+1) + μSTEP[i]*(j+1) - 0.5*(j+1)^2 + 0.5*j
    #  )/(σSTEP[i]*σSLOPE[j])
  end

  return GM
end

function GM12(i::Int, j::Int, d, IT)
    return GM12(i, j, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SLOPE], d.σ[SLOPE])
end
function GM21(i::Int, j::Int, d, IT)
    return GM12(j, i, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SLOPE], d.σ[SLOPE])
end

function GM13(i::Int, j::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64}, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64})
  # step x spike
  GM = Float64[]

  if j > i
    GM = (
      i*μSPIKE[j]*μSTEP[i] - (N-i)*(1-μSTEP[i])*μSPIKE[j] + (1-μSTEP[i])*μSPIKE[j] + (1-μSPIKE[j])*(1-μSTEP[i])
      )/(σSTEP[i]*σSPIKE[j])
  else
    GM = (
      i*μSPIKE[j]*μSTEP[i] - μSTEP[i]*μSPIKE[j] - μSTEP[i]*(1-μSPIKE[j]) - (N-i)*(1-μSTEP[i])*μSPIKE[j]
      )/(σSTEP[i]*σSPIKE[j])
  end

  return GM
end
function GM13(i::Int, j::Int, d, IT)
    return GM13(i, j, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SPIKE], d.σ[SPIKE])
end
function GM31(i::Int, j::Int, d, IT)
    return GM13(j, i, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SPIKE], d.σ[SPIKE])
end

function GM14(i::Int, j::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64}, μSIN::Vector{Float64}, σSIN::Vector{Float64},ω::Float64)
  # step x sine
  GM = Float64[]

  GM = (
    μSTEP[i]*μSIN[j]*(i+1) - 0.5*(μSTEP[i]*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) + 0.5*μSTEP[i]*sin((i+1)*ω)
    - μSTEP[i]*μSIN[j] + 0.5*(μSTEP[i]*sin(ω)*cos(ω))/(cos(ω)-1) - 0.5*μSTEP[i]*sin(ω)
    + (μSTEP[i]*μSIN[j]-μSIN[j])*(N+1) - 0.5*((μSTEP[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*(μSTEP[i]-1)*sin((N+1*ω)) - (μSTEP[i]μSIN[j]-μSIN[j])*(i+1)
    + 0.5*((μSTEP[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1) - 0.5*(μSTEP[i]-1)*sin((i+1)*ω)
  )/(σSTEP[i]*σSIN[j])

  return GM
end
function GM14(i::Int, j::Int, d, IT)
    return GM14(i, j, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SIN], d.σ[SIN],d.fs[j])
end
function GM41(i::Int, j::Int, d, IT)
    return GM14(j, i, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SIN], d.σ[SIN],d.fs[i])
end

function GM15(i::Int, j::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64}, μCOS::Vector{Float64}, σCOS::Vector{Float64},ω::Float64)
  # step x cosine

  ω = 0

  GM = Float64[]

  GM = (
    μSTEP[i]*μCOS[j]*(i+1) + 0.5*μSTEP[i]*cos((i+1)*ω) + 0.5*(μSTEP[i]*sin(ω)*sin((i+1)*ω))/(cos(ω)-1)
    - μSTEP[i]*μCOS[j] - 0.5*μSTEP[i]*cos(ω) - 0.5*(μSTEP[i]*sin(ω)^2)/(cos(ω)-1)
    + (μSTEP[i]*μCOS[j]-μCOS[j])*(N+1) - 0.5*((μSTEP[i]-1)*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
    + 0.5*(μSTEP[i]-1)*sin((N+1)*ω) - (μSTEP[i]*μCOS[j]-μCOS[j])*(i+1)
    + 0.5*((μSTEP[i]-1)*sin(ω)*cos((i+1)*ω))/(cos(ω)-1) - 0.5*(μSTEP[i]-1)*sin((i+1)*ω)
  )/(σSTEP[i]*σCOS[j])

  return GM
end
function GM15(i::Int, j::Int, d, IT)
    return GM15(i, j, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[COS], d.σ[COS],d.fc[j])
end
function GM51(i::Int, j::Int, d, IT)
    return GM15(j, i, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[COS], d.σ[COS],d.fc[i])
end

function GM22(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})
  # slope x slope
  GM = Float64[]
  m = min(i, j)
  M = max(i, j)

  #GM = (
  #  (1/6)*N + m*μ[m]*μ[M] - (1/6)*M-0.5*(N+1)^2 + (1/3)*(N+1)^3 + 0.5*μ[M]*(m+1)^2
  #  - 0.5*μ[M]*(m+1) - μ[M]*(m+1)*m - μ[M]*μ[m]*(m+1) + 0.5*(M+1)^2 - (1/3)*(M+1)^3
  #  + 0.5*m*(M+1)^2 - 0.5*m*(M+1) + 0.5*μ[m]*(M+1)^2 - 0.5*μ[m]*(M+1) + 0.5*M*(M+1)^2
  #  - 0.5*M*(M+1) - 0.5*m*(N+1)^2 + 0.5*m*(N+1) - 0.5*μ[m]*(N+1)^2 + 0.5*μ[m]*(N+1)
  #  -0.5*M*(N+1)^2 + 0.5*M*(N+1) - 0.5*μ[M]*(N+1)^2 + 0.5*μ[M]*(N+1) + (N+1)*m*M
  #  + (N+1)*μ[M]*m + (N+1)*μ[m]*M + μ[M]*(N+1)*μ[m]
  #)/(σ[i]*σ[j])

  GM = ((1/6)*N-(1/6)*M+m*μ[m]*μ[M]-μ[M]*(m+1)*μ[m]-μ[M]*(m+1)*m-(1/2)*(N+1)^2+(1/3)*(N+1)^3+(1/2)*(M+1)^2-(1/3)*(M+1)^3+(1/2)*μ[M]*(m+1)^2-(1/2)*μ[M]*(m+1)+(1/2)*m*(M+1)^2-(1/2)*m*(M+1)+(1/2)*μ[m]*(M+1)^2-(1/2)*μ[m]*(M+1)+(1/2)*M*(M+1)^2-(1/2)*M*(M+1)-(1/2)*m*(N+1)^2+(1/2)*m*(N+1)-(1/2)*μ[m]*(N+1)^2+(1/2)*μ[m]*(N+1)-(1/2)*M*(N+1)^2+(1/2)*M*(N+1)-(1/2)*μ[M]*(N+1)^2-(M+1)*m*M-(M+1)*μ[m]*M+(N+1)*m*M+μ[M]*(N+1)*m+(N+1)*μ[m]*M+μ[M]*(N+1)*μ[m]+(1/2)*μ[M]*(N+1))/(σ[m]*σ[M])

  return GM
end
function GM22(i::Int, j::Int, d, IT)
    return GM22(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE])
end

function GM23(l::Int, p::Int, N::Int, μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64}, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64})
  # slope x spike
  GM = Float64[]

  if p > l

    GM = (l*μSLOPE[l]*μSPIKE[p]+μSPIKE[p]*(N+1)*l+μSPIKE[p]*(N+1)*μSLOPE[l]-(1/2)*μSPIKE[p]*(N+1)^2+(1/2)*μSPIKE[p]*(N+1)-μSPIKE[p]*(l+1)*l-μSPIKE[p]*(l+1)*μSLOPE[l]+(1/2)*μSPIKE[p]*(l+1)^2-(1/2)*μSPIKE[p]*(l+1)+(p-l-μSLOPE[l])*μSPIKE[p]+(p-l-μSLOPE[l])*(1-μSPIKE[p]))/(σSLOPE[l]*σSPIKE[p])

    #GM = (
    #  j*μSPIKE[j]*μSLOPE[i] + μSLOPE[i]*(N+1)*j + μSLOPE[i]*μSPIKE[j]*(N+1) - 0.5*μSLOPE[i]*(N+1)^2 + 0.5*μSLOPE[i]*(N+1)
    #  - μSLOPE[i]*(j+1)*j - μSLOPE[i]*μSPIKE[j]*(j+1) + 0.5*μSLOPE[i]*(j+1)^2 - 0.5*μSLOPE[i]*(j+1)
    #  + (i-j-μSPIKE[j])*μSLOPE[i] + (i-j-μSPIKE[j])*(1-μSLOPE[i])
    #)/(σSLOPE[i]*σSPIKE[j])
  else
    
    GM = (l*μSLOPE[l]*μSPIKE[p]+μSPIKE[p]*(N+1)*l+μSPIKE[p]*(N+1)*μSLOPE[l]-(1/2)*μSPIKE[p]*(N+1)^2+(1/2)*μSPIKE[p]*(N+1)-μSPIKE[p]*(l+1)*l-μSPIKE[p]*(l+1)*μSLOPE[l]+(1/2)*μSPIKE[p]*(l+1)^2-(1/2)*μSPIKE[p]*(l+1)-μSLOPE[l]*μSPIKE[p]-μSLOPE[l]*(1-μSPIKE[p]))/(σSLOPE[l]*σSPIKE[p])

    #GM = (
    #  j*μSLOPE[i]*μSPIKE[j] + μSLOPE[i]*(N+1)*j + μSLOPE[i]*μSPIKE[j]*(N+1) - 0.5*μSLOPE[i]*(N+1)^2 + 0.5*μSLOPE[i]*(N+1)
    #  - μSLOPE[i]*(j+1)*j - μSLOPE[i]*μSPIKE[j]*(j+1) + 0.5*μSLOPE[i]*(j+1)^2 - 0.5*μSLOPE[i]*(j+1) - μSLOPE[i]*μSPIKE[j]
    #  - μSPIKE[j]*(1-μSLOPE[i])
    #)/(σSLOPE[i]*σSPIKE[j])
  end

  return GM
end
function GM23(i::Int, j::Int, d, IT)
    return GM23(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SPIKE], d.σ[SPIKE])
end
function GM32(i::Int, j::Int, d, IT)
    return GM23(j, i, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SPIKE], d.σ[SPIKE])
end

function GM24(i::Int, j::Int, N::Int, μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64}, μSIN::Vector{Float64}, σSIN::Vector{Float64},ω::Float64)
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
  ))/(σSLOPE[i]*σSIN[j])

  return GM
end
function GM24(i::Int, j::Int, d, IT)
    return GM24(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SIN], d.σ[SIN], d.fs[j])
end
function GM42(i::Int, j::Int, d, IT)
    return GM24(j, i, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SIN], d.σ[SIN], d.fs[i])
end

function GM25(i::Int, j::Int, N::Int, μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64}, μCOS::Vector{Float64}, σCOS::Vector{Float64},ω::Float64)
  # slope x cosine
  GM = Float64[]
  ω = 0

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
  )/(σSLOPE[i]*σCOS[j])

  return GM
end
function GM25(i::Int, j::Int, d, IT)
    return GM25(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[COS], d.σ[COS],d.fc[j])
end
function GM52(i::Int, j::Int, d, IT)
    return GM25(j, i, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[COS], d.σ[COS],d.fc[i])
end

function GM33(i::Int, j::Int, N::Int, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64})
  # spike x spike
  GM = Float64[]

  if i != j
    GM = (
      (N-2)*μSPIKE[i]^2 - 2*μSPIKE[i]*(1-μSPIKE[i])
    )/(σSPIKE[i]*σSPIKE[j])
  else
    GM = float(N)
  end

  return GM :: Float64
end
function GM33(i::Int, j::Int, d, IT)
    return GM33(i, j, IT.obs, d.μ[SPIKE], d.σ[SPIKE])
end

function GM34(i::Int, j::Int, N::Int, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64}, μSIN::Vector{Float64}, σSIN::Vector{Float64},ω::Float64)
  # spike x sine
  GM = Float64[]

  GM = (
      μSPIKE[i]*μSIN[j]*(N+1) - 0.5*(μSPIKE[i]*sin(ω)*cos((N+1)*ω))/(cos(ω)-1)
      + 0.5*μSPIKE[i]*sin((N+1)*ω) - μSPIKE[i]*μSIN[j] + 0.5*(μSPIKE[i]*sin(ω)cos(ω))/(cos(ω)-1)
      - 0.5*μSPIKE[i]*sin(ω) + μSPIKE[i]*(sin(i*ω)-μSIN[j]) + (1-μSPIKE[i])*(sin(i*ω)-μSIN[j])
  )/(σSPIKE[i]*σSIN[j])

  return GM
end
function GM34(i::Int, j::Int, d, IT)
    return GM34(i, j, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[SIN], d.σ[SIN], d.fs[j])
end
function GM43(i::Int, j::Int, d, IT)
    return GM34(j, i, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[SIN], d.σ[SIN], d.fs[i])
end

function GM35(i::Int, j::Int, N::Int, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64}, μCOS::Vector{Float64}, σCOS::Vector{Float64}, ω::Float64)
  # spike x cosine
  GM = Float64[]

  GM = (
      μSPIKE[i]*μCOS[j]*(N+1) + 0.5*μSPIKE[i]*cos((N+1)*ω)
      + 0.5*(μSPIKE[i]*sin(ω)*sin((N+1)*ω))/(cos(ω)-1) -μSPIKE[i]*μCOS[j]
      - 0.5*μSPIKE[i]*cos(ω) - 0.5*(μSPIKE[i]*sin(ω)^2)/(cos(ω)-1)
      + μSPIKE[i]*(cos(i*ω)-μCOS[j]) + (1-μSPIKE[i])*(cos(i*ω)-μCOS[j])
  )/(σSPIKE[i]*σCOS[j])

  return GM
end
function GM35(i::Int, j::Int, d, IT)
    return GM35(i, j, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[COS], d.σ[COS],d.fc[j])
end
function GM53(i::Int, j::Int, d, IT)
    return GM35(j, i, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[COS], d.σ[SPIKE],d.fc[i])
end

function GM44(i::Int, j::Int, N::Int, μSIN::Vector{Float64}, σSIN::Vector{Float64},ω::Vector{Float64})
  # sine x sine
  GM = Float64[]

  GM = (
      μSIN[i]*μSIN[j]*(N+1) - 0.5*(sin(ω[j])*cos((N+1)*ω[j])*sin(N+1)*ω[i])/(cos(ω[i])-cos(ω[j]))
      + 0.5*(sin(ω[i])*sin((N+1)*ω[j])*cos((N+1)*ω[i]))/(cos(ω[i])-cos(ω[j]))
      - 0.5*sin((N+1)*ω[j])*sin((N+1)*ω[i]) - 0.5*(μSIN[j]*sin(ω[i])*cos((N+1)*ω[i]))/(cos(ω[i])-1)
      - 0.5*(μSIN[i]*sin(ω[j])*cos((N+1)*ω[j]))/(cos(ω[j])-1) + 0.5*μSIN[j]*sin((N+1)*ω[i])
      + 0.5*μSIN[i]*sin((N+1)*ω[j]) - μSIN[i]*μSIN[j] + 0.5*(sin(ω[i])*sin(ω[j])*cos(j))/(cos(i)-cos(j))
      - 0.5*(sin(ω[j])*sin(ω[i])*cos(ω[i]))/(cos(ω[i])-cos(ω[j])) + 0.5*sin(ω[i])*sin(ω[j])
      + 0.5*(μSIN[j]*sin(ω[i])*cos(ω[i]))/(cos(ω[i])-1) + 0.5*(μSIN[i]*sin(ω[j])*cos(ω[j]))/(cos(ω[j])-1)
      - 0.5*μSIN[j]*sin(ω[i]) - 0.5*μSIN[i]*sin(ω[j])
  )/(σSIN[i]*σSIN[j])

  return GM
end
function GM44(i::Int, j::Int, d, IT)
    return GM44(i, j, IT.obs, d.μ[SIN], d.σ[SIN], d.fs)
end

#CORRIGIG INPUT DOS SENOIDES
function GM45(i::Int, j::Int, N::Int, μSINE::Vector{Float64}, σSIN::Vector{Float64}, μCOS::Vector{Float64}, σCOS::Vector{Float64}, ωs::Vector{Float64}, ωc::Vector{Float64})
  # sine x cosine
  GM = Float64[]

  GM = (
      μSINE[i]*μCOS[j]*(N+1) + 0.5*(sin(ωs[i])*cos((N+1)*ωc[j])*cos((N+1)*ωs[i]))/(cos(ωs[i])-cos(ωc[j]))
      - 0.5*cos((N+1)*ωc[j])*sin((N+1)*ωs[i])
      + 0.5*(sin(ωc[j])*sin((N+1)*ωc[j])*sin((N+1)*ωs[i]))/(cos(ωs[i])-cos(ωc[j]))
      - 0.5*(μCOS[j]*sin(ωs[i])*cos((N+1)*ωs[i]))/(cos(ωs[i])-1) + 0.5*μSINE[i]*cos((N+1)*ωc[j])
      + 0.5*μCOS[j]*sin((N+1)*ωs[i]) + 0.5*(μSINE[i]*sin(ωc[j])*sin((N+1)*ωc[j]))/(cos(ωc[j])-1)
      - μSINE[i]*μCOS[j] - 0.5*(sin(ωs[i])*cos(ωc[j])*cos(ωs[i]))/(cos(ωs[i])-cos(ωc[j]))
      + 0.5*cos(ωc[j])*sin(ωs[i]) - 0.5*(sin(ωc[j])^2*sin(ωs[i]))/(cos(ωs[i])-cos(ωc[j]))
      + 0.5*(μCOS[j]*sin(ωs[i])*cos(ωs[i]))/(cos(ωs[i])-1) - 0.5*μSINE[i]*cos(ωc[j])
      - 0.5*sin(ωs[i])*μCOS[j] - 0.5*(μSINE[i]*sin(ωc[j])^2)/(cos(ωc[j])-1)
  )/(σSIN[i]*σCOS[j])

  return GM
end
function GM45(i::Int, j::Int, d, IT)
    return GM45(i, j, IT.obs, d.μ[SIN], d.σ[SIN], d.μ[COS], d.σ[COS],d.fs,d.fc)
end
function GM54(i::Int, j::Int, d, IT)
    return GM45(j, i, IT.obs, d.μ[SIN], d.σ[SIN], d.μ[COS], d.σ[COS],d.fs,d.fc)
end

function GM55(i::Int, j::Int, N::Int, μCOS::Vector{Float64}, σCOS::Vector{Float64},ω::Vector{Float64})
  # cosine x cosine
  GM = Float64[]

  GM = (
      μCOS[i]*μCOS[j]*(N+1) - 0.5*cos((N+1)*ω[j])*cos((N+1)*ω[i])
      - 0.5*(sin(ω[i])*cos((N+1)*ω[j])*sin((N+1)*ω[i]))/(cos(ω[i])-cos(ω[j]))
      + 0.5*(sin(ω[j])*sin((N+1)*ω[j])*cos((N+1)*ω[i]))/(cos(ω[i])-cos(ω[j]))
      + 0.5*μCOS[j]*cos((N+1)*ω[i]) + 0.5*cos((N+1)*ω[j])*μCOS[i]
      + 0.5*(μCOS[j]*sin(ω[i])*sin((N+1)*ω[i]))/(cos(ω[i])-1)
      + 0.5*(μCOS[i]*sin(ω[j])*sin((N+1)*ω[j]))/(cos(ω[j])-1) - μCOS[i]*μCOS[j]
      + 0.5*cos(ω[i])*cos(ω[j]) + 0.5*(sin(ω[i])^2*cos(ω[j]))/(cos(ω[i])-cos(ω[j]))
      - 0.5*(sin(ω[j])^2*cos(ω[i]))/(cos(ω[i])-cos(ω[j])) - 0.5*μCOS[j]*cos(ω[i])
      - 0.5*μCOS[i]*cos(ω[j]) - 0.5*(μCOS[j]*sin(ω[i])^2)/(cos(ω[i])-1)
      - 0.5*(μCOS[i]*sin(ω[j])^2)/(cos(ω[j])-1)
  )/(σCOS[i]*σCOS[j])

  return GM
end
function GM55(i::Int, j::Int, d, IT)
    return GM55(i, j, IT.obs, d.μ[COS], d.σ[COS],d.fc)
end

GM = Matrix{Function}(5, 5)
GM[1, 1] = GM11
GM[1, 2] = GM12
GM[2, 1] = GM21
GM[1, 3] = GM13
GM[3, 1] = GM31
GM[1, 4] = GM14
GM[4, 1] = GM41
GM[1, 5] = GM15
GM[5, 1] = GM51
GM[2, 2] = GM22
GM[2, 3] = GM23
GM[3, 2] = GM32
GM[2, 4] = GM24
GM[4, 2] = GM42
GM[2, 5] = GM25
GM[5, 2] = GM52
GM[3, 3] = GM33
GM[3, 4] = GM34
GM[4, 3] = GM43
GM[3, 5] = GM35
GM[5, 3] = GM53
GM[4, 4] = GM44
GM[4, 5] = GM45
GM[5, 4] = GM54
GM[5, 5] = GM55

