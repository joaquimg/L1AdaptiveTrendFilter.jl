# Closed-form formulas for the inner products between different kinds of components

# step x step
function GM11(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})

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

# step x slope
function GM12(
    t::Int, l::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64},
    μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64}
    )

  if l > t
    GM = (
      t*μSTEP[t]*μSLOPE[l] - (l-t)*(1-μSTEP[t])*μSLOPE[l] - (N+1)*l-μSLOPE[l]*(N+1)
      + (N+1)*l*μSTEP[t] + μSLOPE[l]*(N+1)*μSTEP[t] + (1/2)*(N+1)^2-(1/2)*N
      - (1/2)*μSTEP[t]*(N+1)^2 + (1/2)*μSTEP[t]*(N+1) + (l+1)*l+μSLOPE[l]*(l+1)
      - (l+1)*l*μSTEP[t] - μSLOPE[l]*(l+1)*μSTEP[t] - (1/2)*(l+1)^2+(1/2)*l
      + (1/2)*μSTEP[t]*(l+1)^2 - (1/2)*μSTEP[t]*(l+1)
      )/(σSTEP[t]*σSLOPE[l])
  else
    GM = (
      - (l+1)*l*μSTEP[t]+(1/2)*t - (1/2)*N+(1/2)*(N+1)^2 + μSLOPE[l]*(N+1)*μSTEP[t]
      - μSLOPE[l]*(l+1)*μSTEP[t] - μSLOPE[l]*(N+1) - (1/2)*μSTEP[t]*(N+1)^2
      + (1/2)*μSTEP[t]*(N+1) + (1/2)*μSTEP[t]*(l+1)^2 - (1/2)*μSTEP[t]*(l+1)
      + μSLOPE[l]*(t+1) + l*μSTEP[t]*μSLOPE[l] - (1/2)*(t+1)^2 + (N+1)*l*μSTEP[t]
      - (N+1)*l+(t+1)*l
      )/(σSTEP[t]*σSLOPE[l])
  end

  return GM
end
function GM12(i::Int, j::Int, d, IT)
    return GM12(i, j, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SLOPE], d.σ[SLOPE])
end
function GM21(i::Int, j::Int, d, IT)
    return GM12(j, i, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SLOPE], d.σ[SLOPE])
end

# step x spike
function GM13(
    i::Int, j::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64},
    μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64}
    )

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

# step x sine
function GM14(
    t::Int, s::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64},
    μSIN::Vector{Float64}, σSIN::Vector{Float64},ω::Float64
    )

  GM = (μSTEP[t]*μSIN[s]*(t+1)-(1/2)*μSTEP[t]*sin(ω)*cos((t+1)*ω)/(cos(ω)-1)
        +(1/2)*μSTEP[t]*sin((t+1)*ω)-μSTEP[t]*μSIN[s]+(1/2)*μSTEP[t]*sin(ω)*cos(ω)/(cos(ω)-1)
        -(1/2)*μSTEP[t]*sin(ω)+(μSTEP[t]*μSIN[s]-μSIN[s])*(N+1)-(1/2)*(-1+μSTEP[t])*sin(ω)*cos((N+1)*ω)/(cos(ω)-1)
        +(-1/2+(1/2)*μSTEP[t])*sin((N+1)*ω)-(μSTEP[t]*μSIN[s]-μSIN[s])*(t+1)
        +(1/2)*(-1+μSTEP[t])*sin(ω)*cos((t+1)*ω)/(cos(ω)-1)-(-1/2+(1/2)*μSTEP[t])*sin((t+1)*ω)
        )/(σSTEP[t]*σSIN[s]
        )

  return GM
end
function GM14(i::Int, j::Int, d, IT)
    return GM14(i, j, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SIN], d.σ[SIN],d.fs[j])
end
function GM41(i::Int, j::Int, d, IT)
    return GM14(j, i, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[SIN], d.σ[SIN],d.fs[i])
end

# step x cosine
function GM15(
    t::Int, c::Int, N::Int, μSTEP::Vector{Float64}, σSTEP::Vector{Float64},
    μCOS::Vector{Float64}, σCOS::Vector{Float64},ω::Float64
    )

  GM = (μSTEP[t]*μCOS[c]*(t+1)+(1/2)*μSTEP[t]*cos((t+1)*ω)+(1/2)*μSTEP[t]*sin(ω)*sin((t+1)*ω)/(cos(ω)-1)
        -μSTEP[t]*μCOS[c]-(1/2)*μSTEP[t]*cos(ω)-(1/2)*μSTEP[t]*sin(ω)^2/(cos(ω)-1)
        +(μSTEP[t]*μCOS[c]-μCOS[c])*(N+1)-(1/2)*(-1+μSTEP[t])*sin(ω)*cos((N+1)*ω)/(cos(ω)-1)
        +(-1/2+(1/2)*μSTEP[t])*sin((N+1)*ω)-(μSTEP[t]*μCOS[c]-μCOS[c])*(t+1)
        +(1/2)*(-1+μSTEP[t])*sin(ω)*cos((t+1)*ω)/(cos(ω)-1)-(-1/2+(1/2)*μSTEP[t])*sin((t+1)*ω)
        )/(σSTEP[t]*σCOS[c]
        )
  
  return GM
end
function GM15(i::Int, j::Int, d, IT)
    return GM15(i, j, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[COS], d.σ[COS],d.fc[j])
end
function GM51(i::Int, j::Int, d, IT)
    return GM15(j, i, IT.obs, d.μ[STEP], d.σ[STEP], d.μ[COS], d.σ[COS],d.fc[i])
end

# slope x slope
function GM22(i::Int, j::Int, N::Int, μ::Vector{Float64}, σ::Vector{Float64})

  m = min(i, j)
  M = max(i, j)

  GM = (
    (1/6)*N - (1/6)*M + m*μ[m]*μ[M] - μ[M]*(m+1)*μ[m] - μ[M]*(m+1)*m
    - (1/2)*(N+1)^2 + (1/3)*(N+1)^3 + (1/2)*(M+1)^2 - (1/3)*(M+1)^3
    + (1/2)*μ[M]*(m+1)^2 - (1/2)*μ[M]*(m+1) + (1/2)*m*(M+1)^2
    - (1/2)*m*(M+1) + (1/2)*μ[m]*(M+1)^2 - (1/2)*μ[m]*(M+1)
    + (1/2)*M*(M+1)^2 - (1/2)*M*(M+1) - (1/2)*m*(N+1)^2 + (1/2)*m*(N+1)
    - (1/2)*μ[m]*(N+1)^2 + (1/2)*μ[m]*(N+1) - (1/2)*M*(N+1)^2
    + (1/2)*M*(N+1) - (1/2)*μ[M]*(N+1)^2 - (M+1)*m*M
    - (M+1)*μ[m]*M + (N+1)*m*M + μ[M]*(N+1)*m + (N+1)*μ[m]*M
    + μ[M]*(N+1)*μ[m] + (1/2)*μ[M]*(N+1)
    )/(σ[m]*σ[M])

  return GM
end
function GM22(i::Int, j::Int, d, IT)
    return GM22(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE])
end

# slope x spike
function GM23(
    l::Int, p::Int, N::Int, μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64},
    μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64}
    )

  if p > l
    GM = (
      l*μSLOPE[l]*μSPIKE[p] + μSPIKE[p]*(N+1)*l + μSPIKE[p]*(N+1)*μSLOPE[l]
      - (1/2)*μSPIKE[p]*(N+1)^2 + (1/2)*μSPIKE[p]*(N+1) - μSPIKE[p]*(l+1)*l
      - μSPIKE[p]*(l+1)*μSLOPE[l] + (1/2)*μSPIKE[p]*(l+1)^2
      - (1/2)*μSPIKE[p]*(l+1) + (p-l-μSLOPE[l])*μSPIKE[p]
      + (p-l-μSLOPE[l])*(1-μSPIKE[p])
      )/(σSLOPE[l]*σSPIKE[p])
  else
    GM = (
      l*μSLOPE[l]*μSPIKE[p] + μSPIKE[p]*(N+1)*l + μSPIKE[p]*(N+1)*μSLOPE[l]
      - (1/2)*μSPIKE[p]*(N+1)^2 + (1/2)*μSPIKE[p]*(N+1) - μSPIKE[p]*(l+1)*l
      - μSPIKE[p]*(l+1)*μSLOPE[l] + (1/2)*μSPIKE[p]*(l+1)^2
      - (1/2)*μSPIKE[p]*(l+1) - μSLOPE[l]*μSPIKE[p] - μSLOPE[l]*(1-μSPIKE[p])
      )/(σSLOPE[l]*σSPIKE[p])
  end

  return GM
end
function GM23(i::Int, j::Int, d, IT)
    return GM23(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SPIKE], d.σ[SPIKE])
end
function GM32(i::Int, j::Int, d, IT)
    return GM23(j, i, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SPIKE], d.σ[SPIKE])
end

# slope x sine
function GM24(
    l::Int, s::Int, N::Int, μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64},
    μSIN::Vector{Float64}, σSIN::Vector{Float64},ω::Float64
    )

  GM = (μSLOPE[l]*μSIN[s]*(l+1)-(1/2)*μSLOPE[l]*sin(ω)*cos((l+1)*ω)/(cos(ω)-1)
        +(1/2)*μSLOPE[l]*sin((l+1)*ω)-μSLOPE[l]*μSIN[s]+(1/2)*μSLOPE[l]*sin(ω)*cos(ω)/(cos(ω)-1)
        -(1/2)*μSLOPE[l]*sin(ω)+(l*μSIN[s]+(1/2)*μSIN[s]+μSLOPE[l]*μSIN[s])*(N+1)
        -(1/2)*μSIN[s]*(N+1)^2+(1/2)*sin(ω)*(N+1)*cos((N+1)*ω)/(cos(ω)-1)-(1/2*(N+1))*sin((N+1)*ω)
        -(1/2)*(μSLOPE[l]+l)*sin(ω)*cos((N+1)*ω)/(cos(ω)-1)
        +(1/2)*(l*cos(ω)+μSLOPE[l]*cos(ω)-l-μSLOPE[l]-1)*sin((N+1)*ω)/(cos(ω)-1)
        -(l*μSIN[s]+(1/2)*μSIN[s]+μSLOPE[l]*μSIN[s])*(l+1)
        +(1/2)*μSIN[s]*(l+1)^2-(1/2)*sin(ω)*(l+1)*cos((l+1)*ω)/(cos(ω)-1)
        +(1/2*(l+1))*sin((l+1)*ω)
        +(1/2)*(μSLOPE[l]+l)*sin(ω)*cos((l+1)*ω)/(cos(ω)-1)
        -(1/2)*(l*cos(ω)+μSLOPE[l]*cos(ω)-l-μSLOPE[l]-1)*sin((l+1)*ω)/(cos(ω)-1)
        )/(σSLOPE[l]*σSIN[s]
        )
  
  return GM
end
function GM24(i::Int, j::Int, d, IT)
    return GM24(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SIN], d.σ[SIN], d.fs[j])
end
function GM42(i::Int, j::Int, d, IT)
    return GM24(j, i, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[SIN], d.σ[SIN], d.fs[i])
end

# slope x cosine
function GM25(
    l::Int, c::Int, N::Int, μSLOPE::Vector{Float64}, σSLOPE::Vector{Float64},
    μCOS::Vector{Float64}, σCOS::Vector{Float64}, ω::Float64
    )

  GM = (μSLOPE[l]*μCOS[c]*(l+1)+(1/2)*μSLOPE[l]*cos((l+1)*ω)
        +(1/2)*μSLOPE[l]*sin(ω)*sin((l+1)*ω)/(cos(ω)-1)-μSLOPE[l]*μCOS[c]
        -(1/2)*μSLOPE[l]*cos(ω)-(1/2)*μSLOPE[l]*sin(ω)^2/(cos(ω)-1)
        +((1/2)*μCOS[c]+l*μCOS[c]+μSLOPE[l]*μCOS[c])*(N+1)
        -(1/2)*μCOS[c]*(N+1)^2+(1/2)*sin(ω)*(N+1)*cos((N+1)*ω)/(cos(ω)-1)
        -(1/2*(N+1))*sin((N+1)*ω)-(1/2)*(μSLOPE[l]+l)*sin(ω)*cos((N+1)*ω)/(cos(ω)-1)
        +(1/2)*(l*cos(ω)+μSLOPE[l]*cos(ω)-l-μSLOPE[l]-1)*sin((N+1)*ω)/(cos(ω)-1)
        -((1/2)*μCOS[c]+l*μCOS[c]+μSLOPE[l]*μCOS[c])*(l+1)+(1/2)*μCOS[c]*(l+1)^2
        -(1/2)*sin(ω)*(l+1)*cos((l+1)*ω)/(cos(ω)-1)+(1/2*(l+1))*sin((l+1)*ω)
        +(1/2)*(μSLOPE[l]+l)*sin(ω)*cos((l+1)*ω)/(cos(ω)-1)-(1/2)*(l*cos(ω)
        +μSLOPE[l]*cos(ω)-l-μSLOPE[l]-1)*sin((l+1)*ω)/(cos(ω)-1)
        )/(σSLOPE[l]*σCOS[c])

  return GM
end
function GM25(i::Int, j::Int, d, IT)
    return GM25(i, j, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[COS], d.σ[COS],d.fc[j])
end
function GM52(i::Int, j::Int, d, IT)
    return GM25(j, i, IT.obs, d.μ[SLOPE], d.σ[SLOPE], d.μ[COS], d.σ[COS],d.fc[i])
end

# spike x spike
function GM33(i::Int, j::Int, N::Int, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64})

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

# spike x sine
function GM34(
    p::Int, s::Int, N::Int, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64},
    μSIN::Vector{Float64}, σSIN::Vector{Float64}, ω::Float64
    )

  GM = (μSPIKE[p]*μSIN[s]*(N+1)-(1/2)*μSPIKE[p]*sin(ω)*cos((N+1)*ω)/(cos(ω])-1)
        +(1/2)*μSPIKE[p]*sin((N+1)*ω)-μSPIKE[p]*μSIN[s]+(1/2)*μSPIKE[p]*sin(ω)*cos(f[s])/(cos(f[s])-1)
        -(1/2)*μSPIKE[p]*sin(ω)+μSPIKE[p]*(sin(p*ω)-μSIN[s])+(1-μSPIKE[p])*(sin(p*ω)-μSIN[s])
        )/(σSPIKE[p]*σSIN[s])

  return GM
end
function GM34(i::Int, j::Int, d, IT)
    return GM34(i, j, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[SIN], d.σ[SIN], d.fs[j])
end
function GM43(i::Int, j::Int, d, IT)
    return GM34(j, i, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[SIN], d.σ[SIN], d.fs[i])
end

# spike x cosine
function GM35(
    p::Int, c::Int, N::Int, μSPIKE::Vector{Float64}, σSPIKE::Vector{Float64},
    μCOS::Vector{Float64}, σCOS::Vector{Float64}, ω::Float64
    )

  GM = (μSPIKE[p]*μCOS[c]*(N+1)+(1/2)*μSPIKE[p]*cos((N+1)*ω)
        +(1/2)*μSPIKE[p]*sin(ω)*sin((N+1)*ω)/(cos(ω)-1)-μSPIKE[p]*μCOS[c]
        -(1/2)*μSPIKE[p]*cos(ω)-(1/2)*μSPIKE[p]*sin(ω)^2/(cos(ω)-1)
        +μSPIKE[p]*(cos(p*ω)-μCOS[c])+(1-μSPIKE[p])*(cos(p*ω)-μCOS[c])
        )/(σSPIKE[p]*σCOS[c])

  return GM
end
function GM35(i::Int, j::Int, d, IT)
    return GM35(i, j, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[COS], d.σ[COS],d.fc[j])
end
function GM53(i::Int, j::Int, d, IT)
    return GM35(j, i, IT.obs, d.μ[SPIKE], d.σ[SPIKE], d.μ[COS], d.σ[SPIKE],d.fc[i])
end

# sine x sine
function GM44(i::Int, j::Int, N::Int, μSIN::Vector{Float64}, σSIN::Vector{Float64},ω::Vector{Float64})

  #ATENTION HERE
  #there still may be a ptoblem from cos(ωi) = cos(ωj) with i != j
  if i == j
    GM = float(N)
  elseif cos(ω[j])-cos(ω[i]) == 0
    #cancelation that need to happen for finiteness
    GM = (μSIN[j]*μSIN[i]*(N+1)
          -(1/2)*sin((N+1)*ω[j])*sin((N+1)*ω[i])-(1/2)*μSIN[j]*sin(ω[i])*cos((N+1)*ω[i])/(-1+cos(ω[i]))
          -(1/2)*μSIN[i]*sin(ω[j])*cos((N+1)*ω[j])/(cos(ω[j])-1)+(1/2)*μSIN[j]*sin((N+1)*ω[i])
          +(1/2)*sin((N+1)*ω[j])*μSIN[i]-μSIN[j]*μSIN[i]
          +(1/2)*sin(ω[i])*sin(ω[j])
          +(1/2)*μSIN[j]*sin(ω[i])*cos(ω[i])/(-1+cos(ω[i]))+(1/2)*μSIN[i]*sin(ω[j])*cos(ω[j])/(cos(ω[j])-1)
          -(1/2)*μSIN[j]*sin(ω[i])-(1/2)*μSIN[i]*sin(ω[j])
          )/(σSIN[j]*σSIN[i])
  else
    GM = (μSIN[j]*μSIN[i]*(N+1)+(1/2)*sin(ω[j])*cos((N+1)*ω[j])*sin((N+1)*ω[i])/(cos(ω[j])-cos(ω[i]))
          -(1/2)*sin(ω[i])*sin((N+1)*ω[j])*cos((N+1)*ω[i])/(cos(ω[j])-cos(ω[i]))
          -(1/2)*sin((N+1)*ω[j])*sin((N+1)*ω[i])-(1/2)*μSIN[j]*sin(ω[i])*cos((N+1)*ω[i])/(-1+cos(ω[i]))
          -(1/2)*μSIN[i]*sin(ω[j])*cos((N+1)*ω[j])/(cos(ω[j])-1)+(1/2)*μSIN[j]*sin((N+1)*ω[i])
          +(1/2)*sin((N+1)*ω[j])*μSIN[i]-μSIN[j]*μSIN[i]-(1/2)*sin(ω[j])*cos(ω[j])*sin(ω[i])/(cos(ω[j])-cos(ω[i]))
          +(1/2)*sin(ω[i])*sin(ω[j])*cos(ω[i])/(cos(ω[j])-cos(ω[i]))+(1/2)*sin(ω[i])*sin(ω[j])
          +(1/2)*μSIN[j]*sin(ω[i])*cos(ω[i])/(-1+cos(ω[i]))+(1/2)*μSIN[i]*sin(ω[j])*cos(ω[j])/(cos(ω[j])-1)
          -(1/2)*μSIN[j]*sin(ω[i])-(1/2)*μSIN[i]*sin(ω[j])
          )/(σSIN[j]*σSIN[i])
  end

  return GM
end
function GM44(i::Int, j::Int, d, IT)
    return GM44(i, j, IT.obs, d.μ[SIN], d.σ[SIN], d.fs)
end


# sine x cosine
function GM45(
    s::Int, c::Int, N::Int, μSIN::Vector{Float64}, σSIN::Vector{Float64},
    μCOS::Vector{Float64}, σCOS::Vector{Float64}, ωs::Vector{Float64}, ωc::Vector{Float64}
    )

  if cos(ωc[c])-cos(ωs[s]) != 0
    GM = (μSIN[s]*μCOS[c]*(N+1)+(1/2)*sin(ωs[s])*cos((N+1)*ωs[s])*cos((N+1)*ωc[c])/(cos(ωs[s])-cos(ωc[c]))
          -(1/2)*sin((N+1)*ωs[s])*cos((N+1)*ωc[c])+(1/2)*sin(ωc[c])*sin((N+1)*ωs[s])*sin((N+1)*ωc[c])/(cos(ωs[s])-cos(ωc[c]))
          +(1/2)*μSIN[s]*cos((N+1)*ωc[c])-(1/2)*μCOS[c]*sin(ωs[s])*cos((N+1)*ωs[s])/(cos(ωs[s])-1)
          +(1/2)*μSIN[s]*sin(ωc[c])*sin((N+1)*ωc[c])/(cos(ωc[c])-1)+(1/2)*sin((N+1)*ωs[s])*μCOS[c]-μSIN[s]*μCOS[c]
          -(1/2)*sin(ωs[s])*cos(ωs[s])*cos(ωc[c])/(cos(ωs[s])-cos(ωc[c]))+(1/2)*sin(ωs[s])*cos(ωc[c])
          -(1/2)*sin(ωc[c])^2*sin(ωs[s])/(cos(ωs[s])-cos(ωc[c]))-(1/2)*μSIN[s]*cos(ωc[c])
          +(1/2)*μCOS[c]*sin(ωs[s])*cos(ωs[s])/(cos(ωs[s])-1)-(1/2)*μSIN[s]*sin(ωc[c])^2/(cos(ωc[c])-1)
          -(1/2)*sin(μSIN[s])*μCOS[c])/(σCOS[c]*σSIN[s])
  else
    GM = (μSIN[s]*μCOS[c]*(N+1)+
          -(1/2)*sin((N+1)*ωs[s])*cos((N+1)*ωc[c])
          +(1/2)*μSIN[s]*cos((N+1)*ωc[c])-(1/2)*μCOS[c]*sin(ωs[s])*cos((N+1)*ωs[s])/(cos(ωs[s])-1)
          +(1/2)*μSIN[s]*sin(ωc[c])*sin((N+1)*ωc[c])/(cos(ωc[c])-1)+(1/2)*sin((N+1)*ωs[s])*μCOS[c]-μSIN[s]*μCOS[c]
          +(1/2)*sin(ωs[s])*cos(ωc[c])
          -(1/2)*μSIN[s]*cos(ωc[c])
          +(1/2)*μCOS[c]*sin(ωs[s])*cos(ωs[s])/(cos(ωs[s])-1)-(1/2)*μSIN[s]*sin(ωc[c])^2/(cos(ωc[c])-1)
          -(1/2)*sin(μSIN[s])*μCOS[c])/(σCOS[c]*σSIN[s])
  end

  return GM
end
function GM45(i::Int, j::Int, d, IT)
    return GM45(i, j, IT.obs, d.μ[SIN], d.σ[SIN], d.μ[COS], d.σ[COS],d.fs,d.fc)
end
function GM54(i::Int, j::Int, d, IT)
    return GM45(j, i, IT.obs, d.μ[SIN], d.σ[SIN], d.μ[COS], d.σ[COS],d.fs,d.fc)
end

# cosine x cosine
function GM55(i::Int, j::Int, N::Int, μCOS::Vector{Float64}, σCOS::Vector{Float64},ω::Vector{Float64})

  if cos(ω[i])-cos(ω[j]) != 0
    GM = (μCOS[j]*μCOS[i]*(N+1)-(1/2)*cos((N+1)*ω[j])*cos((N+1)*ω[i])
          +(1/2)*sin(ω[i])*cos((N+1)*ω[j])*sin((N+1)*ω[i])/(cos(ω[j])-cos(ω[i]))
          -(1/2)*sin(ω[j])*sin((N+1)*ω[j])*cos((N+1)*ω[i])/(cos(ω[j])-cos(ω[i]))
          +(1/2)*μCOS[j]*cos((N+1)*ω[i])+(1/2)*cos((N+1)*ω[j])*μCOS[i]
          +(1/2)*μCOS[j]*sin(ω[i])*sin((N+1)*ω[i])/(-1+cos(ω[i]))
          +(1/2)*μCOS[i]*sin(ω[j])*sin((N+1)*ω[j])/(cos(ω[j])-1)
          -μCOS[j]*μCOS[i]+(1/2)*cos(ω[j])*cos(ω[i])-(1/2)*sin(ω[i])^2*cos(ω[j])/(cos(ω[j])-cos(ω[i]))
          +(1/2)*sin(ω[j])^2*cos(ω[i])/(cos(ω[j])-cos(ω[i]))-(1/2)*μCOS[j]*cos(ω[i])
          -(1/2)*cos(ω[j])*μCOS[i]-(1/2)*μCOS[j]*sin(ω[i])^2/(-1+cos(ω[i]))
          -(1/2)*μCOS[i]*sin(ω[j])^2/(cos(ω[j])-1))/(σCOS[j]*σCOS[i])
  else
    GM = (μCOS[j]*μCOS[i]*(N+1)-(1/2)*cos((N+1)*ω[j])*cos((N+1)*ω[i])
          +(1/2)*μCOS[j]*cos((N+1)*ω[i])+(1/2)*cos((N+1)*ω[j])*μCOS[i]
          +(1/2)*μCOS[j]*sin(ω[i])*sin((N+1)*ω[i])/(-1+cos(ω[i]))
          +(1/2)*μCOS[i]*sin(ω[j])*sin((N+1)*ω[j])/(cos(ω[j])-1)
          -μCOS[j]*μCOS[i]+(1/2)*cos(ω[j])*cos(ω[i])
          -(1/2)*μCOS[j]*cos(ω[i])
          -(1/2)*cos(ω[j])*μCOS[i]-(1/2)*μCOS[j]*sin(ω[i])^2/(-1+cos(ω[i]))
          -(1/2)*μCOS[i]*sin(ω[j])^2/(cos(ω[j])-1))/(σCOS[j]*σCOS[i])
  end
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
