# Closed-form formulas for the inner products between different kinds of components

# step x step
function GMStepStep(m::Int, M::Int, N::Int, μm::Float64, σm::Float64, μM::Float64, σM::Float64)


  return ( m*μm*μM-(M-m)*(1-μm)*μM+(N-M)*(1-μm)*(1-μM) )/(σm*σM) :: Float64
end
function GMStepStep(i::Int, j::Int, d, IT)
  	if i>j
    	return GMStepStep(j, i, IT.obs, d.μ[STEP][j], d.σ[STEP][j], d.μ[STEP][i], d.σ[STEP][i]) :: Float64
    else
    	return GMStepStep(i, j, IT.obs, d.μ[STEP][i], d.σ[STEP][i], d.μ[STEP][j], d.σ[STEP][j]) :: Float64
    end
    #return 0.0 :: Float64
end

# step x slope
function GMStepSlope(
    t::Int, l::Int, N::Int, μSTEP::Float64, σSTEP::Float64,
    μSLOPE::Float64, σSLOPE::Float64
    )

  if l > t
    return (
      t*μSTEP*μSLOPE - (l-t)*(1-μSTEP)*μSLOPE - (N+1)*l-μSLOPE*(N+1)
      + (N+1)*l*μSTEP + μSLOPE*(N+1)*μSTEP + (1/2)*(N+1)^2-(1/2)*N
      - (1/2)*μSTEP*(N+1)^2 + (1/2)*μSTEP*(N+1) + (l+1)*l+μSLOPE*(l+1)
      - (l+1)*l*μSTEP - μSLOPE*(l+1)*μSTEP - (1/2)*(l+1)^2+(1/2)*l
      + (1/2)*μSTEP*(l+1)^2 - (1/2)*μSTEP*(l+1)
      ) :: Float64
  else
    return (
      - (l+1)*l*μSTEP+(1/2)*t - (1/2)*N+(1/2)*(N+1)^2 + μSLOPE*(N+1)*μSTEP
      - μSLOPE*(l+1)*μSTEP - μSLOPE*(N+1) - (1/2)*μSTEP*(N+1)^2
      + (1/2)*μSTEP*(N+1) + (1/2)*μSTEP*(l+1)^2 - (1/2)*μSTEP*(l+1)
      + μSLOPE*(t+1) + l*μSTEP*μSLOPE - (1/2)*(t+1)^2 + (N+1)*l*μSTEP
      - (N+1)*l+(t+1)*l
      ) :: Float64
  end

  #return 0.0 :: Float64
end
function GMStepSlope(i::Int, j::Int, d, IT)
    return GMStepSlope(i, j, IT.obs, d.μ[STEP][i], d.σ[STEP][i], d.μ[SLOPE][j], d.σ[SLOPE][j]) :: Float64
end
function GMSlopeStep(i::Int, j::Int, d, IT)
    return GMStepSlope(j, i, IT.obs, d.μ[STEP][j], d.σ[STEP][j], d.μ[SLOPE][i], d.σ[SLOPE][i]) :: Float64
end

# step x spike
function GMStepSpike(
    i::Int, j::Int, N::Int, μSTEP::Float64, σSTEP::Float64,
    μSPIKE::Float64, σSPIKE::Float64
    )

  if j > i
    return (
      i*μSPIKE*μSTEP - (N-i)*(1-μSTEP)*μSPIKE + (1-μSTEP)*μSPIKE + (1-μSPIKE)*(1-μSTEP)
      )  :: Float64
  else
    return (
      i*μSPIKE*μSTEP - μSTEP*μSPIKE - μSTEP*(1-μSPIKE) - (N-i)*(1-μSTEP)*μSPIKE
      )  :: Float64
  end

  return 0.0 :: Float64
end
function GMStepSpike(i::Int, j::Int, d, IT)
    return GMStepSpike(i, j, IT.obs, d.μ[STEP][i], d.σ[STEP][i], d.μ[SPIKE][j], d.σ[SPIKE][j])  :: Float64
end
function GMSpikeStep(i::Int, j::Int, d, IT)
    return GMStepSpike(j, i, IT.obs, d.μ[STEP][j], d.σ[STEP][j], d.μ[SPIKE][i], d.σ[SPIKE][i])  :: Float64
end

# step x sine
function GMStepSin(
    t::Int, s::Int, N::Int, μSTEP::Float64, σSTEP::Float64,
    μSIN::Float64, σSIN::Float64,ω::Float64
    )

  return (
    μSTEP*μSIN*(t+1)-(1/2)*μSTEP*sin(ω)*cos((t+1)*ω)/(cos(ω)-1)
    +(1/2)*μSTEP*sin((t+1)*ω)-μSTEP*μSIN+(1/2)*μSTEP*sin(ω)*cos(ω)/(cos(ω)-1)
    -(1/2)*μSTEP*sin(ω)+(μSTEP*μSIN-μSIN)*(N+1)-(1/2)*(-1+μSTEP)*sin(ω)*cos((N+1)*ω)/(cos(ω)-1)
    +(-1/2+(1/2)*μSTEP)*sin((N+1)*ω)-(μSTEP*μSIN-μSIN)*(t+1)
    +(1/2)*(-1+μSTEP)*sin(ω)*cos((t+1)*ω)/(cos(ω)-1)-(-1/2+(1/2)*μSTEP)*sin((t+1)*ω)
    ) :: Float64

  #return GM::Float64
end
function GMStepSin(i::Int, j::Int, d, IT)
    return GMStepSin(i, j, IT.obs, d.μ[STEP][i], d.σ[STEP][i], d.μ[SIN][j], d.σ[SIN][j],d.fs[j]) :: Float64
end
function GMSinStep(i::Int, j::Int, d, IT)
    return GMStepSin(j, i, IT.obs, d.μ[STEP][j], d.σ[STEP][j], d.μ[SIN][i], d.σ[SIN][i],d.fs[i]) :: Float64
end

# step x cosine
function GMStepCos(
    t::Int, c::Int, N::Int, μSTEP::Float64, σSTEP::Float64,
    μCOS::Float64, σCOS::Float64,ω::Float64
    )

    return (
    μSTEP*μCOS*(t+1)+(1/2)*μSTEP*cos((t+1)*ω)
    +(1/2)*μSTEP*sin(ω)*sin((t+1)*ω)/(cos(ω)-1)
    -μSTEP*μCOS-(1/2)*μSTEP*cos(ω)-(1/2)*μSTEP*sin(ω)^2/(cos(ω)-1)
    +(μSTEP*μCOS-μCOS)*(N+1)+(-1/2+(1/2)*μSTEP)*cos((N+1)*ω)
    -(1/2)*(1-μSTEP)*sin(ω)*sin((N+1)*ω)/(cos(ω)-1)
    -(μSTEP*μCOS-μCOS)*(t+1)-(-1/2+(1/2)*μSTEP)*cos((t+1)*ω)
    +(1/2)*(1-μSTEP)*sin(ω)*sin((t+1)*ω)/(cos(ω)-1)
    ) :: Float64
end
function GMStepCos(i::Int, j::Int, d, IT)
    return GMStepCos(i, j, IT.obs, d.μ[STEP][i], d.σ[STEP][i], d.μ[COS][j], d.σ[COS][j],d.fc[j]) :: Float64
end
function GMCosStep(i::Int, j::Int, d, IT)
    return GMStepCos(j, i, IT.obs, d.μ[STEP][j], d.σ[STEP][j], d.μ[COS][i], d.σ[COS][i],d.fc[i]) :: Float64
end

# slope x slope
function GMSlopeSlope(m::Int, M::Int, N::Int, μm::Float64, σm::Float64, μM::Float64, σM::Float64)

  return (
    (1/6)*N - (1/6)*M + m*μM*μM - μM*(m+1)*μM - μM*(m+1)*m
    - (1/2)*(N+1)^2 + (1/3)*(N+1)^3 + (1/2)*(M+1)^2 - (1/3)*(M+1)^3
    + (1/2)*μM*(m+1)^2 - (1/2)*μM*(m+1) + (1/2)*m*(M+1)^2
    - (1/2)*m*(M+1) + (1/2)μM*(M+1)^2 - (1/2)*μM*(M+1)
    + (1/2)*M*(M+1)^2 - (1/2)*M*(M+1) - (1/2)*m*(N+1)^2 + (1/2)*m*(N+1)
    - (1/2)*μM*(N+1)^2 + (1/2)*μM*(N+1) - (1/2)*M*(N+1)^2
    + (1/2)*M*(N+1) - (1/2)*μM*(N+1)^2 - (M+1)*m*M
    - (M+1)*μM*M + (N+1)*m*M + μM*(N+1)*m + (N+1)*μM*M
    + μM*(N+1)*μM + (1/2)*μM*(N+1)
    ) :: Float64
end
function GMSlopeSlope(i::Int, j::Int, d, IT)
  	if j>i
    	return GMSlopeSlope(i, j, IT.obs, d.μ[SLOPE][i], d.σ[SLOPE][i], d.μ[SLOPE][j], d.σ[SLOPE][j]) :: Float64
	else
		return GMSlopeSlope(j, i, IT.obs, d.μ[SLOPE][j], d.σ[SLOPE][j], d.μ[SLOPE][i], d.σ[SLOPE][i]) :: Float64
	end
	return 0.0::Float64
end

# slope x spike
function GMSlopeSpike(
    l::Int, p::Int, N::Int, μSLOPE::Float64, σSLOPE::Float64,
                            μSPIKE::Float64, σSPIKE::Float64,
    )

  if p > l
    return (
      l*μSLOPE*μSPIKE + μSPIKE*(N+1)*l + μSPIKE*(N+1)*μSLOPE
      - (1/2)*μSPIKE*(N+1)^2 + (1/2)*μSPIKE*(N+1) - μSPIKE*(l+1)*l
      - μSPIKE*(l+1)*μSLOPE + (1/2)*μSPIKE*(l+1)^2
      - (1/2)*μSPIKE*(l+1) + (p-l-μSLOPE)*μSPIKE
      + (p-l-μSLOPE)*(1-μSPIKE)
      ) :: Float64
  else
    return (
      l*μSLOPE*μSPIKE + μSPIKE*(N+1)*l + μSPIKE*(N+1)*μSLOPE
      - (1/2)*μSPIKE*(N+1)^2 + (1/2)*μSPIKE*(N+1) - μSPIKE*(l+1)*l
      - μSPIKE*(l+1)*μSLOPE + (1/2)*μSPIKE*(l+1)^2
      - (1/2)*μSPIKE*(l+1) - μSLOPE*μSPIKE - μSLOPE*(1-μSPIKE)
      ) :: Float64
  end

  return 0.0::Float64
end
function GMSlopeSpike(i::Int, j::Int, d, IT)
    return GMSlopeSpike(i, j, IT.obs, d.μ[SLOPE][i], d.σ[SLOPE][i], d.μ[SPIKE][j], d.σ[SPIKE][j]) :: Float64
end
function GMSpikeSlope(i::Int, j::Int, d, IT)
    return GMSlopeSpike(j, i, IT.obs, d.μ[SLOPE][j], d.σ[SLOPE][j], d.μ[SPIKE][i], d.σ[SPIKE][i]) :: Float64
end

# slope x sine
function GMSlopeSin(
    l::Int, s::Int, N::Int, μSLOPE::Float64, σSLOPE::Float64,
    μSIN::Float64, σSIN::Float64,ω::Float64
    )

  return (
    μSLOPE*μSIN*(l+1)-(1/2)*μSLOPE*sin(ω)*cos((l+1)*ω)/(cos(ω)-1)
    +(1/2)*μSLOPE*sin((l+1)*ω)-μSLOPE*μSIN+(1/2)*μSLOPE*sin(ω)*cos(ω)/(cos(ω)-1)
    -(1/2)*μSLOPE*sin(ω)+(l*μSIN+(1/2)*μSIN+μSLOPE*μSIN)*(N+1)
    -(1/2)*μSIN*(N+1)^2+(1/2)*sin(ω)*(N+1)*cos((N+1)*ω)/(cos(ω)-1)-(1/2*(N+1))*sin((N+1)*ω)
    -(1/2)*(μSLOPE+l)*sin(ω)*cos((N+1)*ω)/(cos(ω)-1)
    +(1/2)*(l*cos(ω)+μSLOPE*cos(ω)-l-μSLOPE-1)*sin((N+1)*ω)/(cos(ω)-1)
`  -(l*μSIN+(1/2)*μSIN+μSLOPE*μSIN)*(l+1)
    +(1/2)*μSIN*(l+1)^2-(1/2)*sin(ω)*(l+1)*cos((l+1)*ω)/(cos(ω)-1)
    +(1/2*(l+1))*sin((l+1)*ω)
    +(1/2)*(μSLOPE+l)*sin(ω)*cos((l+1)*ω)/(cos(ω)-1)
    -(1/2)*(l*cos(ω)+μSLOPE*cos(ω)-l-μSLOPE-1)*sin((l+1)*ω)/(cos(ω)-1)
    ) :: Float64
end
function GMSlopeSin(i::Int, j::Int, d, IT)
    return GMSlopeSin(i, j, IT.obs, d.μ[SLOPE][i], d.σ[SLOPE][i], d.μ[SIN][j], d.σ[SIN][j], d.fs[j]) :: Float64
end
function GMSinSlope(i::Int, j::Int, d, IT)
    return GMSlopeSin(j, i, IT.obs, d.μ[SLOPE][j], d.σ[SLOPE][j], d.μ[SIN][i], d.σ[SIN][i], d.fs[i]) :: Float64
end

# slope x cosine
function GMSlopeCos(
    l::Int, c::Int, N::Int, μSLOPE::Float64, σSLOPE::Float64,
    μCOS::Float64, σCOS::Float64, ω::Float64
    )

  return (
    μSLOPE*μCOS*(l+1)+(1/2)*μSLOPE*cos((l+1)*ω)
    +(1/2)*μSLOPE*sin(ω)*sin((l+1)*ω)/(cos(ω)-1)
    -μSLOPE*μCOS-(1/2)*μSLOPE*cos(ω)-(1/2)*μSLOPE*sin(ω)^2/(cos(ω)-1)
    +(μSLOPE*μCOS+l*μCOS+(1/2)*μCOS)*(N+1)
    -(1/2)*μCOS*(N+1)^2-(1/2*(N+1))*cos((N+1)*ω)
    -(1/2)*sin(ω)*(N+1)*sin((N+1)*ω)/(cos(ω)-1)
    +(1/2)*(l*cos(ω)+μSLOPE*cos(ω)-l-μSLOPE-1)*cos((N+1)*ω)/(cos(ω)-1)
    +(1/2)*(l+μSLOPE)*sin(ω)*sin((N+1)*ω)/(cos(ω)-1)
    -(μSLOPE*μCOS+l*μCOS+(1/2)*μCOS)*(l+1)
    +(1/2)*μCOS*(l+1)^2+(1/2*(l+1))*cos((l+1)*ω)
    +(1/2)*sin(ω)*(l+1)*sin((l+1)*ω)/(cos(ω)-1)
    -(1/2)*(l*cos(ω)+μSLOPE*cos(ω)-l-μSLOPE-1)*cos((l+1)*ω)/(cos(ω)-1)
    -(1/2)*(l+μSLOPE)*sin(ω)*sin((l+1)*ω)/(cos(ω)-1)
    ) :: Float64
end
function GMSlopeCos(i::Int, j::Int, d, IT)
    return GMSlopeCos(i, j, IT.obs, d.μ[SLOPE][i], d.σ[SLOPE][i], d.μ[COS][j], d.σ[COS][j],d.fc[j]) :: Float64
end
function GMCosSlope(i::Int, j::Int, d, IT)
    return GMSlopeCos(j, i, IT.obs, d.μ[SLOPE][j], d.σ[SLOPE][j], d.μ[COS][i], d.σ[COS][i],d.fc[i]) :: Float64
end

# spike x spike
function GMSpikeSpike(i::Int, j::Int, N::Int, μSPIKEi::Float64, σSPIKEi::Float64)

  if i != j
    return (
      (N-2)*μSPIKEi^2 - 2*μSPIKEi*(1-μSPIKEi)
    ) :: Float64
  else
    return float(N) :: Float64
  end

  return 0.0 :: Float64
end
function GMSpikeSpike(i::Int, j::Int, d, IT)
    return GMSpikeSpike(i, j, IT.obs, d.μ[SPIKE][i], d.σ[SPIKE][i]) :: Float64
end

# spike x sine
function GMSpikeSin(
    p::Int, s::Int, N::Int, μSPIKE::Float64, σSPIKE::Float64,
    μSIN::Float64, σSIN::Float64, ω::Float64
    )

   return (μSPIKE*μSIN*(N+1)-(1/2)*μSPIKE*sin(ω)*cos((N+1)*ω)/(cos(ω)-1)
        +(1/2)*μSPIKE*sin((N+1)*ω)-μSPIKE*μSIN+(1/2)*μSPIKE*sin(ω)*cos(f[s])/(cos(f[s])-1)
        -(1/2)*μSPIKE*sin(ω)+μSPIKE*(sin(p*ω)-μSIN)+(1-μSPIKE)*(sin(p*ω)-μSIN)
        ) :: Float64
end
function GMSpikeSin(i::Int, j::Int, d, IT)
    return GMSpikeSin(i, j, IT.obs, d.μ[SPIKE][i], d.σ[SPIKE][i], d.μ[SIN][j], d.σ[SIN][j], d.fs[j]) :: Float64
end
function GMSinSpike(i::Int, j::Int, d, IT)
    return GMSpikeSin(j, i, IT.obs, d.μ[SPIKE][j], d.σ[SPIKE][j], d.μ[SIN][i], d.σ[SIN][i], d.fs[i]) :: Float64
end

# spike x cosine
function GMSpikeCos(
    p::Int, c::Int, N::Int, μSPIKE::Float64, σSPIKE::Float64,
    μCOS::Float64, σCOS::Float64, ω::Float64
    )

  return (
    μSPIKE*μCOS*(N+1)+(1/2)*μSPIKE*cos((N+1)*ω)
    +(1/2)*μSPIKE*sin(ω)*sin((N+1)*ω)/(cos(ω)-1)
    -μSPIKE*μCOS-(1/2)*μSPIKE*cos(ω)-(1/2)*μSPIKE*sin(ω)^2/(cos(ω)-1)
    +μSPIKE*(cos(p*ω)-μCOS)+(1-μSPIKE)*(cos(p*ω)-μCOS)
    ) :: Float64
end
function GMSpikeCos(i::Int, j::Int, d, IT)
    return GMSpikeCos(i, j, IT.obs, d.μ[SPIKE][i], d.σ[SPIKE][i], d.μ[COS][j], d.σ[COS][j],d.fc[j]) :: Float64
end
function GMCosSpike(i::Int, j::Int, d, IT)
    return GMSpikeCos(j, i, IT.obs, d.μ[SPIKE][j], d.σ[SPIKE][j], d.μ[COS][i], d.σ[COS][i],d.fc[i]) :: Float64
end

# sine x sine
function GMSinSin(i::Int, j::Int, N::Int, μSINi::Float64, σSINi::Float64, μSINj::Float64, σSINj::Float64,ωi::Float64,ωj::Float64)

  #ATENTION HERE
  #there still may be a ptoblem from cos(ωi) = cos(ωj) with i != j
  if i == j
    return float(N) :: Float64
  elseif abs( cos(ωj)-cos(ωi) ) <= 1e-8
    #cancelation that need to happen for finiteness
    return (
      μSINj*μSINi*(N+1)
      -(1/2)*sin((N+1)*ωj)*sin((N+1)*ωi)-(1/2)*μSINj*sin(ωi)*cos((N+1)*ωi)/(-1+cos(ωi))
      -(1/2)*μSINi*sin(ωj)*cos((N+1)*ωj)/(cos(ωj)-1)+(1/2)*μSINj*sin((N+1)*ωi)
      +(1/2)*sin((N+1)*ωj)*μSINi-μSINj*μSINi
      +(1/2)*sin(ωi)*sin(ωj)
      +(1/2)*μSINj*sin(ωi)*cos(ωi)/(-1+cos(ωi))+(1/2)*μSINi*sin(ωj)*cos(ωj)/(cos(ωj)-1)
      -(1/2)*μSINj*sin(ωi)-(1/2)*μSINi*sin(ωj)
      ) :: Float64
  else
    return (
      μSINj*μSINi*(N+1)+(1/2)*sin(ωj)*cos((N+1)*ωj)*sin((N+1)*ωi)/(cos(ωj)-cos(ωi))
      -(1/2)*sin(ωi)*sin((N+1)*ωj)*cos((N+1)*ωi)/(cos(ωj)-cos(ωi))
      -(1/2)*sin((N+1)*ωj)*sin((N+1)*ωi)-(1/2)*μSINj*sin(ωi)*cos((N+1)*ωi)/(-1+cos(ωi))
      -(1/2)*μSINi*sin(ωj)*cos((N+1)*ωj)/(cos(ωj)-1)+(1/2)*μSINj*sin((N+1)*ωi)
      +(1/2)*sin((N+1)*ωj)*μSINi-μSINj*μSINi-(1/2)*sin(ωj)*cos(ωj)*sin(ωi)/(cos(ωj)-cos(ωi))
      +(1/2)*sin(ωi)*sin(ωj)*cos(ωi)/(cos(ωj)-cos(ωi))+(1/2)*sin(ωi)*sin(ωj)
      +(1/2)*μSINj*sin(ωi)*cos(ωi)/(-1+cos(ωi))+(1/2)*μSINi*sin(ωj)*cos(ωj)/(cos(ωj)-1)
      -(1/2)*μSINj*sin(ωi)-(1/2)*μSINi*sin(ωj)
      ) :: Float64
  end

  return 0.0 :: Float64
end
function GMSinSin(i::Int, j::Int, d, IT)
    return GMSinSin(i, j, IT.obs, d.μ[SIN][i], d.σ[SIN][i], d.μ[SIN][j], d.σ[SIN][j], d.fs[i], d.fs[j]) :: Float64
end


# sine x cosine
function GMSinCos(
    s::Int, c::Int, N::Int, μSIN::Float64, σSIN::Float64,
    μCOS::Float64, σCOS::Float64, ωs::Float64, ωc::Float64
    )

  if abs(cos(ωc)-cos(ωs)) >= 1e-5

    return (
      μSIN*μCOS*(N+1)+(1/2)*sin(ωs)*cos((N+1)*ωs)*cos((N+1)*ωc)/(cos(ωs)-cos(ωc))
      -(1/2)*sin((N+1)*ωs)*cos((N+1)*ωc)+(1/2)*sin(ωc)*sin((N+1)*ωs)*sin((N+1)*ωc)/(cos(ωs)-cos(ωc))
      +(1/2)*μSIN*cos((N+1)*ωc)-(1/2)*μCOS*sin(ωs)*cos((N+1)*ωs)/(cos(ωs)-1)
      +(1/2)*μSIN*sin(ωc)*sin((N+1)*ωc)/(cos(ωc)-1)+(1/2)*sin((N+1)*ωs)*μCOS-μSIN*μCOS
      -(1/2)*sin(ωs)*cos(ωs)*cos(ωc)/(cos(ωs)-cos(ωc))+(1/2)*sin(ωs)*cos(ωc)
      -(1/2)*sin(ωc)^2*sin(ωs)/(cos(ωs)-cos(ωc))-(1/2)*μSIN*cos(ωc)+(1/2)*μCOS*sin(ωs)*cos(ωs)/(cos(ωs)-1)
      -(1/2)*μSIN*sin(ωc)^2/(cos(ωc)-1)-(1/2)*sin(ωs)*μCOS
      ) :: Float64

  else
    return (
      μSIN*μCOS*(N+1)-(1/2)*cos(ωs)*cos((N+1)*ωs)^2/sin(ωs)-(1/2)*sin((N+1)*ωs)*cos((N+1)*ωs)
      +(1/2)*(cos(ωs)*μCOS+μSIN*sin(ωs)+μCOS)*cos((N+1)*ωs)/sin(ωs)
      +(1/2)*(-cos(ωs)*μSIN+μCOS*sin(ωs)-μSIN)*sin((N+1)*ωs)/sin(ωs)-μSIN*μCOS
      +(1/2)*cos(ωs)^3/sin(ωs)+(1/2)*cos(ωs)*sin(ωs)-(1/2)*(cos(ωs)*μCOS+μSIN*sin(ωs)+μCOS)*cos(ωs)/sin(ωs)
      +(1/2)*cos(ωs)*μSIN-(1/2)*μCOS*sin(ωs)+(1/2)*μSIN
      ) :: Float64

  end

  return 0.0 :: Float64
end
function GMSinCos(i::Int, j::Int, d, IT)
    return GMSinCos(i, j, IT.obs, d.μ[SIN][i], d.σ[SIN][i], d.μ[COS][j], d.σ[COS][j],d.fs[i],d.fc[j]) :: Float64
end
function GMCosSin(i::Int, j::Int, d, IT)
    return GMSinCos(j, i, IT.obs, d.μ[SIN][j], d.σ[SIN][j], d.μ[COS][i], d.σ[COS][i],d.fs[j],d.fc[i]) :: Float64
end

# cosine x cosine
function GMCosCos(i::Int, j::Int, N::Int, μCOSi::Float64, σCOSi::Float64, μCOSj::Float64, σCOSj::Float64,ωi::Float64,ωj::Float64)

  if i == j
    return float(N)
  elseif abs(cos(ωi)-cos(ωj) ) >= 1e-8
    return (
      μCOSj*μCOSi*(N+1)-(1/2)*cos((N+1)*ωj)*cos((N+1)*ωi)
      +(1/2)*sin(ωi)*cos((N+1)*ωj)*sin((N+1)*ωi)/(cos(ωj)-cos(ωi))
      -(1/2)*sin(ωj)*sin((N+1)*ωj)*cos((N+1)*ωi)/(cos(ωj)-cos(ωi))
      +(1/2)*μCOSj*cos((N+1)*ωi)+(1/2)*cos((N+1)*ωj)*μCOSi
      +(1/2)*μCOSj*sin(ωi)*sin((N+1)*ωi)/(-1+cos(ωi))
      +(1/2)*μCOSi*sin(ωj)*sin((N+1)*ωj)/(cos(ωj)-1)-μCOSj*μCOSi+(1/2)*cos(ωj)*cos(ωi)
      -(1/2)*sin(ωi)^2*cos(ωj)/(cos(ωj)-cos(ωi))
      +(1/2)*sin(ωj)^2*cos(ωi)/(cos(ωj)-cos(ωi))-(1/2)*μCOSj*cos(ωi)
      -(1/2)*cos(ωj)*μCOSi-(1/2)*μCOSj*sin(ωi)^2/(-1+cos(ωi))-(1/2)*μCOSi*sin(ωj)^2/(cos(ωj)-1)
      ) :: Float64
  else
    return (
      μCOSj*μCOSi*(N+1)-(1/2)*cos((N+1)*ωj)*cos((N+1)*ωi)
      +(1/2)*μCOSj*cos((N+1)*ωi)+(1/2)*cos((N+1)*ωj)*μCOSi
      +(1/2)*μCOSj*sin(ωi)*sin((N+1)*ωi)/(-1+cos(ωi))
      +(1/2)*μCOSi*sin(ωj)*sin((N+1)*ωj)/(cos(ωj)-1)
      -μCOSj*μCOSi+(1/2)*cos(ωj)*cos(ωi)
      -(1/2)*μCOSj*cos(ωi)-(1/2)*cos(ωj)*μCOSi
      -(1/2)*μCOSj*sin(ωi)^2/(-1+cos(ωi))
      -(1/2)*μCOSi*sin(ωj)^2/(cos(ωj)-1)
      ) :: Float64
  end
  return 0.0::Float64
end
function GMCosCos(i::Int, j::Int, d, IT)
    return GMCosCos(i, j, IT.obs, d.μ[COS][i], d.σ[COS][i],d.μ[COS][j], d.σ[COS][j],d.fc[i],d.fc[j]) :: Float64
end

# GM = Matrix{Function}(TOTALCOMPONENTS, TOTALCOMPONENTS)
# GM[STEP, STEP]   = GMStepStep
# GM[STEP, SPIKE]  = GMStepSpike
# GM[SPIKE, STEP]  = GMSpikeStep
# GM[STEP, SLOPE]  = GMStepSlope
# GM[SLOPE, STEP]  = GMSlopeStep
# GM[STEP, SIN]    = GMStepSin
# GM[SIN, STEP]    = GMSinStep
# GM[STEP, COS]    = GMStepCos
# GM[COS, STEP]    = GMCosStep
# GM[SPIKE, SPIKE] = GMSpikeSpike
# GM[SPIKE, SLOPE] = GMSpikeSlope
# GM[SLOPE, SPIKE] = GMSlopeSpike
# GM[SPIKE, SIN]   = GMSpikeSin
# GM[SIN, SPIKE]   = GMSinSpike
# GM[SPIKE, COS]   = GMSpikeCos
# GM[COS, SPIKE]   = GMCosSpike
# GM[SLOPE, SLOPE] = GMSlopeSlope
# GM[SLOPE, SIN]   = GMSlopeSin
# GM[SIN, SLOPE]   = GMSinSlope
# GM[SLOPE, COS]   = GMSlopeCos
# GM[COS, SLOPE]   = GMCosSlope
# GM[SIN, SIN]     = GMSinSin
# GM[SIN, COS]     = GMSinCos
# GM[COS, SIN]     = GMCosSin
# GM[COS, COS]     = GMCosCos

function GM2(c1,c2,i::Int, j::Int, d, IT)
  if c1 == STEP && c2 == STEP
    return GMStepStep(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == STEP && c2 == SPIKE
    return GMStepSpike(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SPIKE && c2 == STEP
    return GMSpikeStep(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == STEP && c2 == SLOPE
    return GMStepSlope(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SLOPE && c2 == STEP
    return GMSlopeStep(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == STEP && c2 == SIN
    return GMStepSin(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SIN && c2 == STEP
    return GMSinStep(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == STEP && c2 == COS
    return GMStepCos(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == COS && c2 == STEP
    return GMCosStep(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SPIKE && c2 == SPIKE
    return GMSpikeSpike(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SPIKE && c2 == SLOPE
    return GMSpikeSlope(i::Int, j::Int, d, IT):: Float64
  elseif c1 == SLOPE && c2 == SPIKE
    return GMSlopeSpike(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SPIKE && c2 == SIN
    return GMSpikeSin(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SIN && c2 == SPIKE
    return GMSinSpike(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SPIKE && c2 == COS
    return GMSpikeCos(i::Int, j::Int, d, IT):: Float64
  elseif c1 == COS && c2 == SPIKE
    return GMCosSpike(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SLOPE && c2 == SLOPE
    return GMSlopeSlope(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SLOPE && c2 == SIN
    return GMSlopeSin(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SIN && c2 == SLOPE
    return GMSinSlope(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SLOPE && c2 == COS
    return GMSlopeCos(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == COS && c2 == SLOPE
    return GMCosSlope(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SIN && c2 == SIN
    return GMSinSin(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == SIN && c2 == COS
    return GMSinCos(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == COS && c2 == SIN
    return GMCosSin(i::Int, j::Int, d, IT) :: Float64
  elseif c1 == COS && c2 == COS
    return GMCosCos(i::Int, j::Int, d, IT) :: Float64
  else
    return 0.0 :: Float64
  end
end





