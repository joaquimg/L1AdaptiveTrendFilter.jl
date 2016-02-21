# Computes μ(mean) and σ(standard deviation) for each kind of component

function getStepData(IT,f)
    getStepData(IT)
end
function getStepData(IT)
    N = IT.obs
    Nel = IT.nelements[STEP]
    μ = Vector{Float64}(Nel)
    σ = Vector{Float64}(Nel)

    for i = 1:Nel
        μ[i] = (N-i)/N
        σ[i] = sqrt(( i*μ[i]^2+(N-i)*(1-μ[i])^2 )/N)
    end
    return σ,μ
end

function getSlopeData(IT,f)
    getSlopeData(IT)
end
function getSlopeData(IT)
    N = IT.obs
    Nel = IT.nelements[SLOPE]
    μ = Vector{Float64}(Nel)
    σ = Vector{Float64}(Nel)

    for i = 1:Nel
        μ[i] = (
      0.5*(N+1)^2 - 0.5*N - i*(N+1) - 0.5*(i+1)^2 + 0.5*i + i*(i+1)
      )/N

        σ[i] = sqrt((
      i*μ[i]^2 + (N+1)*i^2 + 2*(N+1)*i*μ[i] + (N+1)*μ[i]^2 - i*(N+1)^2
      + i*(N+1) - μ[i]*(N+1)^2 + μ[i]*(N+1) + (1/3)*(N+1)^3
      - 0.5*(N+1)^2 + (1/6)*N - (i+1)*i^2 - 2*(i+1)*i*μ[i]
      - (i+1)*μ[i]^2 + i*(i+1)^2 - i*(i+1) + μ[i]*(i+1)^2
      - μ[i]*(i+1) - (1/3)*(i+1)^3 + 0.5*(i+1)^2 - (1/6)*i
      )/N)
    end
    return σ,μ
end

function getSpikeData(IT,f)
    getSpikeData(IT)
end
function getSpikeData(IT)
    N = IT.obs
    Nel = IT.nelements[SPIKE]
    μ = Vector{Float64}(Nel)
    σ = Vector{Float64}(Nel)

    for i = 1:Nel
        μ[i] = 1/N
        σ[i] = sqrt( ((N-1)*μ[i]^2 + (1-μ[i])^2)/N )
    end
    return σ,μ
end

function getSineData(IT,f)
    getSineData(IT, f)
end
function getSineData(IT, ω)
    N = IT.obs
    N_frequencies = IT.nelements[4]
    μ = Vector{Float64}(N_frequencies)
    σ = Vector{Float64}(N_frequencies)

    for i = 1:N_frequencies
        sin_ω = sin(ω[i])
        cos_ω = cos(ω[i])
        sin_N1ω = sin((N+1)*ω[i])
        cos_N1ω = cos((N+1)*ω[i])

        μ[i] =  (
      (sin_ω*cos_N1ω)/(cos_ω-1) - sin_N1ω - (sin_ω*cos_ω)/(cos_ω-1) + sin_ω
      )/(2*N)

        σ[i] = sqrt( (
      - μ[i]*sin_ω + (0.5+μ[i]^2)*(N+1) + 0.5*cos_N1ω^2
      - 0.5*(cos_ω*sin_N1ω*cos_N1ω)/sin_ω - (μ[i]*sin_ω*cos_N1ω)/(cos_ω-1)
      + sin_N1ω*μ[i] - 0.5 - μ[i]^2 + (μ[i]*sin_ω*cos_ω)/(cos_ω-1)
      )/N )
    end
    return σ,μ
end

function getCosData(IT,f)
    getCosData(IT, f)
end
function getCosData(IT, ω)
    N = IT.obs
    N_frequencies = IT.nelements[5]
    μ = Vector{Float64}(N_frequencies)
    σ = Vector{Float64}(N_frequencies)

    for i = 1:N_frequencies
        sin_ω = sin(ω[i])
        cos_ω = cos(ω[i])
        sin_N1ω = sin((N+1)*ω[i])
        cos_N1ω = cos((N+1)*ω[i])

        μ[i] =  (
      - 0.5*cos_N1ω - 0.5*(sin_ω*sin_N1ω)/(cos_ω-1) + 0.5*cos_ω
      + 0.5*sin_ω^2/(cos_ω-1)
      )/N

        σ[i] = sqrt( (
      (0.5+μ[i]^2)*(N+1) - 0.5*cos_N1ω^2 - (μ[i]*sin_ω^2)/(cos_ω-1)
      + 0.5*(cos_ω*sin_N1ω*cos_N1ω)/sin_ω + cos_N1ω*μ[i]
      + (μ[i]*sin_ω*sin_N1ω)/(cos_ω-1) -0.5 -μ[i]^2 - cos_ω*μ[i]
      )/N )
    end
    return σ,μ
end

getData = Vector{Function}(5)
getData[STEP] = getStepData
getData[SLOPE] = getSlopeData
getData[SPIKE] = getSpikeData
getData[SIN] = getSineData
getData[COS] = getCosData
