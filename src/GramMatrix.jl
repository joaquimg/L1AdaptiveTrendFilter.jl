#ATENCAO
# considerar mu e sigma como input de cada função
# as outras entradas são i, j (numeros dos elementos)
# pode usar o N como numero total de elementos.
# eu estou usando uma struct chamada iterator (ver em CDtype.jl)
# posso corrigir isso sem problemas
# o principla é escrever as fómulas mesmo

#to complete

GM = Matrix{Function}(2,2)
GM[1,2]=GM12
GM[2,1]=GM12

function GM11(i::Int,j::Int,L::Vector{Float64},U::Vector{Float64},N::Int,
              μ0::Float64,σ0::Float64)
  # step x step
  GM = Float64[]
  m = min(i, j)
  M = max(i, j)

  GM = (
    m*μ[m]*μ[M]-(M-m)*(1-μ[m])*μ[M]+(N-M)*(1-μ[m])*(1-μ[M])
    )/(σ[m]*σ[M])

  return GM
end

function GM12(i::Int,j::Int,N::Int,μp::Vector{Float64},σp::Vector{Float64},
              μs::Vector{Float64},σs::Vector{Float64})
  # step x slope
  GM = Float64[]

  if i > j
    GM = (
      j*μ[j]*μ[i]-(i-j)*(1-μ[j])*μ[i]-μ[i]*(N+1)+μ[i]*μ[j]*(N+1)
      + 0.5*(N+1)^2-0.5*N-0.5*μ[j]*(N+1)^2+0.5*μ[j]*(N+1)
      + μ[i]*(i+1)-μ[i]*μ[j]*(i+1)-0.5*(i+1)^2+0.5*i
      + 0.5*μ[j]*(i+1)^2-0.5*μ[j]*(i+1)
      )/(σ[i]*σ[j])
  else
    GM = (
      i*μ[i]*μ[j]-(i+1)*μ[i]*μ[j]+0.5*μ[j]*(i+1)^2-0.5*μ[j]*(i+1)
      - μ[i]*(N+1)+μ[i]*μ[j]*(N+1)+0.5*(N+1)^2-0.5*N-0.5*μ[j]*(N+1)^2
      + 0.5*μ[j]*(N+1)+μ[i]*(j+1)-0.5*(j+1)^2+0.5*j
      )/(σ[i]*σ[j])
  end

  return GM
end

function GM13(i::Int,j::Int,N::Int,μp::Vector{Float64},σp::Vector{Float64},
              μs::Vector{Float64},σs::Vector{Float64})
  # step x spike
  GM = Float64[]

  if i > j
    GM = (
      j*μ[j]*μ[i]-(N-j)*(1-μ[j])*μ[i]+(1-μ[j])*μ[i]+(1-μ[j])*(1-μ[i])
      )/(σ[i]*σ[j])
  else
    GM = (
      j*μ[j]*μ[i]-μ[i]*μ[j]-μ[j]*(1-μ[i])-(N-j)*(1-μ[j])*μ[i]
      )/(σ[i]*σ[j])
  end

  return GM
end

function GM22(i::Int,j::Int,N::Int,μ::Vector{Float64},σ::Vector{Float64})


    GM = Float64[]
    ab = convert(Float64,abs(i-j))
    m = convert(Float64,minimum(i,j))
    M = convert(Float64,maximum(i,j))

    GM = (μ[i]*μ[j]*(m-1.0) + μ[i]*μ[j]*ab
            +(-μ[M])*ab*(ab+1.0)/2.0
            +(N-M+1.0)*(-μ[m]+ab)*(-μ[M])
            +(-μ[i]-μ[j]+ab)*(N-M+1.0)*(N-M+2.0)/2.0
            +(N-M+1.0)*(N-M+2.0)*(2.0N-2.0M+3.0)/6.0
            )/σ[i]/σ[j]
    return GM
end



