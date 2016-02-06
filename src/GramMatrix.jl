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
function GM12(i::Int,j::Int,N::Int,μp::Vector{Float64},σp::Vector{Float64},
    μs::Vector{Float64},σs::Vector{Float64})
    
    GM = Float64[]
    
    if i > j
        GM = (μp[i]*μs[j]*((i-1.0)+(j-i))
        + (-μs[j]*(j-i)*(j-i+1.0)/2.0)
        + (1.0-μs[j])*(j-i-μp[i])*(N-j+1.0)
        + (1.0-μs[j])*(N-j+1.0)*(N-j+2.0)/2.0
        )/σp[i]/σs[j]
    else
        GM = (μp[i]*μs[j]*((j-1.0)+(i-j))
        + (-μs[j]*(1.0-μs[j])*(i-j)/2.0)
        + (1.0-μs[j]+(i-j))*(-μp[i])*(N-i+1.0)
        + (1.0-μs[j]+(i-j))*(N-i+1.0)*(N-i+2.0)/2.0
        )/σp[i]/σs[j]
    end
    
    return GM
end

function GM12T(i,j,N,μp,σp,μs,σs) = GM12(j,i,N,μp,σp,μs,σs)

function GM11(i::Int,j::Int,L::Vector{Float64},U::Vector{Float64},N::Int,
              μ0::Float64,σ0::Float64)
    
    GM = Float64[]
    
    if i != 0 && j != 0
        
        GM = min(i,j)*U[i]*U[j]+(max(i,j)-min(i,j))*L[i]*U[j]+(N-max(i,j))*L[i]*L[j]
        
    elseif ( i + j == 1 )
        
        GM = ( j*U[j]*((j+1.0)/2.0 - μ0)*(N-j)*L[j]*((j+1.0+N)/2.0 - μ0)  )/σ0
    
    else
        
        GM = N - 1.0
        
    end
    
    return GM
    
end