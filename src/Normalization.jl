

function getSlopeData(N::Int)

    μ = Vector{Float64}(N-1)
    σ = Vector{Float64}(N-1)
    
    #revisar!!!!
    for i = 1:(N)
        μ[i] = (N-i+1.0)*(N-i+2.0)/2.0
        σ[i] = sqrt( μ[i]*(2.0N-2.0i+3.0)/(3.0N) - μ[i]^2.0 )

    end

    return σ,μ

end

function getStepData(N::Int)

    μ = Vector{Float64}(N-1)
    σ = Vector{Float64}(N-1)
    L = Vector{Float64}(N-1)
    U = Vector{Float64}(N-1)

    for i = 1:(N-1)
        μ[i] = (N-i)/N
        σ[i] = sqrt( ( i*μ[i]^2.0+(N-i)*(1.0-μ[i])^2.0 )/N )
        U[i] = -μ[i]/σ[i]
        L[i] = (1.0-μ[i])/σ[i]
    end

    return U,L,σ#,mu

end