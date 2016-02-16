

function getCosData(IT,f)
    N = IT.obs

    nf = IT.nelements[5]

    μ = Vector{Float64}(nf)
    σ = Vector{Float64}(nf)
    sf=0.0
    cf=0.0
    sn1f=0.0
    cn1f=0.0


    for i = 1:nf
        sf = sin(f[i])
        cf = cos(f[i])
        sn1f = sin((N+1)*f[i])
        cn1f = cos((N+1)*f[i])


        μ[i] =  -1/2*cn1f
                -1/2*(sf*sn1f/(cf-1.0))
                +1/2*cf
                +1/2*sf^2/(cf-1.0)

        σ[i] = (1/2 + μ[i]^2)*(N+1) + 1/2*cn1f^2
               -1/2*cf*sn1f*cn1f/sf - μ[i]*sf*cn1f/(cf-1.0)
               +sn1f*μ[i] - 1/2 - μ[i]^2 +μ[i]*sf*cf/(cf-1.0) - μ[i]*sf

    end

    return σ,μ
end


function getSineData(IT,f)
    N = IT.obs

    nf = IT.nelements[4]


    μ = Vector{Float64}(nf)
    σ = Vector{Float64}(nf)
    sf=0.0
    cf=0.0
    sn1f=0.0
    cn1f=0.0


    for i = 1:nf
        sf = sin(f[i])
        cf = cos(f[i])
        sn1f = sin((N+1)*f[i])
        cn1f = cos((N+1)*f[i])


        μ[i] =  1/2*(sf*cn1f/(cf-1))
               -1/2*sn1f
               -1/2*sf*cf/(cf-1)
               +1/2*sf

        σ[i] = (1/2 + μ[i]^2)*(N+1) - 1/2*cn1f^2
               +1/2*cf*sn1f*cn1f/sf + cn1f*μ[i]
               +μ[i]*sf*sn1f/(cf-1.0)
               - 1/2 - μ[i]^2 -cf*μ[i] -μ[i]*sf^2/(cf-1.0)

    end

    return σ,μ
end


function getSpikeData(IT,f)
    getSpikeData(IT)
end
function getSpikeData(IT)

    N = IT.obs

    μ = Vector{Float64}(N)
    σ = Vector{Float64}(N)
    
    #revisar!!!!
    for i = 1:N
        μ[i] = 1/N
        σ[i] = (N-1.0)*μ[i] + (1.0-μ[i])^2 

    end

    return σ,μ
end

function getSlopeData(IT,f)
    getSlopeData(IT)
end
function getSlopeData(IT)

    N = IT.obs

    μ = Vector{Float64}(N)
    σ = Vector{Float64}(N)
    
    #revisar!!!! -- Passar as expressões do arquivo "APorraToda.pdf".
    for i = 1:N
        μ[i] = 1
        σ[i] = 1

    end

    return σ,μ
end

function getStepData(IT,f)
    getStepData(IT) 
end
function getStepData(IT)

    N = IT.obs

    μ = Vector{Float64}(N-1)
    σ = Vector{Float64}(N-1)
    #L = Vector{Float64}(N-1)#remove
    #U = Vector{Float64}(N-1)

    for i = 1:(N-1)
        μ[i] = (N-i)/N
        σ[i] = sqrt( ( i*μ[i]^2.0+(N-i)*(1.0-μ[i])^2.0 )/N )
        #U[i] = -μ[i]/σ[i]
        #L[i] = (1.0-μ[i])/σ[i]
    end

    return σ,μ#U,L,
end


getData = Vector{Function}(5)

getData[STEP] = getStepData
getData[SLOPE] = getSlopeData
getData[SPIKE] = getSpikeData
getData[SIN] = getSineData
getData[COS] = getCosData
