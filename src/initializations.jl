include("CDtypes.jl")

#N=10
#IT = iterator(N,[1, 2],[10, 10, 0, 0, 0])

function initIT_range(N,components, f = Vector{Float64}(0) )
	#components = Int[]
	nelements = zeros(Int,5)
	elements = Vector{Any}(5)
	obs = N

	if STEP in components
		#push!(components,STEP)
		nelements[STEP] = N-1
		elements[STEP] = 1:(N-1)
	end
	if SLOPE in components
		#push!(components,SLOPE)
		if STEP in components
			nelements[SLOPE] = N-1
			elements[SLOPE] = 1:(N-1)
		else
			nelements[SLOPE] = N
			elements[SLOPE] = 1:N
		end
	end
	if SPIKE in components
		#push!(components,SPIKE)
		if ( STEP in components || SLOPE in components )
			nelements[SPIKE] = N-1
			elements[SPIKE] = 1:(N-1)
		else
			nelements[SPIKE] = N
			elements[SPIKE] = 1:N
		end
	end
	if SIN in components
		#push!(components,SIN)
		nelements[SIN] = size(f)[1]
		elements[SIN] = 1:size(f)[1]
	end
	if COS in components
		#push!(components,COS)
		nelements[COS] = size(f)[1]
		elements[COS] = 1:size(f)[1]
	end

	out = iterator(obs,components,elements,nelements,ones(Bool,TOTALCOMPONENTS),sum(nelements),20)

	return out
end

function initSparse(IT)
	BCD = Array{SparseMatrixCSC{Float64,Int},1}[]

	beta_tilde = SparseMatrixCSC{Float64,Int}[]
	beta = SparseMatrixCSC{Float64,Int}[]

	activeSet = SparseMatrixCSC{Bool,Int}[]

	for i in 1:5

		if i in IT.components

			push!(beta_tilde,spzeros(IT.nelements[i],1))
			push!(beta,spzeros(IT.nelements[i],1))
			push!(activeSet,spzeros(Bool,IT.nelements[i],1))

		else

			push!(beta_tilde,spzeros(0,0))
			push!(beta,spzeros(0,0))
			push!(activeSet,spzeros(Bool,0,0))

		end

	end

	return BCD,beta_tilde,beta,activeSet
end

function initDense(IT)
	BCD = Array{Vector{Float64},1}[]

	beta_tilde = Vector{Float64}[]
	beta = Vector{Float64}[]

	activeSet = Vector{Bool}[]

	for i in 1:5

		if i in IT.components

			push!(beta_tilde,zeros(IT.nelements[i]))
			push!(beta,zeros(IT.nelements[i]))
			push!(activeSet,zeros(Bool,IT.nelements[i]))

		else

			push!(beta_tilde,zeros(0))
			push!(beta,zeros(0))
			push!(activeSet,zeros(Bool,0))

		end

	end

	return BCD,beta_tilde,beta,activeSet
end

function initXDY(IT,y,data)

	xdy0 = Vector{Float64}[]

	for i in 1:TOTALCOMPONENTS
		if i in  IT.components

			temp = xdy[i](IT,y,data)
		else
			temp = Vector{Float64}(0)
		end
		push!(xdy0,temp)
	end


	return xdy0
end

#change for no preallocation and simply edit field
#requires changing normalizations.jl (returns)
function initData(IT,fs=[],fc=[])

	d = dataCD()
	for i in IT.components
		d.σ[i],d.μ[i] = getData[i](IT,fs)
	end
	d.fs = fs
	d.fc = fc
	return d
end




