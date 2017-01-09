

function initIT_range(N::Int, components::Vector{Int}, f::Vector{Float64} = Float64[]; MAXITER::Int=100)
	#components = Int[]
	nelements = zeros(Int, TOTALCOMPONENTS)
	elements = Vector{Vector{Int}}(TOTALCOMPONENTS)
	obs = N

	if STEP in components
		#push!(components,STEP)
		nelements[STEP] = N-1
		elements[STEP] = collect( 1:(N-1) )
	else
		nelements[STEP] = 0
		elements[STEP] = Int[]
	end

	if SPIKE in components
		#push!(components,SPIKE)
		if ( STEP in components || SLOPE in components ) # last step and last spike are the same
			nelements[SPIKE] = N-1
			elements[SPIKE] = collect( 1:(N-1) )
		else
			nelements[SPIKE] = N
			elements[SPIKE] = collect( 1:N )
		end
	else
		nelements[SPIKE] = 0
		elements[SPIKE] = Int[]
	end

	if SLOPE in components
		#push!(components,SLOPE)
		nelements[SLOPE] = N-1
		elements[SLOPE] = collect( 1:(N-1) )
	else
		nelements[SLOPE] = 0
		elements[SLOPE] = Int[]
	end

	if SIN in components
		#push!(components,SIN)
		nelements[SIN] = length(f)
		elements[SIN] = collect( 1:length(f))
	else
		nelements[SIN] = 0
		elements[SIN] = Int[]
	end

	if COS in components
		#push!(components,COS)
		nelements[COS] = length(f)
		elements[COS] = collect( 1:length(f))
	else
		nelements[COS] =  0
		elements[COS] = Int[]
	end

	return iterator(obs, components, elements, nelements, ones(Bool,TOTALCOMPONENTS), sum(nelements),MAXITER)
end

function initSparse(IT::iterator)

	BCD = Array{SparseMatrixCSC{Float64,Int},1}[]

	beta_tilde = SparseMatrixCSC{Float64,Int}[]
	beta = SparseMatrixCSC{Float64,Int}[]

	activeSet = SparseMatrixCSC{Bool,Int}[]

	for i in 1:TOTALCOMPONENTS

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

	return BCD, beta_tilde, beta, activeSet
end

function initDense(IT)
	BCD = Array{Vector{Float64},1}[]

	beta_tilde = Vector{Float64}[]
	beta = Vector{Float64}[]

	activeSet = Vector{Bool}[]

	for i in 1:TOTALCOMPONENTS

		if i in IT.components

			push!(beta_tilde, zeros(IT.nelements[i]))
			push!(beta, zeros(IT.nelements[i]))
			push!(activeSet,zeros(Bool,IT.nelements[i]))

		else

			push!(beta_tilde,zeros(0))
			push!(beta, zeros(0))
			push!(activeSet,zeros(Bool,0))

		end

	end

	return BCD, beta_tilde, beta, activeSet
end

function initXDY(IT,y,data)

	xdy0 = Vector{Float64}[]

	for i in 1:TOTALCOMPONENTS
		if in(i,IT.components)

			push!(xdy0, xdy[i](IT,y,data) )

		else
			push!(xdy0, Float64[] )
		end
		#temp)
	end


	return xdy0
end

#change for no preallocation and simply edit field
#requires changing normalizations.jl (returns)
function initData(IT::iterator, fs::Vector{Float64}=Float64[], fc::Vector{Float64}=Float64[])

	d = dataCD()
	for i in IT.components
		d.σ[i],d.μ[i] = getData[i](IT,fs)
	end
	d.fs = fs
	d.fc = fc
	
	return d
end




