include("CDtypes.jl")

#N=10
#IT = iterator(N,[1, 2],[10, 10, 0, 0, 0])

function initIT(N,t = false,l = false,p = false,s = false,c = false,f = [])
	components = Int[]
	nelements = zeros(Int,5)
	obs = N

	if t
		push!(components,1)
		nelements[1] = N-1
	end
	if l
		push!(components,2)
		t ? nelements[2] = N-1 : nelements[2] = N 
	end
	if p
		push!(components,3)
		( t || l ) ? nelements[3] = N-1 : nelements[3] = N
	end
	if s
		push!(components,4)
		nelements[4] = size(f)[1]
	end	
	if c
		push!(components,5)
		nelements[4] = size(f)[1]
	end

	out = iterator(obs,components,nelements,sum(nelements),20)

	return out
end

function initSparse(IT)
	BCD = Array{SparseMatrixCSC{Float64,Int},1}[]
	
	beta_tilde = SparseMatrixCSC{Float64,Int}[]
	beta = SparseMatrixCSC{Float64,Int}[]
	
	activeSet = SparseMatrixCSC{Bool,Int}[]
	
	for i in 1:5
	
		if in(i,IT.components)
	
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
	
	activeSet = Vector{Float64}[]
	
	for i in 1:5
	
		if in(i,IT.components)
	
			push!(beta_tilde,zeros(IT.nelements[i]))
			push!(beta,zeros(IT.nelements[i])
			push!(activeSet,zeros(Bool,IT.nelements[i]))
	
		else
	
			push!(beta_tilde,spzeros(0))
			push!(beta,spzeros(0))
			push!(activeSet,spzeros(Bool,0))
	
		end
	
	end

	return BCD,beta_tilde,beta,activeSet
end

function initXDY(IT,y,data)
	
	xdy = Vector{Float64}[]

	if in(1, IT.components)

		temp = xdy_step(IT,y,data)

		push!(xdy,temp)
	end
	if in(2, IT.components)

		temp = xdy_slope(IT,y,data)

		push!(xdy,temp)
	end
	if in(3, IT.components)

		temp = xdy_spike(IT,y,data)

		push!(xdy,temp)
	end
	if in(4, IT.components)

		temp = xdy_sin(IT,y,data)

		push!(xdy,temp)
	end
	if in(5, IT.components)

		temp = xdy_cos(IT,y,data)

		push!(xdy,temp)
	end

	return xdy
end

function initData()








end





