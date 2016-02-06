include("CDtypes.jl")

N=10
IT = iterator(N,[1, 2],[10, 10, 0, 0, 0])


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