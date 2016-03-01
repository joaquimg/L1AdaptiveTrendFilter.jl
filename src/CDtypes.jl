#try with sparse
const STEP = 1
const SPIKE = 2
const SLOPE = 3
const SIN = 4
const COS = 5
const TOTALCOMPONENTS = 5


type iterator
	obs::Int

	components::Vector{Int}
    #??
	elements::Vector{Vector{Int}}#Vector{UnitRange{Int}}#Union{Vector{UnitRange{Int}},Vector{Vector{Int}}}

	#iterating_elements::Union{Vector{UnitRange{Int}},Vector{Vector{Int}}}
	# inicilaiza igual ao de cima
	# iteraness em vez do de cima
	# depois do primeiro lasso modidifica esse

	nelements::Vector{Int}

	isrange::Vector{Bool}


	ttelements::Int

	maxIter::Int
end


type dataCD

	#y0::Float64

	fs::Vector{Float64}
	fc::Vector{Float64}

	μ::Vector{Vector{Float64}}
	σ::Vector{Vector{Float64}}

	function dataCD()
		temp1 = Vector{Float64}[]
		temp2 = Vector{Float64}[]
		for i in 1:TOTALCOMPONENTS
			push!(temp1,zeros(0))
			push!(temp2,zeros(0))
		end
		new(zeros(0),
			zeros(0),
			temp1,
			temp2)
	end
end
