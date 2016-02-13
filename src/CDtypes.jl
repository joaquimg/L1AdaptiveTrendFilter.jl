#try with sparse
const STEP = 1
const SLOPE = 2
const SPIKE = 3
const SIN = 4
const COS = 5


type iterator
	obs::Int
	
	components::Vector{Int}
??
	elements::Union{Vector{UnitRange{Int}},Vector{Vector{Int}}}
	
	nelements::Vector{Int}
	
	isrange::Vector{Bool}


	ttelements::Int

	maxIter::Int
end


type dataCD
	fs::Vector{Float64}
	fc::Vector{Float64}

	μt::Vector{Float64}
	σt::Vector{Float64}

	μl::Vector{Float64}
	σl::Vector{Float64}

	μp::Vector{Float64}
	σp::Vector{Float64}

	μs::Vector{Float64}
	σs::Vector{Float64}

	μc::Vector{Float64}
	σc::Vector{Float64}

	function dataCD()
		new(zeros(0),
			zeros(0),
			zeros(0),zeros(0),
			zeros(0),zeros(0),
			zeros(0),zeros(0),
			zeros(0),zeros(0),
			zeros(0),zeros(0)
			)
	end
end