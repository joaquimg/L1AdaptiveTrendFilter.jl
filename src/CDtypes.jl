#try with sparse

type iterator
	obs::Int
	
	components::Vector{Int}
	nelements::Vector{Int}

	ttelements::Int

	maxIter::Int
end


type data
	f::Vector{Float64}

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

end