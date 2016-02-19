using Debug


include("CDtypes.jl")
include("initializations.jl")
include("CDmainAlg.jl")
include("GetResult.jl")
include("GramMatrix.jl")
include("InnerProducts.jl")
include("Normalization.jl")
include("utils.jl")

srand(10)
y=rand(10)



BCD=CD(y,[1,2,3], f = [10.0, 11.0])
