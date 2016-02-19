include("CDtypes.jl")
include("initializations.jl")
include("CDmainAlg.jl")
include("GetResult.jl")
include("GramMatrix.jl")
include("InnerProducts.jl")
include("Normalization.jl")
include("utils.jl")

srand(10)
y=rand(500)

BCD=CD(y,[1,2,3,4,5], f = [10.0, 11.0])
