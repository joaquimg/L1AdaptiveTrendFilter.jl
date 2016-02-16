include("CDtypes.jl")
include("CDmainAlg.jl")
include("GetResult.jl")
include("GramMatrix.jl")
include("initializations.jl")
include("InnerProducts.jl")
include("Normalization.jl")
include("utils.jl")

op=1
srand(10)
y=rand(500,1)
BCD=CD(y)
