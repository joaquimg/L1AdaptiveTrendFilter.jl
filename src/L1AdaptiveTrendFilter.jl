
__precompile__()

module L1AdaptiveTrendFilter

  ### exports
  export

  l1_adaptive_trend_filter

  ### include source files
  
  include("CDtypes.jl")
  include("initializations.jl")
  include("CDmainAlg.jl")
  include("GetResult.jl")
  include("GramMatrix.jl")
  include("InnerProducts.jl")
  include("Normalization.jl")
  include("utils.jl")

end