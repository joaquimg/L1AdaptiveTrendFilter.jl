
__precompile__()

module L1AdaptiveTrendFilter

    ### exports
    export
    # main function
    l1_adaptive_trend_filter,
    # input names
    STEP, SPIKE, SLOPE, SIN, COS, TOTALCOMPONENTS

    ### include source files
    
    # types definitions
    include("CDtypes.jl")

    # initializes basic types and vectors (betas, paths...)
    include("initializations.jl")

    # core algorithm
    include("CDmainAlg.jl")

    # computes estimated y from betas
    include("GetResult.jl")

    # functions for all terms in the gram matrix (cross corelations)
    include("GramMatrix.jl")

    # computes innner products between candidates
    include("InnerProducts.jl")

    # computes means and varianeces
    include("Normalization.jl")


    include("utils.jl")

end