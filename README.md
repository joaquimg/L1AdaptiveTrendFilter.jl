# L1AdaptiveTrendFilter.jl

L1AdaptiveTrendFilter.jl is a Julia package for filtering and identifying unknown underlying trends of a given noisy signal. It is originally a supplemental material for the paper "l1 Adaptive Trend Filter via Fast Coordinate Descent".


In order to use it, download the package in julia via `Pkg.clone("git://github.com/joaquimg/L1AdaptiveTrendFilter.jl.git")`. Then simply `using L1AdaptiveTrendFilter`

So far there is only one function: `l1_adaptive_trend_filter` 

The inputs are described below:

