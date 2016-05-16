
# using Gadfly
using DataFrames
# using L1AdaptiveTrendFilter

include("../src/CDtypes.jl")
include("../src/initializations.jl")
include("../src/CDmainAlg.jl")
include("../src/GetResult.jl")
include("../src/GramMatrix.jl")
include("../src/InnerProducts.jl")
include("../src/Normalization.jl")
include("../src/utils.jl")

y = readtable(
  "/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/OpticFiber_CaseStudy/ExperimentalData.csv",
  header=false
  )
y = convert(Array, y[:, 2])
y_original = copy(y)

@time BCD, y_path, β_best, y_best, λ_best, γ_best = l1_adaptive_trend_filter(
  y, [1, 3], numλ=100, numγ=3, verbose=true, upper_bounds=[0, 0,0,0,0]
  )

# components
step = β_best[1]
# spike = β_best[2]
slope = β_best[3]

# # plot (Gadfly)
# draw(
#   SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/OpticFiber_CaseStudy/step.svg", 14inch, 8inch),
#   plot(
#     x=1:length(step), y=step, Geom.point, Geom.line,
#     Guide.xlabel("Hours"), Guide.ylabel("Step size"),
#     Coord.Cartesian(xmin=0,xmax=length(y_original))
#     ))
# # draw(
# #   SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/OpticFiber_CaseStudy/spike.svg", 14inch, 8inch),
# #   plot(x=1:length(spike),y=spike, Geom.point, Geom.line, Guide.xlabel("Hours")))
# draw(
#   SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/OpticFiber_CaseStudy/slope.svg", 14inch, 8inch),
#   plot(x=1:length(slope),y=slope, Geom.point, Geom.line)
#   )
# draw(
#   SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/OpticFiber_CaseStudy/fit.svg", 7inch, 4.5inch),
#   plot(
#     layer(x=1:length(y_original), y=y_original, Geom.line, Theme(default_color=color("red"))),
#     layer(x=1:length(y_best), y=y_best, Geom.line, Theme(default_color=color("blue"))),
#     Guide.xlabel("Hours"), Guide.ylabel("MW"), Guide.title("Wind Power Fit"),
#     Coord.Cartesian(xmin=0,xmax=8000)
#     ))

# write results
writecsv("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/OpticFiber_CaseStudy/y_best.csv", y_best)
writecsv("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/OpticFiber_CaseStudy/step.csv", step)
