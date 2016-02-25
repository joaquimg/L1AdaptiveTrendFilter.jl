#cd("C:/Users/LabOpto/Documents/SmartGit Projects/src")
using Debug
# using PyPlot
using Gadfly
using DataFrames

include("CDtypes.jl")
include("initializations.jl")
include("CDmainAlg.jl")
include("GetResult.jl")
include("GramMatrix.jl")
include("InnerProducts.jl")
include("Normalization.jl")
include("utils.jl")

# get wind power data
y=readtable("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/wind_power.csv")
y = convert(Array , y[:,1])
y_original = copy(y)

# frequencies
f = 2*pi./collect(6:48)
print(f)

@time BCD, β_best_unbiased, β_best_biased, y_best = l1_adaptive_trend_filter(
  y, [1,2,3,4,5], numλ = 200, f = f, lower_bounds=[-10e+7, -10e+7, -10e+7, -10e+7, -10e+7], upper_bounds=[10e+7, 0, 10e+7, 10e+7, 10e+7]
  )

# components
step = β_best_unbiased[1]
spike = β_best_unbiased[2]
slope = β_best_unbiased[3]
seno = β_best_unbiased[4]
cosseno = β_best_unbiased[5]

# plot (Gadfly)
draw(SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/wind_power_case_study/step.svg", 14inch, 8inch), plot(x=1:length(step), y=step, Geom.point, Geom.line))
draw(SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/wind_power_case_study/spike.svg", 14inch, 8inch), plot(x=1:length(spike),y=spike, Geom.point, Geom.line))
draw(SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/wind_power_case_study/slope.svg", 14inch, 8inch), plot(x=1:length(slope),y=slope, Geom.point, Geom.line))
draw(SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/wind_power_case_study/seno.svg", 14inch, 8inch), plot(x=1:length(seno),y=seno, Geom.point, Geom.line))
draw(SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/wind_power_case_study/cosseno.svg", 14inch, 8inch), plot(x=1:length(cosseno), y=cosseno, Geom.point, Geom.line))
draw(SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/wind_power_case_study/fit.svg", 14inch, 8inch), plot(layer(x=1:length(y_original), y=y_original, Geom.point, Geom.line, Theme(default_color=color("red"))),
                                         layer(x=1:length(y_best), y=y_best, Geom.point, Geom.line, Theme(default_color=color("blue")))))
