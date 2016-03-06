#cd("C:/Users/LabOpto/Documents/SmartGit Projects/src")
# using Debug
# using PyPlot
using Gadfly
using DataFrames

include("../src/CDtypes.jl")
include("../src/initializations.jl")
include("../src/CDmainAlg.jl")
include("../src/GetResult.jl")
include("../src/GramMatrix.jl")
include("../src/InnerProducts.jl")
include("../src/Normalization.jl")
include("../src/utils.jl")

# get wind power data
y=readtable("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/wind_power.csv", header=false)
y = convert(Array , y[:,1])
y_original = copy(y)

# frequencies
f = 2*pi./collect(6:48)

@time BCD, y_path, β_best, y_best, λ_best, γ_best = l1_adaptive_trend_filter(
  y, [1, 2, 3, 4, 5], numλ=100, numγ=4, f=f, verbose=true
  )

# components
step = β_best[1]
spike = β_best[2]
slope = β_best[3]
seno = β_best[4]
cosseno = β_best[5]

# plot (Gadfly)
draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/step.svg", 14inch, 8inch),
  plot(
    x=1:length(step), y=step, Geom.point, Geom.line,
    Guide.xlabel("Hours"), Guide.ylabel("Step size"),
    Coord.Cartesian(xmin=0,xmax=length(y_original))
    ))
draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/spike.svg", 14inch, 8inch),
  plot(x=1:length(spike),y=spike, Geom.point, Geom.line, Guide.xlabel("Hours")))
draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/slope.svg", 14inch, 8inch),
  plot(x=1:length(slope),y=slope, Geom.point, Geom.line)
  )
draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/trigonometric.svg", 7inch, 4.5inch),
  plot(
    layer(x=1:length(seno),y=seno, Geom.line, Theme(default_color=color("red"))),
    layer(x=1:length(cosseno), y=cosseno, Geom.line, Theme(default_color=color("blue"))),
    Guide.xlabel("frequency (ω)"), Guide.ylabel("Amplitude"), Guide.title("Trigonometric components"),
    Coord.Cartesian(xmin=0,xmax=length(f))
    ))
draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/fit.svg", 7inch, 4.5inch),
  plot(
    layer(x=1:length(y_original), y=y_original, Geom.line, Theme(default_color=color("red"))),
    layer(x=1:length(y_best), y=y_best, Geom.line, Theme(default_color=color("blue"))),
    Guide.xlabel("Hours"), Guide.ylabel("MW"), Guide.title("Wind Power Fit"),
    Coord.Cartesian(xmin=0,xmax=340)
    ))

# write results
writecsv("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/y_best.csv", y_best)
writecsv("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/cos.csv", cosseno)
writecsv("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/seno.csv", seno)
writecsv("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/step.csv", step)
