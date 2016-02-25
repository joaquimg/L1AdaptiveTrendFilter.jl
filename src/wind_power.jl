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

# frequencies
f = 2*pi./collect(2:100)
print(f)

@time BCD,β1,β2, y_best = CD(y,[1,2,3,4,5], numλ = 100, f = f)

# components
step = β1[1]
spike = β1[2]
slope = β1[3]
seno = β1[4]
cosseno = β1[5]

# plot (Gadfly)
draw(SVG("step.svg", 14inch, 8inch), plot(x=1:length(step), y=step, Geom.point, Geom.line))
draw(SVG("spike.svg", 14inch, 8inch), plot(x=1:length(spike),y=spike, Geom.point, Geom.line))
draw(SVG("slope.svg", 14inch, 8inch), plot(x=1:length(slope),y=slope, Geom.point, Geom.line))
draw(SVG("seno.svg", 14inch, 8inch), plot(x=1:length(seno),y=seno, Geom.point, Geom.line))
draw(SVG("cosseno.svg", 14inch, 8inch), plot(x=1:length(cosseno), y=cosseno, Geom.point, Geom.line))
draw(SVG("y_hat.svg", 14inch, 8inch), plot(x=1:length(y_best), y=y_best, Geom.point, Geom.line))
draw(SVG("y.svg", 14inch, 8inch), plot(x=1:length(y), y=y, Geom.point, Geom.line))
