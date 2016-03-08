# linear constrained quadratica programming filter EXAMPLE


# using JuMP
# using PyPlot
# using Cbc

# include("LCQP_filter.jl")
using Gadfly

include("../src/CDtypes.jl")
include("../src/initializations.jl")
include("../src/CDmainAlg.jl")
include("../src/GetResult.jl")
include("../src/GramMatrix.jl")
include("../src/InnerProducts.jl")
include("../src/Normalization.jl")
include("../src/utils.jl")

op=1
srand(10)
if op == 0
y=rand(500,1)
n=size(1)
    isempty(n)
elseif op == 1
  n=200
A=triu(ones(n,n))
lines=sort!(collect(Set(round((n-1)*rand(20,1)+1, 0))))
    y=A[:,lines]*(40*rand(size(lines)[1],1)+1)+10*(randn(n,1))
draw(
    SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/teste.svg", 7inch, 4.5inch),
    plot(y=y,Geom.line)
     )
end

y=y[:,1]

@time BCD, y_path, β_best, y_best, λ_best, γ_best = l1_adaptive_trend_filter(
  y, [1, 3], numλ=100, numγ=4, verbose=true, lower_bounds=[-Inf, -Inf, -Inf], upper_bounds=[Inf, Inf, 0]
  )

draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/teste2.svg", 7inch, 4.5inch),
  plot(
    layer(x=1:length(y_best), y=y_best, Geom.line, Theme(default_color=color("blue"))),
    layer(x=1:length(y), y=y, Geom.line, Theme(default_color=color("red"))),
    Guide.xlabel("Hours"), Guide.ylabel("MW"), Guide.title("Wind Power Fit"),
    Coord.Cartesian(xmin=0,xmax=340)
    ))

# components
step = β_best[1]
# spike = β_best[2]
slope = β_best[3]

draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/step.svg", 14inch, 8inch),
  plot(
    x=1:length(step), y=step, Geom.point, Geom.line,
    Guide.xlabel("Hours"), Guide.ylabel("Step size"),
    Coord.Cartesian(xmin=0,xmax=length(y))
    ))
# draw(
#   SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/spike.svg", 14inch, 8inch),
#   plot(x=1:length(spike),y=spike, Geom.point, Geom.line, Guide.xlabel("Hours")))
draw(
  SVG("/Users/mariosouto/Dropbox/SAM/L1_Adaptive_Trend_Filter/WindPower_CaseStudy/slope.svg", 14inch, 8inch),
  plot(x=1:length(slope),y=slope, Geom.point, Geom.line)
  )
