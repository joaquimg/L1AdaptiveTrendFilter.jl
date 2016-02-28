# linear constrained quadratica programming filter EXAMPLE


using JuMP
using PyPlot
using Cbc

include("LCQP_filter.jl")


op=1
srand(10)
if op == 0
y=rand(500,1)
n=size(1)
    isempty(n)
elseif op == 1 
A=triu(ones(5000,5000))
lines=sort!(collect(Set(int(4999*rand(20,1)+1))))
    y=A[:,lines]*(40*rand(size(lines)[1],1)+1)+(randn(5000,1))
plot(y)
end


x,w,u,a=paramfree_l1(y)


fig = figure()
p1=plot(y,color="b");
p2=plot(w,color="r");



