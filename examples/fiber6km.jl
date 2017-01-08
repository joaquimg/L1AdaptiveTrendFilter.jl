push!(LOAD_PATH, "/Users/joaquimdiasgarcia/Repositories")
cd("/Users/joaquimdiasgarcia/Repositories/L1AdaptiveTrendFilter.jl/examples")

using L1AdaptiveTrendFilter

f=open("fiber6km_data.txt")
data=readlines(f)
close(f)
data_float = map(x->parse(Float64, x),data)

println("Begin")

@time BCD, y_path, b_best, y_best, b_best, b_best = l1_adaptive_trend_filter(data_float, [STEP, SLOPE], numλ=100, numγ=3, verbose=true)
