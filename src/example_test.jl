cd("C:/Users/LabOpto/Documents/SmartGit Projects/src")

using PyPlot

include("CDtypes.jl")
include("initializations.jl")
include("CDmainAlg.jl")
include("GetResult.jl")
include("GramMatrix.jl")
include("InnerProducts.jl")
include("Normalization.jl")
include("utils.jl")

y=readcsv("C:/Users/LabOpto/Documents/SmartGit Projects/TestBench.csv")
y = y[:,1]

run = 6

if run == 1
	##TESTE STEP + SPIKE
	@time BCD,β1,β2 = CD(y,[1,2], numλ = 100)
	a = β1[1] + β1[2]
	b = β2[1] + β2[2]
	plot(a)
	plot(b)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia1.csv",a)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia2.csv",b)
end

if run == 2
	##TESTE SLOPE
	@time BCD,β1,β2 = CD(y,[3], numλ = 100)
	a = β1[3]
	b = β2[3]
	plot(a)
	plot(b)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia1.csv",a)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia2.csv",b)
end

if run == 3
	##TESTE STEP + SPIKE + SLOPE
	@time BCD,β1,β2 = CD(y,[1,2,3], numλ = 100)
	a = β1[1] + β1[2] + β1[3]
	b = β2[1] + β2[2] + β2[3]
	plot(a)
	plot(b)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia1.csv",a)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia2.csv",b)
end

if run == 4
	##TESTE STEP + SLOPE
	@time BCD,β1,β2 = CD(y,[1,3], numλ = 100)
	a = β1[1] + β1[3]
	b = β2[1] + β2[3]
	plot(a)
	plot(b)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia1.csv",a)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia2.csv",b)
end

if run == 5
	##TESTE SPIKE + SLOPE
	@time BCD,β1,β2 = CD(y,[2,3], numλ = 100)
	a = β1[2] + β1[3]
	b = β2[2] + β2[3]
	plot(a)
	plot(b)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia1.csv",a)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia2.csv",b)
end

if run == 6
	## WARNING!!!!
	##TESTE SINE
	f = 2*pi./collect(1:100)
	@time BCD,β1,β2 = CD(y,[4], numλ = 100, f = f)
	a = β1[4]
	b = β2[4]
	plot(a)
	plot(b)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia1.csv",a)
	writecsv("C:/Users/LabOpto/Documents/SmartGit Projects/ResultsJulia2.csv",b)
end