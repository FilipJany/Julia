#include("MainGA_Parallel.jl")
include("NewGA.jl")
include("../Other/Constraints.jl")

#fromIteration = 1
#toIteration = 10000
#populationSize = 100
#best = @time Simulate(populationSize, fromIteration, toIteration)
ga = @time GeneticAlgorithm(Connections, 1000, 100)
println(ga)
println(Cost(ga))

#println("Best path: \n\t")
#for i in 1:length(best)
#    print(best[i].id, "->")
#end
#println("with cost: ", Cost(best), "\n")


