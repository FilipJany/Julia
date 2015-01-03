include("MainSA_Parallel.jl")

include("../Other/Constraints.jl")

fromIteration = 1
toIteration = 10000
alpha = 0.93
startTemp = 1000.0
minTemp = 0.0001
bestPath = @time Simulate(alpha, startTemp, minTemp, fromIteration, toIteration)

#print("Best path: \n\t", best, "\nwith cost: ", Cost(best), "\n")
#best = Inf
#bestPath = []
#for i in 1:1000
#    sa = SA(CreateCityListFromFile("../Other/Problems/eil51.tsp"), 0.95, 1000.0, 0.00001, 10000)
    #sa = SA(Connections, 0.95, 1000.0, 0.00001, 10000)
#    print(Cost(sa), "\n")
#    if(Cost(sa) < best)
#        best = Cost(sa)
#        bestPath = sa
#    end
#end

print("Best path: \n\t")
for i in 1:length(bestPath)
    print(bestPath[i].id, "->")
end
println("\nwith cost: ", Cost(bestPath))
