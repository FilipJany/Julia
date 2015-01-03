include("MainIterative_Parallel.jl")

include("../Other/Tools.jl")
include("../Other/Constraints.jl")

fromIteration = 1
toIteration = 100

best = @time Simulate(fromIteration, toIteration)

print("Best path: \n\t", best, "\nwith cost: ", Cost(best), "\n")
