include("Tools.jl")
include("Constraints.jl")
#printLengths(Connections)
#optimal_path = greedyTSP(Connections[1], Connections)
#println(optimal_path)
SA(Connections[1], Connections)


