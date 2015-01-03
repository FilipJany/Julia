include("TO.jl")
include("../Other/Tools.jl")
include("../Other/Constraints.jl")



@everywhere function summElems()
    if(nprocs() > 1)
        cities = CreateCityListFromFile("../Other/Problems/eil51.tsp")
        results::Array{City, 1} = []
        workArray::Array{Int64, 1} = ManageWork(cities, nprocs()-1)
        println(workArray)
        results = TwoOpt(cities, 10, workArray)
    else
        println("Please add more processors for work!")
    end
    return results
end

