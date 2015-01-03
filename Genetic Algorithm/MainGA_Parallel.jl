include("GA.jl")
include("../Other/Constraints.jl")
include("../Other/Tools.jl")

export master, Simulate

function Simulate(populationSize::Int64, fromIteration::Int64, toIteration::Int64)

    cities = CreateCityListFromFile("../Other/Problems/eil51.tsp")

    results::Array{Array{City, 1}, 1} = []

    function GetData()
        if(length(results) == 0)
            return cities
        else
            return results[length(results)]
        end
    end

    midPoints = [0, 30, 50, 80]

    for i in 1:length(midPoints) - 1
        push!(results, master(GetData(), populationSize, midPoints[i], midPoints[i + 1]))
    end
    best = GetBest(results)
    return best
end

function master(startingPath::Array{City, 1}, populationSize::Int64, fromIteration::Int64, toIteration::Int64)
    results::Array{Array{City, 1}, 1} = []
    @sync for p in 1:nprocs()
        if p != myid() || nprocs() == 1
            @async for i in 1:1
                push!(results, worker(startingPath, populationSize, fromIteration, toIteration))
            end
        end
    end
    return GetBest(results)
end

function worker(startingPath::Array{City, 1}, populationSize::Int64, fromIteration::Int64, toIteration::Int64)
    return GeneticAlgorithm(startingPath, (toIteration-fromIteration), populationSize)
end

function GetBest(results::Array{Array{City, 1}, 1})
    bestTour::Array{City, 1} = results[1]
    currentBest::Float64 = Inf
    for i in 1:length(results)
        currentBest = Cost(results[i])
        if(currentBest <= Cost(bestTour))
            bestTour = results[i]
        end
    end
    return bestTour
end
