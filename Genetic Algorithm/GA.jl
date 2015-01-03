include("../Other/Constraints.jl")

export GeneticAlgorithm

function GeneticAlgorithm(availableCities::Array{City, 1}, maxSteps::Int64, populationSize::Int64)
    i::Int64 = 0
    bestDistance::Float64 = Inf
    bestIteration::Float64 = Inf
    population::Array{Array{City, 1}, 1} = GenerateRandomPopulation(availableCities, populationSize)
    bestTour::Array{City, 1} = []
    while i < maxSteps
        population = Evolve(population)
        if(IsAnyBetterThanBest(bestDistance, population))
            bestDistance = Cost(GetFittestFromPopulation(population))
            bestIteration = i
            bestTour = GetFittestFromPopulation(population)
        end
        i += 1
    end
    return bestTour
end

function IsAnyBetterThanBest(best, recentPopulation)
    recentBest = Cost(GetFittestFromPopulation(recentPopulation))
    if(recentBest < best)
        return true;
    end
    return false;
end

function Mutate(child::Array{City, 1})
    mutationPoint1::Int64 = abs(rand(Int64) % length(child)) + 1
    mutationPoint2::Int64 = abs(rand(Int64) % length(child)) + 1
    temp = child[mutationPoint1]
    child[mutationPoint1] = child[mutationPoint2]
    child[mutationPoint2] = temp
    return child
end

function Reproduce(x, y)
    child::Array{City, 1} = []
    startPos::Int64 = abs(rand(Int64) % length(x)) + 1
    endPos::Int64 = abs(rand(Int64) % length(x)) + 1
    for i = 1 : length(x)
        if(startPos < endPos && i > startPos && i < endPos)
            push!(child, x[i])
        elseif (startPos > endPos)
            if(!(i < startPos && i > endPos))
                push!(child, x[i])
            end
        end
    end
    for j = 1 : length(y)
        if(!isCityInPath(y[j], child))
            push!(child, y[j])
        end
    end
    return child
end

function Cost(tour::Array{City, 1})
    summ::Float64 = 0.0
    for i = 1 : length(tour)-1
        summ += GetLength(tour[i], tour[i+1])
    end
    summ += GetLength(tour[length(tour)], tour[1])
    return summ
end

function FitnessFunction(tour::Array{City, 1})
    return 1 / Cost(tour)
end

function GetFittestFromPopulation(population)
    fittest::Array{City, 1} = population[1]
    for i = 2 : length(population)
        if((1/FitnessFunction(fittest)) > (1/FitnessFunction(population[i])))
            fittest = population[i]
        end
    end
    return fittest
end

function Selection(population)
    k::Int64 = abs(rand(Int64) % (length(population) - 1)) + 1
    tournament::Array{Array{City, 1}, 1} = []
    for i = 1 : k
        id = abs(rand(Int64) % (length(population) - 1)) + 1
        push!(tournament, population[id])
    end
    return GetFittestFromPopulation(tournament)
end

function Evolve(population::Array{Array{City, 1}, 1})
    mutationProbability::Float64 = 0.5
    newPopulation::Array{Array{City, 1}, 1} = []
    for i = 1 : length(population)
        parent1::Array{City, 1} = Selection(population)
        parent2::Array{City, 1} = Selection(population)
        child::Array{City, 1} = Reproduce(parent1, parent2)#todo
        if(rand(Float64) % 1.0 < mutationProbability)
            child = Mutate(child)#todo
        end
        push!(newPopulation, child)
    end
    return newPopulation
end

function GenerateRandomPopulation(availableCities::Array{City, 1}, populationSize::Int64)
    population::Array{Array{City, 1}, 1} = []
    for i = 1 : populationSize
        push!(population, GetRandomTour(availableCities, length(availableCities)))
    end
    return population
end

function GetRandomTour(availableCities, tourLength)
    tour = []
    i = 0
    while i < length(availableCities)
        id = abs(rand(Int64) % tourLength ) + 1
        if(!isCityInPath(availableCities[id], tour))
            push!(tour, availableCities[id])
            i += 1
        end
    end
    return tour
end

function isCityInPath(city::City, path)
    for i = 1 : length(path)
        if(path[i] == city)
            return true;
        end
    end
    return false;
end

function GetLength(p::City, q::City)
    dx = p.x - q.x
    dy = p.y - q.y
    return floor(sqrt(dx*dx + dy*dy))
end
