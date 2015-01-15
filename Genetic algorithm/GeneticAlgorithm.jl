#File: SimulatedAnnealing.jl
#Author: Filip Jany 194208 Wrocław University of Technology
#Purpose: Engeneering Thesis (2015)
#Description: Implementation of Generic Algorithm and it's helping functions
include("../Other/Constraints.jl")

export GeneticAlgorithm

#Parallel implementation of iterative Genetic Algorithm.
#Input:
#       availableCities: an array of cities to visit
#       maxSteps : number of iteartions that algoritm should perform
#       populationSize : size of population
#       mutProb : probability of mutation
#Output:
#   Best(with smallest cost) of calculated tours
@everywhere function GeneticAlgorithm(availableCities::Array{City, 1}, maxSteps::Int64, populationSize::Int64, mutProb::Float64)
    i::Int64 = 0
    results = []
    bestDistance::Float64 = Inf
    bestIteration::Float64 = Inf
    population::Array{Array{City, 1}, 1} = GenerateRandomPopulation(availableCities, populationSize)
    bestTour::Array{City, 1} = []
    workArray::Array{Int64, 1} = ManageWork(maxSteps)
    @sync begin
        for p in 1:nprocs()
            if p != myid() || nprocs() == 1
                @async begin
                    if(length(workArray) <= 2)
                        push!(results, remotecall_fetch(p, WorkerTask, workArray[1], workArray[2], population, mutProb))
                    else
                        push!(results, remotecall_fetch(p, WorkerTask,  workArray[p-1], workArray[p], population, mutProb))
                    end
                end
            end
        end
    end
    return GetBest(results)
end

#Serial implementation of iterative Genetic Algorithm.
#Input:
#       availableCities: an array of cities to visit
#       maxSteps : number of iteartions that algoritm should perform
#       populationSize : size of population
#       mutProb : probability of mutation
#Output:
#   bestTout : best(with smallest cost) of calculated tours
function GeneticAlgorithmSerial(availableCities::Array{City, 1}, maxSteps::Int64, populationSize::Int64, mutProb::Float64)
    population::Array{Array{City, 1}, 1} = GenerateRandomPopulation(availableCities, populationSize)
    bestDistance::Float64 = Inf
    bestTour::Array{City, 1} = []
    i::Int64 = 0
    for i in 1:maxSteps
        population = Evolve(population, length(population), mutProb)
        if(IsAnyBetterThanBest(bestDistance, population))
            bestDistance = Cost(GetFittestFromPopulation(population))
            bestTour = GetFittestFromPopulation(population)
        end
    end
    return bestTour
end

#Finds and returns best (with smaller cost) tour from set of results.
#Input:
#   results: an array of results
#Output:
#   best: best (smallest cost) solution
@everywhere function GetBest(results)
    best::Array{City, 1} = results[1]
    for i in 1:length(results)
        if(Cost(results[i]) < Cost(best))
            best = results[i]
        end
    end
    return best
end

@everywhere function ManageWork(maxSteps::Int64)
    workArray::Array{Int64, 1} = []
    np::Int64 = nprocs()
    stepsPerProcess = floor(maxSteps/np)
    current::Int64 = 1
    for i in 1:np-1
        current = stepsPerProcess*i
        push!(workArray, current)
    end
    if(current < maxSteps)
        push!(workArray, maxSteps)
    end
    return workArray
end

#Function specifies a worker task.
#Input:
#       fromIter: number of iteration to start with
#       toIter: number of iteration to end with
#       population : array of tours(chromosomes)
#       mutProb : mutation probablity
#Output:
#   bestTour: best found tour
@everywhere function WorkerTask(fromIter::Int64, toIter::Int64, population::Array{Array{City, 1}, 1}, mutProb::Float64)
    bestDistance::Float64 = Inf
    bestTour::Array{City, 1} = []
    i::Int64 = 0
    for i in fromIter:toIter
        population = Evolve(population, length(population), mutProb)
        if(IsAnyBetterThanBest(bestDistance, population))
            bestDistance = Cost(GetFittestFromPopulation(population))
            bestTour = GetFittestFromPopulation(population)
        end
    end
    return bestTour
end

#Function determines if there is a better soution than given in population
#Input:
#       best: tour to be compared with
#       recentPopulation : population where function will search for better tour
#Output:
#   true if better found, false otherwise
@everywhere function IsAnyBetterThanBest(best, recentPopulation)
    recentBest = Cost(GetFittestFromPopulation(recentPopulation))
    if(recentBest < best)
        return true;
    end
    return false;
end

#Function mutates given tour by swapping two randomly choosen cities
#Input:
#       child : tour (chromosome) to mutate
#Output:
#   child : mutated input tour
@everywhere function Mutate(child::Array{City, 1})
    mutationPoint1::Int64 = abs(rand(Int64) % length(child)) + 1
    mutationPoint2::Int64 = abs(rand(Int64) % length(child)) + 1
    temp = child[mutationPoint1]
    child[mutationPoint1] = child[mutationPoint2]
    child[mutationPoint2] = temp
    return child
end

#Implementation of OX operator to reproduce parents
#Input:
#       parent1 : first parent that genes will come from
#       parent2 : second parent that genes will come from
#Output:
#   A set of two children generated from parent's genes
@everywhere function Reproduce(parent1::Array{City, 1}, parent2::Array{City, 1})
    cutPoint1::Int64 = rand(1:length(parent1))
    cutPoint2::Int64 = rand(1:length(parent1))

    child1 = InitChild(length(parent1))
    child2 = InitChild(length(parent1))

    child1Base::Array{City, 1} = []
    child2Base::Array{City, 1} = []

    if (cutPoint1 > cutPoint2)
        temp = cutPoint1
        cutPoint1 = cutPoint2
        cutPoint2 = temp
    end

    for i in cutPoint1:cutPoint2
        child1[i] = parent2[i]
        child2[i] = parent1[i]
    end

    citiesToInsert1::Array{City, 1} = []
    citiesToInsert2::Array{City, 1} = []
    #Get rest of genes
    for i in cutPoint2+1:length(parent2)
        if(!IsAlreadyInChild(parent2[i].id, child2))
            push!(citiesToInsert1, parent2[i])
        end
    end
    for i in 1:cutPoint2
        if(!IsAlreadyInChild(parent2[i].id, child2))
            push!(citiesToInsert1, parent2[i])
        end
    end

    for i in cutPoint2+1:length(parent1)
        if(!IsAlreadyInChild(parent1[i].id, child1))
            push!(citiesToInsert2, parent1[i])
        end
    end
    for i in 1:cutPoint2
        if(!IsAlreadyInChild(parent1[i].id, child1))
            push!(citiesToInsert2, parent1[i])
        end
    end

    #Put remaining genes in children
    j::Int64 = 1
    for i in cutPoint2+1:length(parent1)
       child1[i] = citiesToInsert2[j]
       j += 1
    end
    for i in 1:cutPoint1-1
        child1[i] = citiesToInsert2[j]
        j += 1
    end
    j = 1
    for i in cutPoint2+1:length(parent1)
       child2[i] = citiesToInsert1[j]
       j += 1
    end
    for i in 1:cutPoint1-1
        child2[i] = citiesToInsert1[j]
        j += 1
    end
    return (child1, child2)
end

#Function calculates cost of given tour.
#Input:
#   tour: an array of cities
#Output:
#   summ: cost of a tour
@everywhere function Cost(tour::Array{City, 1})
    summ::Float64 = 0.0
    for i = 1 : length(tour)-1
        summ += GetLength(tour[i], tour[i+1])
    end
    summ += GetLength(tour[length(tour)], tour[1])
    return summ
end

#Function determines fittness of given tour
#Input:
#   tour: an array of cities
#Output:
#   cost of a tour (its fittness)
@everywhere function FitnessFunction(tour::Array{City, 1})
    return Cost(tour)
end

#Function finds and returns best tour from population
#Input:
#   population : a set of tours (chromosomes)
#Output:
#   fittest : best solution that has been found
@everywhere function GetFittestFromPopulation(population)
    fittest::Array{City, 1} = population[1]
    for i = 2 : length(population)
        if(FitnessFunction(fittest) > FitnessFunction(population[i]))
            fittest = population[i]
        end
    end
    return fittest
end

#Function selects two parents to reproduce by using tournament selection
#Input:
#   population : a set of tours possible to choose
#Output:
#   Fittest tour from tournament
@everywhere function Selection(population)
    k::Int64 = abs(rand(Int64) % (length(population) - 1)) + 1
    tournament::Array{Array{City, 1}, 1} = []
    for i = 1 : k
        id = abs(rand(Int64) % (length(population) - 1)) + 1
        push!(tournament, population[id])
    end
    return GetFittestFromPopulation(tournament)
end

#Function is responsible for evolving population into new one
#Input:
#   population : set of tours
#   populationSize : size of population
#   mutationProbability : probability of mutation
#Output:
#   newPopulation : newly ctreated population
@everywhere function Evolve(population::Array{Array{City, 1}, 1}, populationSize::Int64, mutationProbability::Float64)
    newPopulation::Array{Array{City, 1}, 1} = []
    for i = 1 : floor(length(population)/2) + 1
        parent1::Array{City, 1} = Selection(population)
        parent2::Array{City, 1} = Selection(population)
        r = Reproduce(parent1, parent2)
        child1::Array{City, 1} = r[1]
        child2::Array{City, 1} = r[2]
        if(rand(Float64) % 1.0 < mutationProbability)
            child1 = Mutate(child1)
            child2 = Mutate(child1)
        end
        if(length(newPopulation) < populationSize)
            push!(newPopulation, child1)
        end
        if(length(newPopulation) < populationSize)
            push!(newPopulation, child2)
        end
    end
    return newPopulation
end

#Function generates random population
#Input:
#   availableCities : an array of cities that salesman should visit
#   populationSize : size of population
#Output:
#   population : newly generated population
@everywhere function GenerateRandomPopulation(availableCities::Array{City, 1}, populationSize::Int64)
    population::Array{Array{City, 1}, 1} = []
    for i = 1 : populationSize
        push!(population, GetRandomTour(availableCities, length(availableCities)))
    end
    return population
end

#Function generates random tour
#Input:
#   availableCities : an array of cities that salesman should visit
#   tourLength : length that new tour should have
#Output:
#   tour : newly generated tour
@everywhere function GetRandomTour(availableCities, tourLength)
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

#Function determines if city is alredy in given tour
#Input:
#   city : city to be checked
#   path : a tour to be checked with
#Output:
#   true if city already is in tour, false otherwise
@everywhere function isCityInPath(city::City, path)
    for i = 1 : length(path)
        if(path[i] == city)
            return true;
        end
    end
    return false;
end

#Function determines the distance between two cities (Euclidean)
#Input:
#   p : first city
#   q : second city
#Output:
#   Calculated distance
@everywhere function GetLength(p::City, q::City)
    dx = p.x - q.x
    dy = p.y - q.y
    return floor(sqrt(dx*dx + dy*dy))
end

#Function determines if city is alredy in given child
#Input:
#   city : city to be checked
#   child : a tour to be checked with
#Output:
#   true if city already is in tour, false otherwise
@everywhere function IsAlreadyInChild(city::Int64, child::Array{City, 1})
    for i in 1:length(child)
        if(child[i].id == city)
            return true
        end
    end
    return false
end

#Function inits the child
#Input:
#  length : length of tour
#Output:
#   child : a tour that consist of cities with id -1
@everywhere function InitChild(length::Int64)
    unknownCity::City = City(-1, -1, -1)
    child::Array{City, 1} = []
    for i in 1:length
        push!(child, unknownCity)
    end
    return child
end
