#Module provides several tools(methods) for Iterative TSP algorithm
include("Constraints.jl")

export greedyTSP, SA, GeneticAlgorithm

#Genetic algorithm begin

function GeneticAlgorithm(availableCities, startCity, maxSteps::Int64, populationSize::Int64)
    i = 0
    bestDistance = 999.0
    bestIteration = 0
    bestTour = City[]
    print("Initial distance : ")
    print(1/FitnessFunction(availableCities))
    population = InitPopulation(availableCities, populationSize, startCity)
    while i <= maxSteps
        population = Evolve(population, availableCities[startCity])
        print("\niteration: ", i, " distance: ", GetTourDistance(GetFittestFromPopulation(population)))
        if(IsAnyBetterThanBest(bestDistance, population))
            bestDistance = GetTourDistance(GetFittestFromPopulation(population))
            bestIteration = i
            bestTour = GetFittestFromPopulation(population)
        end
        i = i + 1
    end
    print("\n")
    print("Shortest path in ", maxSteps, " iteration: ")

    print("\nShortest path after ", maxSteps, " iterations: ")
    for i = 1 : length(bestTour) - 1
        print(bestTour[i].id, "-> ")
    end
    print(bestTour[length(bestTour)].id)
    print(". \nDistance : ", bestDistance, " achived in iteration #", bestIteration, ".\n")

end

function IsAnyBetterThanBest(best, recentPopulation)
    recentBest = GetTourDistance(GetFittestFromPopulation(recentPopulation))
    if(recentBest < best)
        return true;
    end
    return false;
end

function Evolve(population, startCity)
    mutationProbability = 0.8
    newPopulation = {}
    for i = 1 : length(population)
        parent1 = Selection(population)
        parent2 = Selection(population)
        child = Reproduce(parent1, parent2, startCity)
        if(GetRandomFloatNumber(1.0) < mutationProbability)
            child = Mutate(child)
        end
        push!(newPopulation, child)
    end
    return newPopulation
end

function InitPopulation(availableCities, pSize::Int64, startCity)
    path = GetRandomRoute(startCity, availableCities)
    tours = {path}
    for i = 1 : pSize
        push!(tours, GetRandomRoute(startCity, availableCities))
    end
    return tours
end

function GetFittestFromPopulation(population)
    fittest = population[1]
    for i = 2 : length(population)
        if((1/FitnessFunction(fittest)) > (1/FitnessFunction(population[i])))
            fittest = population[i]
        end
    end
    return fittest
end

function GetRandomRoute(startCity, cities)
    path = City[]
    push!(path, cities[startCity])
    i = 1
    while i < length(cities)
        randomIndex = abs(rand(Int64)%(length(cities) - 1)) + 2
        if(!isCityInPath(cities[randomIndex], path))
            push!(path, cities[randomIndex])
            i = i + 1
        end
    end
    return path
end

function Selection(population)
    k = GetRandomIntNumber(length(population))
    tournament = {}
    for i = 1 : k
        push!(tournament, population[GetRandomIntNumber(length(population))])
    end
    return GetFittestFromPopulation(tournament)
end

function GetCityIndexInPopulation(population, tour)
    for i = 1 : length(population)
        if(population[i] == tour)
            return i
        end
    end
    return -1
end

function GetCityFromPopulation(population, cityId)
    for i = 1 : length(population)
        if(population[i].id == cityId)
            return population[i]
        end
    end
    return City(0, 0, 0)
end

function GetTourDistance(tour)
    totalDistance = 0.0
    for i = 2 : length(tour)
        totalDistance = totalDistance + getLength(tour[i-1], tour[i])
    end
    totalDistance = totalDistance + getLength(tour[length(tour)], tour[1])
    return totalDistance
end

function FitnessFunction(tour)
    return 1 / GetTourDistance(tour)
end

function Reproduce(x, y, startCity)
    child = {startCity}
    startPos = GetRandomIntNumber(length(x))
    endPos = GetRandomIntNumber(length(x))
    for i = 2 : length(x)
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

function Mutate(child)
    mutationPoint1 = GetRandomIntNumber(length(child)-1)+1
    mutationPoint2 = GetRandomIntNumber(length(child)-1)+1
    temp = child[mutationPoint1]
    child[mutationPoint1] = child[mutationPoint2]
    child[mutationPoint2] = temp
    return child
end

function GetRandomIntNumber(range::Int64)
    return abs(rand(Int64) % (range - 1)) + 1
end

function GetRandomFloatNumber(range)
    return rand(Float64) % range
end

function GetIdListFromCities(cityList)
    idList = Int64[]
    for i = 1 : length(cityList)
        push!(idList, cityList[i].id)
    end
    return idList
end

function PreparePopulation(population)
    cityIdList = {GetIdListFromCities(population[1])}
    for i = 2 : length(population)
        push!(cityIdList, GetIdListFromCities(population[i]))
    end
    return cityIdList
end


#Genetic algo stop

#SimAnn begin

function SA(startCity::City, connections)
    temp = 1000.0
    eps = 0.0001
    coolingRate = 0.997
    useNew = false
    currentSolution = GenerateRandomTour(startCity, connections)
    workingSolution = CopyTour(currentSolution)
    bestSolution = CopyTour(currentSolution)
    while temp > eps
        random1 = GetRandomIntNumber(length(currentSolution)) + 1
        random2 = GetRandomIntNumber(length(currentSolution)) + 1
        t = workingSolution[random1]
        workingSolution[random1] = workingSolution[random2]
        workingSolution[random2] = t
        energy = GetTourDistance(workingSolution)
        if(energy <= GetTourDistance(currentSolution))
            useNew = true
        else
            randomFloat = rand(Float64) % 1.0
            delta = energy - GetTourDistance(currentSolution)
            calc = exp(-delta / temp)
            if(calc > randomFloat)
                useNew = true
            end
        end
        if(useNew == true)
            useNew = false
            currentSolution = CopyTour(workingSolution)
            if(GetTourDistance(currentSolution) < GetTourDistance(bestSolution))
                bestSolution = CopyTour(currentSolution)
            end
        else
            bestSolution = CopyTour(currentSolution)
        end
        #print("\nCurrent sol energy: ", GetTourDistance(currentSolution))
        #print("\nWorking sol energy: ", GetTourDistance(workingSolution))
        #print("\nBest sol energy: ", GetTourDistance(bestSolution))
        temp *= coolingRate
    end
    PrintTour(bestSolution, "Best Solution: ")
    print("Best Distance: ", GetTourDistance(bestSolution), "\n")
end

function PrintTour(tour, descr)
    print("\n", descr)
    for i = 1 : length(tour) - 1
        print(tour[i].id, "-> ")
    end
    print(tour[length(tour)].id, "\n")
end

function Fitness(x, xo, temp)
    delta = (x - xo) / temp
    delta *= -1
    return min(1, exp(delta))
end

function CopyTour(tour)
    newTour = {}
    for i = 1 : length(tour)
        push!(newTour, tour[i])
    end
    return newTour
end

function GetNewTour(tour)
    newTour = {}
    for i = 1 : length(tour)
        push!(newTour, tour[i])
    end
    tourPos1 = GetRandomIntNumber(length(tour)) + 1
    tourPos2 = GetRandomIntNumber(length(tour)) + 1

    temp = newTour[tourPos1]
    newTour[tourPos1] = newTour[tourPos2]
    newTour[tourPos2] = temp
    return newTour
end


function GenerateRandomTour(startCity, connections)
    newTour = {}
    push!(newTour, startCity)
    i = 2
    while i <= length(connections)
        randomIndex = GetRandomIntNumber(length(connections)) + 1
        randomCity = connections[randomIndex]
        if(!isCityInPath(randomCity, newTour))
            push!(newTour, randomCity)
            i += 1
        end
    end
    return newTour
end



function getStartingPath(startCity::City, pathLength::Int64, connections)
    path = City[]
    push!(path, startCity)
    i = 2
    while (i <= pathLength)
        city = connections[abs(rand(Int64) % pathLength) + 1]
        if(!isCityInPath(city, path))
            push!(path, city)
            i = i + 1
        end
    end
    return path
end

function getTotalDistance(startCity::City, path)
    dist = 0
    for i = 1 : length(path) - 1
        dist = dist + getLength(path[i], path[i+1])
    end
    lastCity = path[length(path)]
    dist = dist + getLength(path[length(path)], path[getCityIndex(startCity, path)])
    return dist
end

function computeNextPath(path, fullyRandom::Bool)
    nPath = path
    for i = 1 : length(path)
        if (fullyRandom)
            x = abs(rand(Int64) % (length(path) - 1)) + 1
            y = abs(rand(Int64) % (length(path) - 1)) + 1
        else
            x = abs(rand(Int64) % (length(path) - 1)) + 2
            y = abs(rand(Int64) % (length(path) - 1)) + 2
        end

        temp = nPath[x]
        nPath[x] = nPath[y]
        nPath[y] = temp
    end
    return nPath
end

#SimAnn End
#Greedy begin
function greedyTSP(startCity::City, connections)
    visited = City[]
    notVisited = connections
    totalLength = 0
    index = getCityIndex(startCity, connections)
    if(index > 0)
        splice!(notVisited, index)
    else
        return "Error: Index of start city not found..."
    end
    push!(visited, startCity)
    u = startCity
    while(length(notVisited) > 0)
        v = closestCity(u, connections)
        totalLength = totalLength + getLength(u, v)
        push!(visited, v)
        index = getCityIndex(v, connections)
        if(index > 0)
            splice!(notVisited, index)
        else
            return "Error: Index of selected city not found..."
        end
        u = v
    end
    print("Total length: ")
    println(totalLength)
    return visited
end

function closestCity(fromCity::City, connections)
    min = Inf
    city = City(0, 0, 0)
    for i = 1 : length(connections)
        if(getLength(fromCity, connections[i]) < min)
            min = getLength(fromCity, connections[i])
            city = connections[i]
        end
    end
    return city
end

function getCityIndex(city::City, connections)
    for i = 1 : length(connections)
        if(connections[i] == city)
            return i
        end
    end
    return -1
end

function getLength(from::City, to::City)
    return sqrt((from.x - to.x)^2 + (from.y - to.y)^2)
end

function printLengths(connections)
    println("Printing all distances:")
    for i = 1 : length(connections)
        for j = 1 : length(connections)
            print("distance from: ")
            print(connections[i])
            print("to: ")
            print(connections[j])
            print(" = ")
            println(getLength(connections[i], connections[j]))
        end
    end
    println("***********************")
end

function generateRandomPath(startCity::Int64, pathLength::Int64)
     path = Int64[]
     push!(path, startCity)
     i = 1
     while(i < pathLength)
         cCity = abs(rand(Int64) % pathLength)
         if(!isCityInPath(cCity, path))
             push!(path, cCity)
             i = i + 1
         end
     end
     return path
end

function isCityInPath(city::City, path)
    for i = 1 : length(path)
        if(path[i] == city)
            return true;
        end
    end
    return false;
end

function P(key, value)
    print("\n********\n")
    print(key)
    print(": ")
    print(value)
    print("\n********\n")
end
