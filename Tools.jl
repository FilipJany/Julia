#Module provides several tools(methods) for Iterative TSP algorithm
include("Constraints.jl")

export greedyTSP, SA

#Genetic algorithm begin

function GeneticAlgorithm(population, steps::Int64)
    i = 0
    mutationBorder = 0.5
    while i < steps
        newPopulation = PreparePopulation(population)
        for j = 1 : length(population)
           x = Selection(population)
           y = Selection(population)
           child = Reproduce(x, y)
           if(GetRandomNumber() < mutationBorder)
               child = Mutate(child)
            end
            push!(newPopulation, child)
        end
    end
end

function Reproduce(x, y)
    cutPoint = GetRandomIntNumber(length(x))
    child = Int64[]
    for i = 1 : cutPoint
        push!(child, x[i])
    end
    for i = cutPoint : length(y)
        push!(child, y[i])
    end
    return child
end

function Mutate(child)
    mutationPoint = GetRandomIntNumber(length(child) - 1)
    temp = child[mutationPoint]
    child[mutationPoint] = child[mutationPoint + 1]
    child[mutationPoint + 1] = temp
    return child
end

function InitPopulation(pSize::Int64, connections)
    tours = Connections[]
    path = getStartingPath(startCity, length(connections), connections)
    for i = 1 : pSize
        push!(tours, computeNextPath(path, true)
    end
    return tours
end

function GetRandomIntNumber(range::Int64)
    return (rand(Int64) % (range - 1)) + 1
end

function GetRandomFloatNumber(range::Int64)
    return rand(Float64) % range
end

function GetIdListFromCities(cityList)
    idList = Int64[]
    for i = 1 : length(cityList)
        push!(idList, cityList.id)
    end
end

function PreparePopulation(population)
    startCityIdList = GetIdListFromCities(population)
    populationInt = [startCityIdList]
    return populationInt
end


#Genetic algo stop

#SimAnn begin

function SA(startCity::City, connections)
    print("Starting city : ")
    print(startCity)
    print("\n")
    sPath = getStartingPath(startCity, length(connections), connections)
    print("Starting Path : ")
    print(sPath)
    print("\n")
    nPath = City[]
    iter = -1
    P = 0.0
    alpha = 0.999
    temp = 600.0
    eps = 0.001
    delta = 0.0
    dist = getTotalDistance(startCity, sPath)
    while (temp > eps)
        iter = iter + 1
        nPath = computeNextPath(sPath, false)
        delta = getTotalDistance(startCity, nPath) - dist
        if(delta < 0)
            sPath = nPath
            dist = delta + dist
        else
            P = rand(Float64)
            if(P < exp(-delta/temp))
                sPath = nPath
                dist = delta + dist
            end
        end
        temp = temp * alpha
        if(iter % 100 == 0)
            print("Iteration : ")
            print(iter)
            print(" distance : ")
            print(dist)
            print("\n")
        end
    end
    print("Iteration : ")
    print(iter)
    print(" distance : ")
    print(dist)
    print("\n")
    print("Path : ")
    print(sPath)
    print("\n")
    return dist
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
        print("Distance from: ")
        print(u)
        print(" to: ")
        print(v)
        print(" = ")
        println(getLength(u, v))
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
