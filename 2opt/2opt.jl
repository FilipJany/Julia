#File: 2opt.jl
#Author: Filip Jany 194208 Wrocław University of Technology
#Purpose: Engeneering Thesis (2015)
#Description: Implementation of 2opt algorithm and it's helping functions
include("../Other/Constraints.jl")

#There was a need to declare separate type for parallel implementation of algorithm.
#Type stands for representation of swap's coordinates, i.e. Coordinates(1, 2) tells to
#   perform swap operation on cities at indexes 1 and 2.
@everywhere type Coordinates
    x::Int64
    y::Int64
end

#Parallel implementation of iterative 2-opt algorithm.
#Input:
#       cities: an array of cities to visit
#       maxRetries: maximal number of consecutive retries
#       workArray: an array with numbers, bounds of each process work
#Output:
#   bestT: an array containing shortest path
@everywhere function TwoOpt(cities::Array{City, 1}, maxRetries::Int64, workArray::Array{Int64, 1})
    T = cities
    improve::Int64 = 0
    bestT = T
    iter::Int64 = 0
    iterFile::Float64 = 0.0
    results = []
    while(improve < maxRetries)
        bestDistance = Cost(T)
        @sync begin
            for p in 1:nprocs()
                if p != myid() || nprocs() == 1
                    @async begin
                        if(length(workArray) <= 2)
                            push!(results, remotecall_fetch(p, WorkerTask, T, workArray[1], workArray[2]))
                        else
                            push!(results, remotecall_fetch(p, WorkerTask, T, workArray[p-1], workArray[p]))
                        end
                    end
                end
            end
        end
        currentlyBest = GetBest(results)
        if(Cost(currentlyBest) < bestDistance)
            bestDistance = Cost(currentlyBest)
            T = currentlyBest
            bestT = currentlyBest
            improve = 0
        end
        iter += 1
        improve += 1
    end
    return bestT
end

#Serial implementation of iterative 2-opt algorithm.
#Input:
#       cities: an array of cities to visit
#       maxRetries: maximal number of consecutive retries
#Output:
#   bestT: an array containing shortest path
@everywhere function TwoOptSerial(aCities::Array{City, 1}, maxRetries::Int64)
    T = aCities
    improve::Int64 = 0
    bestT = T
    iter::Int64 = 0
    iterFile::Float64 = 0.0
    while(improve < maxRetries)
        bestDistance = Cost(T)
        for(i = 1:length(T)-1)
            for(j = i+1:length(T))
                newT = Swap(T, i, j)
                newDistance = Cost(newT)
                if(newDistance < bestDistance)
                    bestDistance = newDistance
                    T = newT
                    bestT = newT
                    improve = 0
                end
            end
        end
        iter += 1
        improve += 1
    end
    return bestT
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

#Function specifies a worker task.
#Input:
#       tour: an array of cities to visit
#       from: number of iteration to start with
#       to: number of iteration to end with
#Output:
#   bestTour: best found tour
@everywhere function WorkerTask(tour::Array{City, 1}, from::Int64, to::Int64)
    bestTour::Array{City, 1} = tour
    for(k = from:to)
        coord = GetSwapCoordinates(k)
        newT = Swap(tour, coord.x, coord.y)
        if(Cost(newT) < Cost(bestTour))
            bestTour = newT
        end
    end
    return bestTour
end

#Function is responsible for generating bounds of work for given processes.
#Input:
#       cities: an array of cities to visit
#       np: number of processes
#Output:
#   points: array containing bounds of work for each process
@everywhere function ManageWork(cities::Array{City, 1}, np::Int64)
    n::Int64 = length(cities)
    allSwapsPerProcess::Int64 = int((n * (n-1))/(2*np))
    points::Array{Int64, 1} = [0]
    current::Int64 = 1
    for i in 1:np-1
        current = allSwapsPerProcess*i
        push!(points, current)
    end
    allSwaps = int((n * (n-1))/2)
    if(current < allSwaps)
        push!(points, allSwaps-1)
    end
    return points
end

#Function is responsible for returning coordinates of cities that should be swaped in given iteration.
#Input:
#   swapNumber: number of iteration
#Output:
#   Coordinates(x, y): Coordinates type where 'x' and 'y' are suitable city ids.
@everywhere function GetSwapCoordinates(swapNumber::Int64)
    a = 1
    b = 1
    c = -2 * swapNumber
    delta = (-b)^2 - 4*a*c
    i = int(floor(((-b) + sqrt(delta))/2*a))
    iNumber = (i * (i + 1)/2)
    iNextNumber = ((i+1) * (i + 2)/2)
    j = 1
    for k in iNumber:iNextNumber
        if(iNumber != swapNumber)
            iNumber += 1
            j += 1
        end
    end
    return Coordinates(i+2, j)
end


#Function performs a swap operation on given tour.
#Input:
#       tour:   an array of cities
#       city1Index:     index of first city to swap
#       city2Index:     index of secon city to swap
#Output:
#   newTour:    tour with swapped cities
@everywhere function Swap(tour::Array{City, 1}, city1Index::Int64, city2Index::Int64)
    newTour::Array{City, 1} = []
    #Exchange indexes if second is lesser than first
    if(city2Index < city1Index)
        for i = 1 : city2Index - 1
            push!(newTour, tour[i])
        end
        dec = 0
        for i = city2Index : city1Index
            push!(newTour, tour[city1Index - dec])
            dec += 1
        end
        for i = city1Index+1 : length(tour)
            push!(newTour, tour[i])
        end
    else
        for i = 1 : city1Index - 1
            push!(newTour, tour[i])
        end
        dec = 0
        for i = city1Index : city2Index
            push!(newTour, tour[city2Index - dec])
            dec += 1
        end
        for i = city2Index+1 : length(tour)
            push!(newTour, tour[i])
        end
    end
    return newTour
end

#Function calculates cost of given tour.
#Input:
#   tour: an array of cities
#Output:
#   summ: cost of a tour
@everywhere function Cost(tour::Array{City, 1})
    summ::Float64 = 0.0
    for i = 1 : length(tour)-1
        summ += Distance(tour[i], tour[i+1])
    end
    summ += Distance(tour[length(tour)], tour[1])
    return summ
end

#Function determines the distance between two cities (Euclidean)
#Input:
#   p : first city
#   q : second city
#Output:
#   Calculated distance
@everywhere function Distance(p::City, q::City)
    dx = p.x - q.x
    dy = p.y - q.y
    return floor(sqrt(dx*dx + dy*dy))
end
