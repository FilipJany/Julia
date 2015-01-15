#File: SimulatedAnnealing.jl
#Author: Filip Jany 194208 Wrocław University of Technology
#Purpose: Engeneering Thesis (2015)
#Description: Implementation of Simulated Annealing algorithm and it's helping functions
include("../Other/Constraints.jl")

export SA

#Parallel implementation of iterative Simulated Annealing algorithm.
#Input:
#       aCities: an array of cities to visit
#       alpha : cooling rate value
#       startTemp : starting temperature
#       iterations : number of iterations to perform
#Output:
#   Best(with smallest cost) of calculated tours
function SA(aCities::Array{City, 1}, alpha::Float64, startTemp::Float64, iterations::Int64)
    T = startTemp
    x0 = aCities
    xopt = x0
    iter::Int64 = 0
    workArray::Array{Int64, 1} = ManageWork(iterations, nprocs()-1)
    results = []
    @sync begin
        for p in 1:nprocs()
            if p != myid() || nprocs() == 1
                @async begin
                    if(length(workArray) <= 2)
                        push!(results, remotecall_fetch(p, WorkerTask, x0, T, alpha, workArray[p-1], workArray[p]))
                    else
                        push!(results, remotecall_fetch(p, WorkerTask, x0, T, alpha, workArray[p-1], workArray[p]))
                    end
                end
            end
        end
    end
    return GetBest(results)
end

#Serial implementation of iterative Simulated Annealing algorithm.
#Input:
#       aCities: an array of cities to visit
#       alpha : cooling rate value
#       startTemp : starting temperature
#       iterations : number of iterations to perform
#Output:
#   bestTour : best(with smallest cost) of calculated tours
function SASerial(aCities::Array{City, 1}, alpha::Float64, startTemp::Float64, iterations::Int64)
    T = startTemp
    x0 = aCities
    xopt = x0
    iter::Int64 = 0
    bestTour::Array{City, 1} = x0
    for i in 1:iterations
        x = Select(x0)
        if(rand(0.0:1.0) < P(x0, x, T))
            x0 = x
            if(Cost(x0) < Cost(bestTour))
                bestTour = x0
            end
        end
        T *= alpha
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

#Function specifies a worker task.
#Input:
#       tour: an array of cities to visit
#       T : current temperature
#       alpha : cooling rate
#       fromIter: number of iteration to start with
#       toIter: number of iteration to end with
#Output:
#   bestTour: best found tour
@everywhere function WorkerTask(tour::Array{City, 1}, T::Float64, alpha::Float64, fromIter::Int64, toIter::Int64)
    bestTour::Array{City, 1} = tour
    x0 = tour
    temp::Float64 = T
    for i in fromIter:toIter
        x = Select(x0)
        if(rand(0.0:1.0) < P(x0, x, temp))
            x0 = x
            if(Cost(x0) < Cost(bestTour))
                bestTour = x0
            end
        end
        temp *= alpha
    end
    return bestTour
end

#Function is responsible for generating bounds of work for given processes.
#Input:
#       cities: an array of cities to visit
#       np: number of processes
#Output:
#   points: array containing bounds of work for each process
@everywhere function ManageWork(iterations::Int64, np::Int64)
    workArray::Array{Int64, 1} = [0]
    itersPerProcess::Int64 = int(iterations/np)
    current::Int64 = 1
    for i in 1:np-1
        current = itersPerProcess*i
        push!(workArray, current)
    end
    if(current < iterations)
        push!(workArray, iterations)
    end
    return workArray
end

#Function calculates the acceptance function
#Input
#   x0 : previous solution
#   x : current solution
#   T : temperature
#Output
#   Proper value
@everywhere  function P(x0, x, T)
    value::Float64 = (Cost(x) - Cost(x0)) / T
    return min(1, exp(-value))
end

#Function selects two random cities to swap
#Input
#   tour : an array with cities
#Output
#   New - swapped tour
@everywhere function Select(tour::Array{City, 1})
    i::Int64 = abs(rand(Int64) % length(tour)) + 1
    j::Int64 = abs(rand(Int64) % length(tour)) + 1
    while(i == j)
        j = abs(rand(Int64) % length(tour)) + 1
    end
    return Swap(tour, i, j)
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

#Function calculates cost of given tour.
#Input:
#   tour: an array of cities
#Output:
#   summ: cost of a tour
@everywhere function Distance(p::City, q::City)
    dx::Float64 = q.x - p.x
    dy::Float64 = q.y - p.y
    #println("DISTANCE: dx:", dx, " dy: ", dy)
    return sqrt(dx^2 + dy^2)
end
