include("../Other/Constraints.jl")

@everywhere type Coordinates
    x::Int64
    y::Int64
end

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

@everywhere function GetBest(results)
    best::Array{City, 1} = results[1]
    for i in 1:length(results)
        if(Cost(results[i]) < Cost(best))
            best = results[i]
        end
    end
    return best
end

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
end

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

@everywhere function Cost(tour::Array{City, 1})
    summ::Float64 = 0.0
    for i = 1 : length(tour)-1
        summ += Distance(tour[i], tour[i+1])
    end
    summ += Distance(tour[length(tour)], tour[1])
    return summ
end

@everywhere function Distance(p::City, q::City)
    dx = p.x - q.x
    dy = p.y - q.y
    return floor(sqrt(dx*dx + dy*dy))
end
