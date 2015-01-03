include("../Other/Constraints.jl")

@everywhere function TO(aCities::Array{City, 1}, maxRetries::Int64, startFrom::Int64, stopAt::Int64)
    #file = open("fittness.txt", "w")
    T = aCities
    improve::Int64 = 0
    bestT = T
    iter::Int64 = 0
    iterFile::Float64 = 0.0
    while(improve < maxRetries)
        bestDistance = Cost(T)
        for(i = startFrom:stopAt-1)
            for(j = i+1:stopAt)
                newT = Swap(T, i, j)
                newDistance = Cost(newT)
                #if(iterFile < 1000.0)
                #    stringToWrite = string(iterFile, " ", newDistance, "\n")
                #    write(file, stringToWrite)
                #    iterFile += 1.0
                #end

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
    #close(file)
    #print("Best tour:\n")
    #for i = 1 : length(bestT)
    #    print(bestT[i].id, ", ")
    #end
    #print("\n")
    return bestT
end

@everywhere function Swap(tour::Array{City, 1}, city1Index::Int64, city2Index::Int64)
    newTour::Array{City, 1} = []
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
