include("../Other/Constraints.jl")

export SA

function SA(aCities::Array{City, 1}, alpha::Float64, startTemp::Float64, minTemp::Float64, iterations::Int64)
    file = open("fittness.txt", "w")
    T = startTemp
    x0 = aCities
    xopt = x0
    iter::Int64 = 0
    while(T >= minTemp && iter <= iterations)
       x = Select(x0)
       write(file, string(T, " ", Cost(x), "\n"))
       if(Cost(x) < Cost(xopt))
           xopt = x
       else
           random = rand(Float64) % 1.0
           acc::Float64 = P(x0, x, T)
           if(random < acc)
               x0 = x
               if(Cost(x0) < Cost(xopt))
                   xopt = x0
               end
           end
       end
       T *= alpha
       iter += 1
       #println("T: ", T, " iteration: ", iter)
    end
    close(file)
    return xopt
end

function P(x0, x, T)
    value::Float64 = (Cost(x) - Cost(x0)) / T
    return min(1, exp(-value))
end

function Select(tour::Array{City, 1})
    i::Int64 = abs(rand(Int64) % length(tour)) + 1
    j::Int64 = abs(rand(Int64) % length(tour)) + 1
    while(i == j)
        j = abs(rand(Int64) % length(tour)) + 1
    end
    return Swap(tour, i, j)
end

function Swap(tour::Array{City, 1}, city1Index::Int64, city2Index::Int64)
    #println("Will swap city at Index1: ", city1Index, " and city at: ", city2Index)
    #for i in 1:length(tour)
    #    print(tour[i].id, " ")
    #end
    #print("\n")
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
    #for i in 1:length(newTour)
    #    print(newTour[i].id, " ")
    #end
    #print("\n")
    return newTour
end

function Cost(tour::Array{City, 1})
    summ::Float64 = 0.0
    for i = 1 : length(tour)-1
        summ += Distance(tour[i], tour[i+1])
    end
    summ += Distance(tour[length(tour)], tour[1])
    return summ
end

function CostWithLog(tour::Array{City, 1})
    summ::Float64 = 0.0
    for i = 1 : length(tour)-1
        summ += Distance(tour[i], tour[i+1])
        println("Distance from: ")
        println("\t", tour[i].id)
        println("\t", "X: ", tour[i].x, "Y: ", tour[i].y)
        println(" to: ")
        println("\t",tour[i + 1].id)
        println("\t", "X: ", tour[i+1].x, "Y: ", tour[i+1].y)
        println(" = ", Distance(tour[i], tour[i+1]), ". ", summ, " in total.")
    end
    summ += Distance(tour[length(tour)], tour[1])
    return summ
end

function Distance(p::City, q::City)
    dx::Float64 = q.x - p.x
    dy::Float64 = q.y - p.y
    #println("DISTANCE: dx:", dx, " dy: ", dy)
    return sqrt(dx^2 + dy^2)
end

