#File: RunGeneticAlgorithm.jl
#Author: Filip Jany 194208 Wrocław University of Technology
#Purpose: Engeneering Thesis (2015)
#Description: A script that helps to process parameters given while execution from terminal.
#   Its other task is to choose which version of algorithm should be executed.
include("GeneticAlgorithm.jl")
include("../Other/Tools.jl")

parametersNumber = 4

if(length(ARGS) != parametersNumber)
    println("Usage:")
    println("Please provide ", parametersNumber, " parameters:\n",
            "\tfirst for psth to text file with TSP data,\n",
            "\tsecond for population size,\n",
            "\tthird for iterations number,\n",
            "\tfourth for mutation probability.")
    println("Example: julia -p 8 RunSimulatedAnnealing.jl ../TSPData.txt 100 10000 0.75")
    println("Execution terminated")
else
    dataFile = ARGS[1]
    populationSize = int(ARGS[2])
    iterations = int(ARGS[3])
    mutationProbability = float(ARGS[4])

    if(nprocs() > 1 && isfile(dataFile) && iterations > 1)
        println("Running Simulated Annealing in parallel mode with: ", nprocs()-1, " slaves and one master.")
        cities = CreateCityListFromFile(dataFile)
        println("Starting algorithm, please wait for results.")
        result = @time GeneticAlgorithm(cities, iterations, populationSize, mutationProbability)
        println("Calculated tour:")
        for i in 1 : length(result)
            print(result[i].id, "->")
        end
        print("\n")
        println("Result cost: ", Cost(result))
    elseif(isfile(dataFile) && iterations > 1)
        println("Running Simulated Annealing in serial mode.")
        cities = CreateCityListFromFile(dataFile)
        result = @time GeneticAlgorithmSerial(cities, iterations, populationSize, mutationProbability)
        println("Calculated tour:")
        for i in 1 : length(result)
            print(result[i].id, "->")
        end
        print("\n")
        println("Result cost: ", Cost(result))
    elseif(iterations <= 1)
        println("To small iterations number!")
    elseif(!isfile(dataFile))
        println("File broken or does not exist!")
    else
        println("Unknown error!")
    end
end
