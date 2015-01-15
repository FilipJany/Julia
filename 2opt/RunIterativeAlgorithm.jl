#File: RunIterativeAlgorithm.jl
#Author: Filip Jany 194208 Wrocław University of Technology
#Purpose: Engeneering Thesis (2015)
#Description: A script that helps to process parameters given while execution from terminal.
#   Its other task is to choose which version of algorithm should be executed.
include("2opt.jl")
include("../Other/Tools.jl")

parametersNumber = 2

if(length(ARGS) != parametersNumber)
    println("Usage:")
    println("Please provide ", parametersNumber, " parameters:\n",
            "\tfirst for path to text file with TSP data,\n",
            "\tsecond for maximal consecutive minimal repetitions (greater than 1).")
    println("Example: julia -p 8 RunIterativeAlgorithm.jl ../TSPData.txt 10")
    println("Execution terminated")
else
    maxConsecutiveValues = int(ARGS[2])
    dataFile = ARGS[1]
    if(nprocs() > 1 && isfile(dataFile) && maxConsecutiveValues > 1)
        println("Running 2opt in parallel mode with: ", nprocs()-1, " slaves and one master.")
        cities = CreateCityListFromFile(dataFile)
        println("Starting algorithm, please wait for results.")
        result = TwoOpt(cities, maxConsecutiveValues, ManageWork(cities, nprocs()-1))
        println("Calculated tour:")
        for i in 1 : length(result)
            print(result[i].id, "->")
        end
        print("\n")
        println("Result cost: ", Cost(result))
    elseif(isfile(dataFile) && isfile(dataFile) && maxConsecutiveValues > 1)
        println("Running 2opt in serial mode.")
        cities = CreateCityListFromFile(dataFile)
        result = TwoOptSerial(cities, maxConsecutiveValues)
        println("Calculated tour:")
        for i in 1 : length(result)
            print(result[i].id, "->")
        end
        print("\n")
        println("Result cost: ", Cost(result))
    elseif(maxConsecutiveValues <= 1)
        println("To small consecutive repetitions number!")
    elseif(!isfile(dataFile))
        println("File broken or does not exist!")
    else
        println("Unknown error!")
    end
end
