#File: RunIterativeAlgorithm.jl
#Author: Filip Jany 194208 Wrocław University of Technology
#Purpose: Engeneering Thesis (2015)
#Module provides several tools(methods) for TSP algorithms, i.e reading from file
include("Constraints.jl")

export ReadFromFile, CreateCityListFromFile

#Function is responsible for reading text file lines
#Input
#   fileName : path to a file
#Output
#   lines : array of file's lines
function ReadFromFile(fileName::String)
    f = open(fileName)
    lines = readlines(f)
    print(length(lines), " have been successfully read.\n")
    return lines
end

#Function is responsible for creating available cities list.
#Input
#   fileName : path to a file
#Output
#   cities : an array of available cities
function CreateCityListFromFile(fileName::String)
    regex = r"^{.+\d+}|\d+.+\d+.+\d+"
    lines = ReadFromFile(fileName)
    cities = City[]
    for i = 1 : length(lines)
        currentEntry = lines[i]
        currentMatch = match(regex, currentEntry)
        if(currentMatch != nothing)
            push!(cities, GetCityFromString(currentMatch.match))
        end
    end
    return cities
end

#Function is responsible for creating a single city entry from single line reded from file
#Input
#   entry : single line readed from file
#Output
#   A City type with proper values
function GetCityFromString(entry::String)
    array = split(entry)
    return City(float(convert(ASCIIString, array[2])), float(convert(ASCIIString, array[3])), parseint(convert(ASCIIString, array[1])))
end






















