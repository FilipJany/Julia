#Module provides several tools(methods) for Iterative TSP algorithm
include("Constraints.jl")

export ReadFromFile, CreateCityListFromFile

function ReadFromFile(fileName::String)
    f = open(fileName)
    lines = readlines(f)
    print(length(lines), " have been successfully read.\n")
    return lines
end

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

function GetCityFromString(entry::String)
    array = split(entry)
    return City(float(convert(ASCIIString, array[2])), float(convert(ASCIIString, array[3])), parseint(convert(ASCIIString, array[1])))
end






















