#File: RunIterativeAlgorithm.jl
#Author: Filip Jany 194208 Wrocław University of Technology
#Purpose: Engeneering Thesis (2015)
#Description: A script that contains a definition of program's constraints
export City

#Declader type that represnts a city in tour it consists of:
#   x - x coordinate
#   y - y coordinate
#   id - city id
@everywhere type City
    x
    y
    id
end


