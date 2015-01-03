include("Tools.jl")

cities = CreateCityListFromFile("../Other/Problems/berlin52.tsp")
for i in 1:length(cities)
    println("ID: ", cities[i].id, " X: ", cities[i].x, " Y: ", cities[i].y)
end
