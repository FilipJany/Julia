export Connections, City, Population, Ption, Tour

@everywhere type City
    x
    y
    id
end

@everywhere c1 = City(-4, 4, 1)
@everywhere c2 = City(-1, 7, 2)
@everywhere c3 = City(4, 7, 3)
@everywhere c4 = City(7, 4, 4)
@everywhere c5 = City(5, -3, 5)
@everywhere c6 = City(-4, -3, 6)
@everywhere c7 = City(2, 3, 7)
@everywhere c8 = City(1, -2, 8)
@everywhere c9 = City(11, 8, 9)

@everywhere Connections = City[
                c7,
                c1,
                c2,
                c5,
                c3,
                c8,
                c4,
                c6]
                #c9 ]


