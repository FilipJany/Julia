export Connections, City, Population, Ption, Tour

type City
    x
    y
    id
end

type Ption
    tours
end



c1 = City(-4, 4, 1)
c2 = City(-1, 7, 2)
c3 = City(4, 7, 3)
c4 = City(7, 4, 4)
c5 = City(5, -3, 5)
c6 = City(-4, -3, 6)
c7 = City(2, 3, 7)
c8 = City(1, -2, 8)
c9 = City(11, 8, 9)

Connections = City[
                c1,
                c2,
                c3,
                c4,
                c5,
                c6,
                c7,
                c8,
                c9 ]

Tour = Connections[]
Population = Ption(Tour)

