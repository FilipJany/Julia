include("TO.jl")
include("Task.jl")
include("../Other/Constraints.jl")

res = @time summElems()

for i in 1:length(res)
    print(res[i], "->")
end
print("\n", Cost(res), "\n")

