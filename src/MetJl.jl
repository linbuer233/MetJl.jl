module MetJl

include("./calc/calc.jl")
export calc
"""
hi = hello_world()
A simple function to return "Hello, World!"
"""
function hello_world()
    return "Hello, World!"
end

function compute(num)
    return num * num
end

end
