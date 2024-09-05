using MetJl
using Test

@testset "MetJl.jl" begin
    @test MetJl.hello_world() == "Hello, World!"
end
