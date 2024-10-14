using MetJl
using Unitful
using Test

@testset "MetJl.jl" begin
    @test MetJl.hello_world() == "Hello, World!"
    @test try
        calc.dewpoint_from_rh(20u"K", 0.5)
        true
    catch e
        false
    end
    @test try
        calc.dewpoint_from_rh([20, 20, 20] * u"K", [0.5, 0.5, 0.5])
        true
    catch e
        false
    end
    @test try
        calc.dewpoint_from_sh([1000, 1000, 1000] * u"hPa", [10, 10, 10]u"g/kg")
        true
    catch e
        false
    end
    @test try
        calc.dewpoint_from_sh(1000 * u"hPa", 10u"g/kg")
        true
    catch e
        false
    end
    @test try
        calc.dewpoint(500 * u"hPa")
        true
    catch e
        false
    end
    @test try
        calc.vapor_pressure_from_dewpoint(300 * u"K")
        true
    catch e
        false
    end
    @test try
        calc.wind_chill_index(10u"m/s", 300 * u"K")
        true
    catch e
        false
    end
end
