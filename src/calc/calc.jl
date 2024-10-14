module calc
"""
Provides tools for calculating meteorological elements
"""

include("basefun.jl")
using Unitful

export c
function c(num::typeof(1u"m/s"))
    return num
end

"""
T  在函数里转换单位
"""
# dewpoint
function dewpoint_from_rh(T_air::Unitful.Temperature, RH::Real)
    T_air = uconvert(u"°C", T_air)
    T_air = ustrip(T_air)
    logval = (7.5 * T_air) / (237.3 + T_air)
    e_s = 6.11 * 10^logval * u"hPa"
    e = e_s * RH
    while e < e_s
        T_air = T_air - 0.001
        logval = (7.5 * T_air) / (237.3 + T_air)
        e_s = 6.11 * 10^logval * u"hPa"
    end
    return uconvert.(u"K", T_air * u"°C")
end
function dewpoint_from_rh(T_air, RH)
    dewpoint_from_rh.(T_air, RH)
end

function dewpoint_from_sh(pressure, q)
    q = uconvert.(u"g/kg", q)
    pressure = uconvert.(u"hPa", pressure)
    w = q ./ (1 .- q)
    e = w .* pressure ./ (0.622 .+ w)
    return dewpoint.(e)
end

function dewpoint(vapor_pressure)
    vapor_pressure = uconvert.(u"hPa", vapor_pressure)
    logvalue = log10.(ustrip.(vapor_pressure) ./ 6.112)
    dewPoint = (237.3 .* logvalue) ./ (7.5 .- logvalue) .* u"°C"
    return uconvert.(u"K", dewPoint)
end

function vapor_pressure_from_dewpoint(dewPoint)
    dewPoint = uconvert.(u"°C", dewPoint)
    logvalue = (7.5 .* ustrip(dewPoint)) ./ (237.3 .+ ustrip(dewPoint))
    vapor_pressure = 6.112 .* [10] .^ logvalue .* u"hPa"
    return vapor_pressure
end


# wind_chill_index
function wind_chill_index(wind_speed, T_air)
    wind_speed = ustrip(uconvert.(u"m/s", wind_speed))
    T_air = ustrip(uconvert.(u"°C", T_air))
    WCI_up = (10 .* sqrt.(wind_speed) .- wind_speed .+ 10.5)
    WCI_down = (33 .- T_air)
    return (WCI_up .* WCI_down) * u"kg*cal/m^2/h"
end

end
