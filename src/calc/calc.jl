module calc
# Provides tools for calculating meteorological elements

include("basefun.jl")
using Unitful

# T  在函数里转换单位

# dewpoint
"""
calculates the dewPoint
# Arguments:
# T_air - temperature of the air [°C or K]
# RH - relative humidity [%]
#
# Returns:
# dewpoint - dewpoint temperature [°C or K]
"""
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

function rh_from_dewpont(T_air, dewpoint)
    T_air = uconvert.(u"°C", T_air)
    dewpoint = uconvert.(u"°C", dewpoint)
    e = vapor_pressure_from_dewpoint.(dewpoint)
    e_s = vapor_pressure_from_dewpoint.(T_air)
    return e ./ e_s
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
    vapor_pressure = 6.112 .* (10) .^ logvalue .* u"hPa"
    return vapor_pressure
end

function q_from_dewpoint(pressure, dewpoint)
    pressure = uconvert.(u"hPa", pressure)
    dewpoint = uconvert.(u"K", dewpoint)
    e = vapor_pressure_from_dewpoint.(dewpoint)
    return 0.622 .* e ./ (pressure .- 0.378 .* e)
end

function q_from_vapor_pressure(pressure, vapor_pressure)
    pressure = uconvert.(u"hPa", pressure)
    vapor_pressure = uconvert.(u"hPa", vapor_pressure)
    e = vapor_pressure
    return 0.622 .* e ./ (pressure .- 0.378 .* e)
end

function LCL(
    pressure::Unitful.Pressure,
    T_air::Unitful.Temperature,
    dewpoint::Unitful.Temperature,
)
    pressure = uconvert.(u"hPa", pressure)
    T_air = uconvert.(u"K", T_air)
    dewpoint = uconvert.(u"K", dewpoint)
    e_z = vapor_pressure_from_dewpoint.(dewpoint)
    θ = theta.(T_air, pressure)
    q_z = q_from_vapor_pressure(pressure, e_z)
    P_L = 300 * u"hPa"
    local T_L = 1
    while true
        T_L = θ .* (P_L ./ (1000 * u"hPa")) .^ 0.286
        e_L = vapor_pressure_from_dewpoint.(T_L)
        q_L = q_from_vapor_pressure.(P_L, e_L)
        if q_L >= q_z
            break
        end
        P_L += 1 * u"hPa"
    end
    return T_L, P_L
end

function LCL(pressure, T_air, dewpoint)
    return LCL.(pressure, T_air, dewpoint)
end

function theta(T_air, Pressure)
    T_air = ustrip(uconvert(u"K", T_air))
    Pressure = ustrip(uconvert(u"hPa", Pressure))
    return T_air .* (1000 ./ Pressure) .^ (0.286) .* u"K"
end

# height and pressure values
function add_height_to_pressure(pressure, height)
    pressure = uconvert.(u"hPa", pressure)
    height = uconvert.(u"m", height)
    pressure_level_height = basefun.pressure_to_height_std.(pressure)
    return basefun.height_to_pressure_std(pressure_level_height .+ height)
end

function add_pressure_to_height(height, pressure)
    height = uconvert.(u"m", height)
    pressure = uconvert.(u"hPa", pressure)
    height_level_pressure = basefun.height_to_pressure_std.(height)
    return basefun.pressure_to_height_std(height_level_pressure .- pressure)
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
