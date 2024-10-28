# Tools that provide basic meteorological elements and basic functions
module basefun

using Unitful


# Constants
const g = 9.81 * u"m/s^2"
const L = 287.1 * u"J/(kg*K)"
const M = 0.0289652 * u"kg/mol"
const R = 8.314462618 * u"J/(mol*K)"
const Md = 0.02896546 * u"kg/mol"
const Rd = R / Md
const t0 = 288.15 * u"K"
const p0 = 1023.5 * u"hPa"
const Γ = 6.5 * u"K" / u"km"
const L_vaper = 2.5 * 10^6 * u"J/kg"
const γ_d = 9.8 * u"K" / u"km"
const ϵ = 0.62197 * u"kg/kg"
const Cp = 1004 * u"J/kg/K"  #干空气定压比热
const Cpd = 1005 * u"J/kg/K" # 湿空气定压比热
const Cw = 4218 * u"J/K/kg"

function height_to_pressure_std(height)
    height = uconvert(u"m", height)
    return p0 * (1 - (Γ / t0) * height)^(g / (Rd * Γ))
end

function pressure_to_height_std(pressure)
    pressure = uconvert(u"hPa", pressure)
    downval = (t0 / Γ)
    logval = (1 - (pressure / p0)^(Rd * Γ / g))
    return downval * logval
end


function wind_speed(u, v)
    return sqrt(u^2 + v^2)
end

function wind_direction(u, v)
    if u >= 0 && v >= 0
        return atand(u / v) * u"°"
    end
    if u > 0 && v <= 0
        return 90 - atand(v / u) * u"°"
    end
    if u <= 0 && v < 0
        return 180 + atand(u / v) * u"°"
    end
    if u < 0 && v >= 0
        return 270 - atand(v / u) * u"°"
    end
end

function wind_components(speed, direction)
    u = speed * sind(direction)
    v = speed * cosd(direction)
    return u, v
end

end
