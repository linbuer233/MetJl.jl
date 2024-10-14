# Tools that provide basic meteorological elements and basic functions
module basefun

using Unitful

const t0 = 288.15 * u"K"
const p0 = 1023.5 * u"hPa"
const Γ = 6.5 * u"K" / u"km"

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
