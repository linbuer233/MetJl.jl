module calc
# Provides tools for calculating meteorological elements

include("basefun.jl")
using Unitful

# dewpoint
"""
dewpoint(vapor_pressure)

Calculates the dew point temperature from the given vapor pressure.

This function uses the formula to convert the vapor pressure to the dew point
temperature in Celsius, and then converts the result to Kelvin.

Arguments:
- vapor_pressure: The vapor pressure in a Unitful pressure unit.

Returns:
- Unitful.Temperature: The dew point temperature in Kelvin.

Example:
```julia
julia> calc.dewpoint(6.11u"hPa")
273.1455029235403 K
```
"""
function dewpoint(vapor_pressure)
    vapor_pressure = uconvert.(u"hPa", vapor_pressure)
    logvalue = log10.(ustrip.(vapor_pressure) ./ 6.112)
    dewPoint = (237.3 .* logvalue) ./ (7.5 .- logvalue) .* u"°C"
    return uconvert.(u"K", dewPoint)
end

"""
dewpoint_from_rh(T_air::Unitful.Temperature, RH::Real)

Calculates the dew point temperature given the air temperature and relative humidity.

This function uses the Magnus-Tetens approximation to estimate the saturation vapor pressure
and iteratively adjusts the dew point temperature until the calculated vapor pressure matches
the saturation vapor pressure corresponding to the dew point.

Arguments:
- T_air::Unitful.Temperature: The current air temperature.
- RH::Real: The relative humidity as a real number between 0 and 1.

Returns:
- Unitful.Temperature: The dewpoint temperature in Kelvin.

Example:
```julia
julia> calc.dewpoint_from_rh(300u"K",0.5)
288.701999999987 K

julia> calc.dewpoint_from_rh([300,300,300]*u"K",[0.5,0.5,0.5])
3-element Vector{Quantity{Float64, 𝚯, Unitful.FreeUnits{(K,), 𝚯, nothing}}}:
 288.701999999987 K
 288.701999999987 K
 288.701999999987 K
```
Note:
- The function assumes that the input temperature is in a Unitful temperature unit and the relative humidity is a real number between 0 and 1.
- The function uses an iterative approach to find the dew point temperature, which may not converge for very low temperatures or high humidities.
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

"""
rh_from_dewpoint(T_air, dewpoint)

Calculates the relative humidity (RH) given the air temperature and dew point temperature.

This function uses the formula for saturation vapor pressure as a function of temperature
to calculate the vapor pressures at the air temperature and dew point temperature, then
computes the relative humidity as the ratio of these two vapor pressures.

Arguments:
- T_air: The current air temperature in a Unitful temperature unit.
- dewpoint: The dew point temperature in a Unitful temperature unit.

Returns:
- A Unitful.Quantity: The relative humidity as a real number between 0 and 1.

Example:
```julia
julia> calc.rh_from_dewpont(300u"K",288u"K")
0.4779147464796312

julia> calc.rh_from_dewpont([300,300,300]u"K",[288,288,288]u"K")
3-element Vector{Float64}:
 0.4779147464796312
 0.4779147464796312
 0.4779147464796312
```
"""
function rh_from_dewpont(T_air, dewpoint)
    T_air = uconvert.(u"°C", T_air)
    dewpoint = uconvert.(u"°C", dewpoint)
    e = vapor_pressure_from_dewpoint.(dewpoint)
    e_s = vapor_pressure_from_dewpoint.(T_air)
    return e ./ e_s
end

"""
dewpoint_from_sh(pressure, q)

Calculates the dew point temperature from the given pressure and specific humidity.

This function uses the saturation vapor pressure formula to calculate the dew point
temperature. It first converts the input pressure to hectopascals (hPa) and the specific
humidity to grams per kilogram (g/kg). Then, it calculates the vapor pressure using the
relationship between specific humidity, vapor pressure, and pressure. Finally, it uses
the `dewpoint` function to find the dew point temperature corresponding to the calculated
vapor pressure.

Arguments:
- pressure::Unitful.Pressure: The atmospheric pressure.
- q::Unitful.Quantity: The specific humidity as a ratio (g/kg).

Returns:
- Unitful.Temperature: The dew point temperature.

Example:
```julia
julia> calc.dewpoint_from_sh(1013.25u"hPa", 0.01u"g/kg")
212.51310407992628 K
```
"""
function dewpoint_from_sh(pressure, q)
    q = uconvert.(u"g/kg", q)
    pressure = uconvert.(u"hPa", pressure)
    w = q ./ (1 .- q)
    e = w .* pressure ./ (0.622 .+ w)
    return dewpoint.(e)
end

"""
vapor_pressure_from_dewpoint(dewPoint)

Calculates the saturation vapor pressure from the dew point temperature.

This function uses the Tetens formula to estimate the saturation vapor pressure
from the dew point temperature. The Tetens formula is an approximation that is
valid for temperatures between 0°C and 50°C.

Arguments:
- dewPoint: The dew point temperature in a Unitful temperature unit.

Returns:
- Unitful.Pressure: The saturation vapor pressure in hectopascals (hPa).

Example:
```julia
julia> calc.vapor_pressure_from_dewpoint(20u"°C")
23.397012933719008 hPa
```
"""
function vapor_pressure_from_dewpoint(dewPoint)
    dewPoint = uconvert.(u"°C", dewPoint)
    logvalue = (7.5 .* ustrip(dewPoint)) ./ (237.3 .+ ustrip(dewPoint))
    vapor_pressure = 6.112 .* (10) .^ logvalue .* u"hPa"
    return vapor_pressure
end

"""
    q_from_dewpoint(pressure, dewpoint)

Calculates the specific humidity (q) given the atmospheric pressure and dew point temperature.

This function uses the dew point temperature to calculate the saturation vapor pressure (e) using the
`vapor_pressure_from_dewpoint` function, and then computes the specific humidity using the formula:

    q = 0.622 * e / (pressure - 0.378 * e)

where:
- q is the specific humidity (in kg/kg),
- e is the saturation vapor pressure (in hPa),
- pressure is the atmospheric pressure (in hPa),
- dewpoint is the dew point temperature (in Kelvin).

Arguments:
- pressure::Unitful.Pressure: The atmospheric pressure.
- dewpoint::Unitful.Temperature: The dew point temperature.

Returns:
- q Unitful.Quantity: The specific humidity (in g/kg).

Example:
```julia
julia> calc.q_from_dewpoint(983u"hPa", 298.15u"K")
20.303761468426334 g kg⁻¹
```
"""
function q_from_dewpoint(pressure, dewpoint)
    pressure = uconvert.(u"hPa", pressure)
    dewpoint = uconvert.(u"K", dewpoint)
    e = vapor_pressure_from_dewpoint.(dewpoint)
    q = (0.622 .* e ./ (pressure .- 0.378 .* e)) .* u"kg/kg"
    return uconvert.(u"g/kg", q)
end
"""
    q_from_vapor_pressure(pressure, vapor_pressure)

Calculates the specific humidity (`q`) from the given atmospheric pressure and vapor pressure.

This function uses the formula to convert vapor pressure to specific humidity, which is the ratio of the mass of water vapor to the total mass of the air. The formula takes into account the molecular weight difference between water vapor and dry air, and the assumption that the air behaves as an ideal gas.

Arguments:
- pressure::Unitful.Pressure: The atmospheric pressure in which the vapor pressure is measured.
- vapor_pressure::Unitful.Pressure: The vapor pressure of water in the air.

Returns:
- q Unitful.Quantity: The specific humidity in grams per kilogram (`g/kg`).

Example:
```julia
julia> calc.q_from_vapor_pressure(1013.25u"hPa", 2.3u"hPa")
1.4131049133102709 g kg⁻¹
"""
function q_from_vapor_pressure(pressure, vapor_pressure)
    pressure = uconvert.(u"hPa", pressure)
    vapor_pressure = uconvert.(u"hPa", vapor_pressure)
    e = vapor_pressure
    q = (0.622 .* e ./ (pressure .- 0.378 .* e)) .* u"kg/kg"
    return uconvert.(u"g/kg", q)
end

"""
    LCL(pressure::Unitful.Pressure, T_air::Unitful.Temperature, dewpoint::Unitful.Temperature)

Calculates the lifting condensation level (LCL) of a parcel of air.

The LCL is the pressure level at which a parcel of air becomes saturated and begins to condense
when lifted dry adiabatically from the surface. This function uses an iterative approach to find
the LCL by comparing the specific humidity of the parcel at a given pressure level to the
specific humidity at the surface.

Arguments:
- pressure: The pressure at the surface in Unitful pressure units.
- T_air: The temperature at the surface in Unitful temperature units.
- dewpoint: The dew point temperature at the surface in Unitful temperature units.

Returns:
- A tuple of two Unitful quantities: the temperature and pressure at the LCL.

Example:
```julia
julia> calc.LCL(1013.25u"hPa", 25u"°C", 20u"°C")
(291.9970210003215 K, 942 hPa)
```
"""
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

"""
    theta(T_air, Pressure)

Calculates the potential temperature (θ) given the air temperature and pressure.

This function uses the formula for potential temperature, which is a measure of the temperature of air after it has been brought adiabatically to a reference pressure, usually 1000 hPa. The formula is:

θ = T_air * (1000 / Pressure)^0.286

Arguments:
- T_air: The current air temperature in a Unitful temperature unit.
- Pressure: The atmospheric pressure in a Unitful pressure unit.

Returns:
- A Unitful.Quantity: The potential temperature in Kelvin.

Example:
```julia
julia> calc.theta(25u"°C", 800u"hPa")
317.79793975452117 K
```
"""
function theta(T_air, Pressure)
    T_air = ustrip.(uconvert.(u"K", T_air))
    Pressure = ustrip.(uconvert.(u"hPa", Pressure))
    return T_air .* (1000 ./ Pressure) .^ (0.286) .* u"K"
end

# height and pressure values
"""
add_height_to_pressure(pressure, height)

Adjusts the pressure level by adding a specified height, assuming standard atmospheric conditions.

This function uses a standard atmospheric model to convert the input pressure to a height,
adds the specified height to this pressure level, and then converts the new height back to
a pressure level using the same standard model.

Arguments:
- pressure: The pressure level in a Unitful pressure unit.
- height: The height to be added in a Unitful length unit.

Returns:
- Unitful.Pressure: The new pressure level after adding the specified height, in hPa.

Example:
```julia
julia> calc.add_height_to_pressure(1000u"hPa",800u"m")
908.303207892483 hPa
```
"""
function add_height_to_pressure(pressure, height)
    pressure = uconvert.(u"hPa", pressure)
    height = uconvert.(u"m", height)
    pressure_level_height = basefun.pressure_to_height_std.(pressure)
    return basefun.height_to_pressure_std.(pressure_level_height .+ height)
end

"""
add_pressure_to_height(height, pressure)

Adjusts a given height by a certain pressure difference to find the new height at which the adjusted pressure would occur.

This function uses standard atmospheric pressure to height conversion functions to determine the new height.
It first converts the input height to meters and pressure to hectopascals (hPa), then finds the pressure level corresponding to the input height.
Finally, it subtracts the given pressure from this pressure level to find the new height at which this adjusted pressure would occur.

Arguments:
- height: The initial height in a Unitful length unit.
- pressure: The pressure difference to be subtracted in a Unitful pressure unit.

Returns:
- A Unitful.Length: The new height at which the adjusted pressure would occur.

Example:
```julia
julia> add_pressure_to_height(1000u"m", 10u"hPa")
1.0911907787454445 km
```
"""
function add_pressure_to_height(height, pressure)
    height = uconvert.(u"m", height)
    pressure = uconvert.(u"hPa", pressure)
    height_level_pressure = basefun.height_to_pressure_std.(height)
    return basefun.pressure_to_height_std.(height_level_pressure .- pressure)
end

# wind_chill_index
"""
wind_chill_index(wind_speed, T_air)

Calculates the wind chill index (WCI) based on the wind speed and air temperature.

This function implements the formula for the wind chill index, which estimates the
perceived temperature felt by the human body due to the cooling effect of the wind.
The formula is:

WCI = (10 * sqrt(wind_speed) - wind_speed + 10.5) * (33 - T_air)

Arguments:
- wind_speed: The wind speed in a Unitful velocity unit.
- T_air: The air temperature in a Unitful temperature unit.

Returns:
- A Unitful.Quantity: The wind chill index with units of kg*cal/m^2/h.

Example:
```julia
julia> calc.wind_chill_index(20u"m/s", 290u"K")
568.824956732432 cal kg m⁻² s⁻¹
```
"""
function wind_chill_index(wind_speed, T_air)
    wind_speed = ustrip.(uconvert.(u"m/s", wind_speed))
    T_air = ustrip.(uconvert.(u"°C", T_air))
    WCI_up = (10 .* sqrt.(wind_speed) .- wind_speed .+ 10.5)
    WCI_down = (33 .- T_air)
    return (WCI_up .* WCI_down) * u"kg*cal/m^2/s"
end

"""
cape_cin(T_air, dewpoint, pressure, T_stra)
"""
function cape_cin(
    T_air::typeof([300, 300, 300, 300] * u"K"),
    dewpoint::typeof([300, 300, 300, 300] * u"K"),
    pressure::typeof([1000, 900, 800, 700] * u"hPa"),
    T_parcel_profile::typeof([300, 300, 300, 300] * u"K"),
)
    # Calculate parcel's profile
end

"""
TODO 计算层结曲线
"""
function parcel_profile(
    pressure::typeof([1000, 900, 800, 700] * u"hPa"),
    T_start::Unitful.Temperature,
    dewpoint_start::Unitful.Temperature,
)

end

"""
downdraft_cape
"""
function dcape()

end

"""
涡度
"""
function vo()

end

"""
散度
"""
function div()

end

"""
垂直风切
"""
function vertical_wind_shear(down_winds_u, down_winds_v, up_winds_u, up_winds_v)
    down_winds_u = ustrip.(uconvert.(u"m/s", down_winds_u))
    down_winds_v = ustrip.(uconvert.(u"m/s", down_winds_v))
    up_winds_u = ustrip.(uconvert.(u"m/s", up_winds_u))
    up_winds_v = ustrip.(uconvert.(u"m/s", up_winds_v))
    return sqrt.((up_winds_u - down_winds_u) .^ 2 .+ (up_winds_v .- down_winds_v) .^ 2) .*
           u"m/s"
end

"""
垂直速度
"""
function vertical_speed()

end


end
