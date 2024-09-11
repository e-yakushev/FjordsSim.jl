"""
`sea_water_flux` to calculate flux of gases at the sea water interface.
now is mad efor oxygen, CO2 and CH4 should be added. Dependence on wind should be replaced 
by that with assuming the bubbles effect. 
"""
#- - - - - - - - - - - - - - - - - - - - - - 
# Function to calculate oxygen saturation in seawater
#- - - - - - - - - - - - - - - - - - - - - - 
# Coefficients from Garcia and Gordon (1992)
const A1 = -173.4292
const A2 = 249.6339
const A3 = 143.3483
const A4 = -21.8492
const A5 = -0.033096
const A6 = 0.014259
const B1 = -0.035274
const B2 = 0.001429
const B3 = -0.00007292
const C1 = 0.0000826

function oxygen_saturation(T::Float64, S::Float64, P::Float64)::Float64

    T_kelvin = T + 273.15  # Convert temperature to Kelvin

    # Calculate the natural logarithm of oxygen saturation concentration
    ln_O2_sat =
        A1 +
        A2 * (100 / T_kelvin) +
        A3 * log(T_kelvin / 100) +
        A4 * T_kelvin / 100 +
        A5 * (T_kelvin / 100)^2 +
        A6 * (T_kelvin / 100)^3 +
        S * (B1 + B2 * (T_kelvin / 100) + B3 * (T_kelvin / 100)^2) +
        C1 * S^2

    # Oxygen saturation concentration in µmol/kg
    O2_sat = exp(ln_O2_sat) * 44.66

    # Pressure correction factor (Weiss, 1970)
    P_corr = 1.0 + P * (5.6e-6 + 2.0e-11 * P)

    # Adjusted oxygen saturation with pressure correction
    return (O2_sat * P_corr)
end

#- - - - - - - - - - - - - - - - - - - - - - 
# Sc, Schmidt number for O2  following Wanninkhof 2014
#- - - - - - - - - - - - - - - - - - - - - - 
@inline function OxygenSchmidtNumber(T::Float64)::Float64
    return ((1920.4 - 135.6 * T + 5.2122 * T^2 - 0.10939 * T^3 + 0.00093777 * T^4))
    # can be replaced by PolynomialParameterisation{4}((a, b, c, d, e)) i.e.:
    #    a = 1953.4, b = - 128.0, c = 3.9918, d = -0.050091, e = 0.00093777  
    # Sc = PolynomialParameterisation{4}((a, b, c, d, e))
end

#- - - - - - - - - - - - - - - - - - - - - - 
# WindDependence, [mmol m-2s-1], Oxygen Sea Water Flux 
#- - - - - - - - - - - - - - - - - - - - - - 
function WindDependence(windspeed::Float64)::Float64
    return (0.251 * windspeed^2.0) #ko2o=0.251*windspeed^2*(Sc/660)^(-0.5)  Wanninkhof 2014
end

#- - - - - - - - - - - - - - - - - - - - - - 
# OxygenSeaWaterFlux, [mmol m-2s-1], Oxygen Sea Water Flux 
#- - - - - - - - - - - - - - - - - - - - - - 
function OxygenSeaWaterFlux(
    T::Float64,
    S::Float64,
    P::Float64,
    O₂::Float64,
    windspeed::Float64,
)::Float64
    return (
        WindDependence(windspeed) *
        (OxygenSchmidtNumber(T) / 660.0)^(-0.5) *
        (O₂ - oxygen_saturation(T, S, P)) *
        0.24 / 86400.0        # 0.24 is to convert from [cm/h] to [m/day]  * 0.24  / 86400.0
    )
end
