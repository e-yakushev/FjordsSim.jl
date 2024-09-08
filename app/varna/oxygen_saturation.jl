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

# Function to calculate oxygen saturation in seawater
function oxygen_saturation(T::Float64, S::Float64, P::Float64 = 0.0)::Float64
    T_kelvin = T + 273.15  # Convert temperature to Kelvin

    # Calculate the natural logarithm of oxygen saturation concentration
    ln_O2_sat = A1 + A2 * (100 / T_kelvin) + A3 * log(T_kelvin / 100) +
                A4 * T_kelvin / 100 + A5 * (T_kelvin / 100)^2 + A6 * (T_kelvin / 100)^3 +
                S * (B1 + B2 * (T_kelvin / 100) + B3 * (T_kelvin / 100)^2) + C1 * S^2

    # Oxygen saturation concentration in µmol/kg
    O2_sat = exp(ln_O2_sat) * 44.66

    # Pressure correction factor (Weiss, 1970)
    P_corr = 1.0 + P * (5.6e-6 + 2.0e-11 * P)

    # Adjusted oxygen saturation with pressure correction
    O2_sat_corrected = O2_sat * P_corr

    return O2_sat_corrected
end

# Example usage:
#T = 15.0  # Temperature in °C
#S = 35.0  # Salinity in PSU
#P = 0.0   # Pressure in dbar (surface pressure is 0)

#O2_sat = oxygen_saturation(T, S, P)
#println("Oxygen saturation concentration: $O2_sat µmol/kg")
