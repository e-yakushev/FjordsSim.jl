using Oceananigans
using ClimaOcean
using ClimaOcean.DataWrangling.JRA55: compute_bounding_nodes, infer_longitudinal_topology
using ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    default_gravitational_acceleration,
    TemperatureDependentAirViscosity,
    MomentumRoughnessLength,
    ScalarRoughnessLength,
    SimilarityScales

import ClimaOcean.DataWrangling.JRA55: compute_bounding_indices

""" Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth """
@inline PAR⁰(x, y, t) =
    60 * (1 - cos((t + 15days) * 2π / 365days)) * (1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days)^2))) + 2

# ClimaOcean v0.5.4 fix for the custom longitude and latitude
# this is called from set! and uses grid to find the locations,
# which are 1 index more than necessary
function compute_bounding_indices(longitude::Nothing, latitude::Nothing, grid, LX, LY, λc, φc)
    λbounds = compute_bounding_nodes(longitude, grid, LX, λnodes)
    φbounds = compute_bounding_nodes(latitude, grid, LY, φnodes)

    i₁, i₂ = compute_bounding_indices(λbounds, λc)
    j₁, j₂ = compute_bounding_indices(φbounds, φc)
    TX = infer_longitudinal_topology(λbounds)

    # to prevent taking larger than grid areas
    i₁ = (i₂ - i₁ >= grid.Nx) ? (i₂ - grid.Nx + 1) : i₁
    j₁ = (j₂ - j₁ >= grid.Ny) ? (j₂ - grid.Ny + 1) : j₁

    return i₁, i₂, j₁, j₂, TX
end

function regional_roughness_lengths(FT = Oceananigans.defaults.FloatType)
    momentum = MomentumRoughnessLength(
        FT;
        gravitational_acceleration = default_gravitational_acceleration,
        maximum_roughness_length = 1.0, # An estimate?
        air_kinematic_viscosity = TemperatureDependentAirViscosity(FT),
        gravity_wave_parameter = 0.011,
        laminar_parameter = 0.11,
    )
    temperature = ScalarRoughnessLength(FT)
    water_vapor = ScalarRoughnessLength(FT)
    return SimilarityScales(momentum, temperature, water_vapor)
end
