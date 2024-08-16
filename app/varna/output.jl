using Oceananigans
using Oceananigans.Units: hour, days, meters
using JLD2
using CairoMakie

using OceanBioME, Oceananigans, Printf
using OceanBioME: Boundaries, GasExchange
using OceanBioME.Boundaries.Sediments: sinking_flux
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

include("../../src/FjordsSim.jl")
using .FjordsSim:
       OXYDEP

year=365days
stoptime = 730
#depth_extent = 10meters
#grid = RectilinearGrid(size = (1, 1, 10), extent = (500meters, 500meters, 10meters), topology = (Bounded, Bounded, Bounded))
grid = RectilinearGrid(size = (1, 1, 12), extent = (500meters, 500meters, 67meters), topology = (Bounded, Bounded, Bounded))

include("images.jl")