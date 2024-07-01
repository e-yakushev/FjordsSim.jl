## Packages and modules
import Oceananigans.Biogeochemistry: biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers,
    update_tendencies!

using CSV
using CUDA: @allowscalar
using CairoMakie
using DataFrames
using DelimitedFiles
using FileIO
using Interpolations
using JLD2
using NCDatasets
using OceanBioME
using OceanBioME.Boundaries.Sediments: sinking_flux
using OceanBioME.SLatissimaModel: SLatissima
using OceanBioME: Boundaries, GasExchange
using Oceananigans
using Oceananigans.Fields: ConstantField, FunctionField
using Oceananigans.Forcings
using Oceananigans.Operators: Δzᵃᵃᶜ, ℑxyᶜᶠᵃ, ℑxyᶠᶜᵃ
using Oceananigans.Units
using Oceananigans: architecture
using Printf
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Statistics

## 
include("../../src/Oxydep.jl")
using .OXYDEPModel

## Setup
arch = CPU()
save_interval = 30minutes

## Grid
z_levels = -reverse([0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20])
z_middle = -reverse([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 
                     14.5, 15.5, 16.5, 17.5, 18.5, 19.5])
Nx = 119
Ny = 42
Nz = length(z_levels)-1
dx = 200  # m
dy = 50  # m
underlying_grid = RectilinearGrid(topology=(Bounded, Bounded, Bounded), size=(Nx, Ny, Nz),
                                        x = (0, dx*Nx), y = (0, dy*Ny),
                                        z = z_levels, halo=(7,7,7))

filepath_topo = joinpath(homedir(), "data_Varna", "Varna_topo.jld2")
@load filepath_topo depth
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)

## Biogeochemistry
const year = years = 365days
# Surface PAR
# Setting up idealised functions for PAR and diffusivity 
# (details here can be ignored but these are typical of the North Atlantic)
@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) * 
                        (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

## Oxydep
biogeochemistry = OXYDEP(; grid, 
                          surface_photosynthetically_active_radiation = PAR⁰,
                          particles = nothing)
