# Copyright 2024 The FjordsSim Authors.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

## Packages and modules
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

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers,
    update_tendencies!

include("../../src/FjordsSim.jl")
include("setup.jl")

using .FjordsSim
using .FjordsSim: ImmersedBoundaryGrid

## Setup
# save_interval = 30minutes
setup_grid = SetupGridPredefinedFromFile(args_grid...)
grid = ImmersedBoundaryGrid(setup_grid)

## Biogeochemistry
const year = years = 365days
# Surface PAR
# Setting up idealised functions for PAR and diffusivity 
# (details here can be ignored but these are typical of the North Atlantic)
@inline PAR⁰(x, y, t) =
    60 *
    (1 - cos((t + 15days) * 2π / year)) *
    (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

## Oxydep
biogeochemistry = OXYDEP(; grid, surface_photosynthetically_active_radiation = PAR⁰)

## Initial Conditions
# z_ini = -reverse([0.5, 1, 2, 3, 4, 6, 8, 9, 10, 12, 16, 20])
# tprof = reverse([20, 20, 20, 20, 18, 15, 14, 13, 12, 12, 12, 12])
# itp = LinearInterpolation(z_ini, tprof)
# tprof_target = itp(z_middle)
# 
# T₀ = Array{Float64}(undef, Nx, Ny, Nz)
# for i = 1:Nx
#     for j = 1:Ny
#         T₀[i, j, :] = tprof_target
#     end
# end
# 
# sprof = reverse([14, 14.1, 14.2, 14.5, 14.8, 15, 15.1, 15.2, 15.3, 15.4, 15.5, 15.5])
# itps = LinearInterpolation(z_ini, sprof)
# sprof_target = itps(z_middle)
# 
# S₀ = Array{Float64}(undef, Nx, Ny, Nz)
# for i = 1:Nx
#     for j = 1:Ny
#         S₀[i, j, :] = sprof_target
#     end
# end

## Physics
# const surface_νz = 1e-2
# const background_νz = 1e-4
# const background_κz = 1e-5
# 
# @inline νz(x, y, z, t) = ifelse(z > -15, surface_νz, background_νz)
# 
# horizontal_viscosity = HorizontalScalarDiffusivity(ν = 1e3)
# vertical_viscosity =
#     VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = background_κz)
# convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 0.1)
# κ_skew = 900.0      # [m² s⁻¹] skew diffusivity
# κ_symmetric = 900.0 # [m² s⁻¹] symmetric diffusivity
# 
# gent_mcwilliams_diffusivity =
#     IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric, slope_limiter = FluxTapering(1e-2))
# 
# closure =
#     (vertical_viscosity, horizontal_viscosity, convective_adjustment, gent_mcwilliams_diffusivity)

## Boundary confitions
reference_density = 1000.0
reference_heat_capacity = 3991.0
reference_salinity = 15

## Wind stress 
# https://en.wikipedia.org/wiki/Wind_stress
Cd = 0.0025  # 0.0015
ρₐᵢᵣ = 1.225
Ntimes = 12
uwind = [10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
vwind = [8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
τˣ = Array{Float64}(undef, Nx, Ny, Ntimes)
τʸ = Array{Float64}(undef, Nx, Ny, Ntimes)
for i = 1:Nx
    for j = 1:Ny
        τˣ[i, j, :] = -ρₐᵢᵣ * Cd .* (uwind .^ 2) ./ reference_density
        τʸ[i, j, :] = ρₐᵢᵣ * Cd .* (vwind .^ 2) ./ reference_density
    end
end

## Time dependent fluxes
const Nyears = 1.0
const Nmonths = 12
const thirty_days = 30days

@inline current_time_index(time, tot_months) =
    mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
@inline next_time_index(time, tot_months) =
    mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
@inline cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / thirty_days, 1) * (u₂ - u₁)

Δz_top = @allowscalar Δzᵃᵃᶜ(1, 1, grid.Nz, grid.underlying_grid)
Δz_bottom = @allowscalar Δzᵃᵃᶜ(1, 1, 1, grid.underlying_grid)

@inline function surface_wind_stress(i, j, grid, clock, fields, τ)
    time = clock.time
    n₁ = current_time_index(time, Nmonths)
    n₂ = next_time_index(time, Nmonths)

    @inbounds begin  # inbounds for faster performance, doesn't check the boundaries of the array
        τ₁ = τ[i, j, n₁]
        τ₂ = τ[i, j, n₂]
    end

    return cyclic_interpolate(τ₁, τ₂, time)
end

## Bottom drag
Δz_top = @allowscalar grid.Δzᵃᵃᶜ[Nz]

# Linear bottom drag:
μ = 0.003 # Non dimensional

@inline speedᶠᶜᶜ(i, j, k, grid, fields) =
    @inbounds sqrt(fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)^2)
@inline speedᶜᶠᶜ(i, j, k, grid, fields) =
    @inbounds sqrt(fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)^2)

@inline u_bottom_drag(i, j, grid, clock, fields, μ) =
    @inbounds -μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) =
    @inbounds -μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) =
    @inbounds -μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields)
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) =
    @inbounds -μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields)

drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form = true, parameters = μ)
drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form = true, parameters = μ)

no_slip_bc = ValueBoundaryCondition(0)

u_immersed_bc = ImmersedBoundaryCondition(
    bottom = drag_u,
    west = no_slip_bc,
    east = no_slip_bc,
    south = no_slip_bc,
    north = no_slip_bc,
)

v_immersed_bc = ImmersedBoundaryCondition(
    bottom = drag_v,
    west = no_slip_bc,
    east = no_slip_bc,
    south = no_slip_bc,
    north = no_slip_bc,
)

## Constructing BC
u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ)
v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ)

u_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τˣ);
v_wind_stress_bc = FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τʸ);

u_bcs = FieldBoundaryConditions(
    top = u_wind_stress_bc,
    bottom = u_bottom_drag_bc,
    immersed = u_immersed_bc,
)

v_bcs = FieldBoundaryConditions(
    top = v_wind_stress_bc,
    bottom = v_bottom_drag_bc,
    immersed = v_immersed_bc,
)

## BGC boundary conditions
# for the bottom boundary
O2_suboxic = 30.0  # OXY threshold for oxic/suboxic switch (mmol/m3)
Trel = 10000.0      # Relaxation time for exchange with the sediments (s/m)
b_ox = 15.0        # difference of OXY in the sediment and water, 
b_NUT = 15.0        # NUT in the sediment, (mmol/m3)  
b_DOM_ox = 10.0    # OM in the sediment (oxic conditions), (mmol/m3) 
b_DOM_anox = 20.0   # OM in the sediment (anoxic conditions), (mmol/m3)  
bu = 0.7           # Burial coeficient for lower boundary (0<Bu<1), 1 - for no burying, (nd)

@inline F_ox(conc, threshold) = (0.5 + 0.5 * tanh(conc - threshold))
@inline F_subox(conc, threshold) = (0.5 - 0.5 * tanh(conc - threshold))

#---OXY----------------------
OXY_top = GasExchange(; gas = :O₂)
@inline OXY_bottom_cond(i, j, grid, clock, fields) = @inbounds -(
    F_ox(fields.O₂[i, j, 1], O2_suboxic) * b_ox +
    F_subox(fields.O₂[i, j, 1], O2_suboxic) * (0.0 - fields.O₂[i, j, 1])
) / Trel
OXY_bottom = FluxBoundaryCondition(OXY_bottom_cond, discrete_form = true)

#---NUT----------------------
@inline NUT_bottom_cond(i, j, grid, clock, fields) = @inbounds (
    F_ox(fields.O₂[i, j, 1], O2_suboxic) * (b_NUT - fields.NUT[i, j, 1]) +
    F_subox(fields.O₂[i, j, 1], O2_suboxic) * (0.0 - fields.NUT[i, j, 1])
) / Trel
NUT_bottom = FluxBoundaryCondition(NUT_bottom_cond, discrete_form = true) #ValueBoundaryCondition(10.0)

#---PHY----------------------
# w_PHY = biogeochemical_drift_velocity(biogeochemistry, Val(:PHY)).w[1, 1, 1]
# @inline PHY_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_PHY * fields.PHY[i, j, 1]
# PHY_bottom = FluxBoundaryCondition(PHY_bottom_cond, discrete_form = true)

#---HET----------------------
# w_HET = biogeochemical_drift_velocity(biogeochemistry, Val(:HET)).w[1, 1, 1]
# @inline HET_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_HET * fields.HET[i, j, 1]
# HET_bottom = FluxBoundaryCondition(HET_bottom_cond, discrete_form = true)

#---POM----------------------
# w_POM = biogeochemical_drift_velocity(biogeochemistry, Val(:POM)).w[1, 1, 1]
# @inline POM_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_POM * fields.POM[i, j, 1]
# POM_bottom = FluxBoundaryCondition(POM_bottom_cond, discrete_form = true)

#---DOM----------------------
DOM_top = ValueBoundaryCondition(0.0)
@inline DOM_bottom_cond(i, j, grid, clock, fields) = @inbounds (
    F_ox(fields.O₂[i, j, 1], O2_suboxic) * (b_DOM_ox - fields.DOM[i, j, 1]) +
    F_subox(fields.O₂[i, j, 1], O2_suboxic) * 2.0 * (b_DOM_anox - fields.DOM[i, j, 1])
) / Trel
DOM_bottom = FluxBoundaryCondition(DOM_bottom_cond, discrete_form = true) #, parameters = (; O2_suboxic, b_DOM_ox, Trel),)

oxy_bcs = FieldBoundaryConditions(top = OXY_top, bottom = OXY_bottom)
nut_bcs = FieldBoundaryConditions(bottom = NUT_bottom)
dom_bcs = FieldBoundaryConditions(top = DOM_top, bottom = DOM_bottom)
# pom_bcs = FieldBoundaryConditions(bottom = POM_bottom)
# phy_bcs = FieldBoundaryConditions(bottom = PHY_bottom)
# het_bcs = FieldBoundaryConditions(bottom = HET_bottom)

boundary_conditions = (
    u = u_bcs,
    v = v_bcs,
    O₂ = oxy_bcs,
    NUT = nut_bcs,
    DOM = dom_bcs,
    # POM = pom_bcs,
    # PHY = phy_bcs,
    # HET = het_bcs,
)

## River forcing
λ = 1 / (1minute)  # Relaxation timescale [s⁻¹].

# Temperature and salinity of the meltwater outflow.
T_source = 15
S_source = 0

# Index of the point source at the middle of the southern wall.
source_index = (1, 13, Nz)

# # Point source
# @inline T_point_source(i, j, k, grid, time, U, C, p) =
#     @inbounds ifelse((i, j, k) == p.source_index, -p.λ * (C.T[i, j, k] - p.T_source), 0)

# @inline S_point_source(i, j, k, grid, time, U, C, p) =
#     @inbounds ifelse((i, j, k) == p.source_index, -p.λ * (C.S[i, j, k] - p.S_source), 0)

# Point source
T_point_source(i, j, k, grid, clock, model_fields) =
    @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.T[i, j, k] - T_source), 0)

S_point_source(i, j, k, grid, clock, model_fields) =
    @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.S[i, j, k] - S_source), 0)

# params = (source_index=source_index, T_source=T_source, S_source=S_source, λ=λ)
params = (T_source = T_source, S_source = S_source, λ = λ)

# Tforcing = Forcing(T_point_source, parameters=params)
Tforcing = Forcing(T_point_source, field_dependencies = :T, discrete_form = true)
Sforcing = Forcing(S_point_source, field_dependencies = :S, discrete_form = true)

## Model
model = HydrostaticFreeSurfaceModel(;
    grid,
    closure,
    biogeochemistry,
    buoyancy = SeawaterBuoyancy(),
    # boundary_conditions,
    # forcing = (T = Tforcing, S = Sforcing),
    momentum_advection = VectorInvariant(),
    tracer_advection = WENO(underlying_grid),
    tracers = (:NUT, :PHY, :HET, :POM, :DOM, :O₂, :T, :S),
)

## Set initial conditions
set!(model, T = T₀, S = S₀, NUT = 10.0, PHY = 0.01, HET = 0.05, O₂ = 350.0, DOM = 1.0)

## Simulation
Δt = 1seconds
stop_time = 24hours
simulation = Simulation(model; Δt, stop_time)
progress(sim) = @info "Time : $(prettytime(sim.model.clock.time)),
                        max(|u|): $(maximum(abs, sim.model.velocities.u)),
                        max(S): $(maximum(model.tracers.S)),
                        Δt: $(prettytime(sim.Δt))"
simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))
u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
O₂ = model.tracers.O₂
NUT = model.tracers.NUT
DOM = model.tracers.DOM
POM = model.tracers.POM
PHY = model.tracers.PHY
HET = model.tracers.HET

output_prefix = joinpath(homedir(), "data_Varna", "simulation_snapshots")
simulation.output_writers[:surface_fields] = JLD2OutputWriter(
    model,
    (; u, v, w, T, S, O₂, NUT, DOM, POM, PHY, HET),
    schedule = TimeInterval(2minutes),
    filename = "$output_prefix.jld2",
    with_halos = true,
    overwrite_existing = true,
)

# checkpoints_prefix = joinpath(homedir(), "data_Varna", "simulation_chkpnt")
# simulation.output_writers[:checkpointer] =
    # Checkpointer(model; schedule = IterationInterval(200), prefix = checkpoints_prefix)
    # Checkpointer(model; schedule = IterationInterval(200), prefix = "model_checkpoint")


## run simulation
run!(simulation)  # , pickup = true)
