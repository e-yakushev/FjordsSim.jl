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

# using CSV
# using CUDA: @allowscalar
# using CairoMakie
# using DataFrames
# using DelimitedFiles
# using FileIO
# using JLD2: JLD2OutputWriter
# using NCDatasets
using Oceananigans:
    SeawaterBuoyancy,
    VectorInvariant,
    WENO,
    set!,
    Simulation,
    TimeInterval,
    prettytime,
    Callback,
    IterationInterval,
    JLD2OutputWriter,
    run!
using Oceananigans.Models: HydrostaticFreeSurfaceModel
# using Oceananigans.Forcings
using Oceananigans.Units
# using Printf
# using Statistics

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers,
    update_tendencies!

include("../../src/FjordsSim.jl")
include("setup.jl")

using .FjordsSim:
    ImmersedBoundaryGrid,
    OXYDEP,
    SetupGridPredefinedFromFile,
    initial_conditions_temp_salt_3d_predefined,
    turbulence_closures_a,
    # bgh_oxydep_boundary_conditions,
    # rivers_forcing,
    physics_boundary_conditions

## Grid
setup_grid = SetupGridPredefinedFromFile(args_grid...)
grid = ImmersedBoundaryGrid(setup_grid)

## Biogeochemistry
# (details here can be ignored but these are typical of the North Atlantic)
const year = years = 365days
@inline PAR⁰(x, y, t) =
    60 *
    (1 - cos((t + 15days) * 2π / year)) *
    (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
biogeochemistry = OXYDEP(; args_oxydep..., grid, surface_photosynthetically_active_radiation = PAR⁰)

## The turbulence closure
closure = turbulence_closures_a()

## Boundary conditions
# physics
u_bcs, v_bcs = physics_boundary_conditions()

# BGC boundary conditions
# todo

boundary_conditions = (u = u_bcs, v = v_bcs)

## River forcing
# forcing = rivers_forcing()

## Model
model = HydrostaticFreeSurfaceModel(;
    grid,
    closure,
    biogeochemistry,
    buoyancy = SeawaterBuoyancy(),
    boundary_conditions,
    # forcing = forcing,
    momentum_advection = VectorInvariant(),
    tracer_advection = WENO(grid.underlying_grid),
    tracers = (:NUT, :PHY, :HET, :POM, :DOM, :O₂, :T, :S),
)

## Set initial conditions
T₀, S₀ = initial_conditions_temp_salt_3d_predefined(setup_grid)
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
