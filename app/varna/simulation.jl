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

using Printf
using FileIO
using JLD2
using ClimaOcean
using Oceananigans

include("setup.jl")

using .FjordsSim: progress, safe_execute

## Model Setup
# sim_setup = setup_varna_3d()
# sim_setup = setup_varna_3d_Lobster()
sim_setup = setup_varna_3d_OXYDEP()

grid = sim_setup.grid_callable!(sim_setup)
buoyancy = sim_setup.buoyancy
closure = sim_setup.closure
tracer_advection = sim_setup.tracer_advection
momentum_advection = sim_setup.momentum_advection
tracers = sim_setup.tracers
free_surface = sim_setup.free_surface_callable(sim_setup.free_surface_args...)
coriolis = sim_setup.coriolis
forcing = safe_execute(sim_setup.forcing_callable)(sim_setup.forcing_args...)
boundary_conditions = safe_execute(sim_setup.bc_callable)(sim_setup.bc_args...)
biogeochemistry = safe_execute(sim_setup.biogeochemistry_callable)(sim_setup.biogeochemistry_args...)

## Model
ocean_model = HydrostaticFreeSurfaceModel(;
    grid,
    buoyancy,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    free_surface,
    coriolis,
    forcing,
    boundary_conditions,
    biogeochemistry
)

## Simulation
Δt = 1seconds
ocean_sim = Simulation(ocean_model; Δt)

## Set initial conditions
set!(ocean_model, T = 10, S = 15)

## Coupled model / simulation
sea_ice = nothing
atmosphere = safe_execute(sim_setup.atmosphere_callable)(sim_setup.atmosphere_args...)
radiation = sim_setup.radiation
coupled_model = OceanSeaIceModel(ocean_sim, sea_ice; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model; Δt)

## Callbacks
coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

## Set up output writers
surface_prefix = joinpath(homedir(), "data_Varna", "surface_snapshots")
ocean_sim.output_writers[:surface] = JLD2OutputWriter(
    ocean_model, merge(ocean_model.tracers, ocean_model.velocities);
    schedule = TimeInterval(1hour),
    filename = "$surface_prefix.jld2",
    indices=(:, :, grid.Nz),
    overwrite_existing = true,
    array_type=Array{Float32}
)

profile_prefix = joinpath(homedir(), "data_Varna", "profile_snapshots")
ocean_sim.output_writers[:profile] = JLD2OutputWriter(
    ocean_model, merge(ocean_model.tracers, ocean_model.velocities);
    schedule = TimeInterval(1hour),
    filename = "$profile_prefix.jld2",
    indices=(:, 18, :),
    overwrite_existing = true,
    array_type=Array{Float32}
)

## Spinning up the simulation
# We use an adaptive time step that maintains the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) equal to 0.1.
ocean_sim.stop_time = 10days
coupled_simulation.stop_time = 10days
conjure_time_step_wizard!(ocean_sim; cfl=0.1, max_Δt=1.5minutes, max_change=1.01)
run!(coupled_simulation)
nothing #hide

## Running the simulation
# This time, we set the CFL in the time_step_wizard to be 0.25 as this is the maximum recommended CFL to be
# used in conjunction with Oceananigans' hydrostatic time-stepping algorithm ([two step Adams-Bashfort](https://en.wikipedia.org/wiki/Linear_multistep_method))
ocean_sim.stop_time = 355days
coupled_simulation.stop_time = 355days
conjure_time_step_wizard!(ocean_sim; cfl=0.2, max_Δt=1.5minutes, max_change=1.01)
run!(coupled_simulation)
nothing #hide
