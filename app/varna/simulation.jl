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

using FileIO
using JLD2
using ClimaOcean
using Oceananigans
using Oceananigans.OutputReaders: InMemory
# using Oceananigans.Models: HydrostaticFreeSurfaceModel
# using Oceananigans.Units

# import Oceananigans.Biogeochemistry:
#     biogeochemical_drift_velocity,
#     required_biogeochemical_auxiliary_fields,
#     required_biogeochemical_tracers,
#     update_tendencies!

include("../../src/FjordsSim.jl")
include("setup.jl")

using .FjordsSim:
    # ImmersedBoundaryGrid,
    # OXYDEP,
    # SetupGridPredefinedFromFile,
    initial_conditions_temp_salt_3d_predefined,
    turbulence_closures_a,
    # bgh_oxydep_boundary_conditions,
    rivers_forcing,
    physics_boundary_conditions

## Grid
filepath_topo = joinpath(args_grid.datadir, args_grid.filename)
@load filepath_topo depth
Nx, Ny = size(depth)
Nz = size(args_grid.z_middle)[1]
setup_grid = (; args_grid..., Nx = Nx, Ny = Ny, Nz = Nz)
z_faces = exponential_z_faces(; Nz=Nz, depth=20)
underlying_grid = LatitudeLongitudeGrid(setup_grid.arch;
     size=(setup_grid.Nx, setup_grid.Ny, setup_grid.Nz),
     halo=(7, 7, 7),
     z=z_faces,
     latitude=(43.177, 43.214),
     longitude=(27.640, 27.947))
# depth[-setup_grid.minimum_depth .<= depth .< 0] .= -setup_grid.minimum_depth
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)

## Biogeochemistry
# (details here can be ignored but these are typical of the North Atlantic)
# const year = years = 365days
# @inline PAR⁰(x, y, t) =
#     60 *
#     (1 - cos((t + 15days) * 2π / year)) *
#     (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
# biogeochemistry = OXYDEP(; grid, args_oxydep..., TS_forced = false, surface_photosynthetically_active_radiation = PAR⁰)

backend = InMemory()
atmosphere = JRA55_prescribed_atmosphere(setup_grid.arch; backend, grid)
# radiation  = Radiation(arch)

## The turbulence closure
closure = turbulence_closures_a()

## Boundary conditions
# physics
u_bcs, v_bcs = physics_boundary_conditions(setup_grid.arch, setup_grid.Nx, setup_grid.Ny)

# BGC boundary conditions
# todo

boundary_conditions = (u = u_bcs, v = v_bcs)

## River forcing
forcing = rivers_forcing(setup_grid.Nz)

## Model
model = HydrostaticFreeSurfaceModel(;
    grid,
    closure,
    # biogeochemistry,
    buoyancy = SeawaterBuoyancy(),
    boundary_conditions,
    forcing = forcing,
    momentum_advection = VectorInvariant(),
    tracer_advection = WENO(grid.underlying_grid),
    # tracers = (:NUT, :PHY, :HET, :POM, :DOM, :O₂, :T, :S),
)

## Set initial conditions
T₀, S₀ = initial_conditions_temp_salt_3d_predefined(setup_grid)
set!(model, T = T₀, S = S₀)  # , NUT = 10.0, PHY = 0.01, HET = 0.05, O₂ = 350.0, DOM = 1.0)

## Simulation
Δt = 0.5seconds
stop_time = 90days
simulation = Simulation(model; Δt, stop_time)
# conjure_time_step_wizard!(simulation; cfl = 0.1, max_Δt = 1, max_change = 1.01)

progress(sim) = @info "Time : $(prettytime(sim.model.clock.time)),
    max(|u|): $(maximum(abs, sim.model.velocities.u)),
    max(S): $(maximum(model.tracers.S)),
    Δt: $(prettytime(sim.Δt)),
    walltime: $(prettytime(sim.run_wall_time))"
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
# O₂ = model.tracers.O₂
# NUT = model.tracers.NUT
# DOM = model.tracers.DOM
# POM = model.tracers.POM
# PHY = model.tracers.PHY
# HET = model.tracers.HET

output_prefix = joinpath(homedir(), "data_Varna", "simulation_snapshots")
simulation.output_writers[:surface_fields] = JLD2OutputWriter(
    model,
    (; u, v, w, T, S),
    schedule = TimeInterval(1hour),
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
