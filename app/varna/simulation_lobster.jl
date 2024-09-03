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
using Oceananigans.OutputReaders: InMemory

include("../../src/FjordsSim.jl")
include("setup.jl")

using .FjordsSim:
    rivers_forcing_with_NO₃,
    biogeochemical_simulation

## Grid
filepath_topo = joinpath(args_grid.datadir, args_grid.filename)
@load filepath_topo depth
Nx, Ny = size(depth)
Nz = args_grid.Nz
setup_grid = (; args_grid..., Nx = Nx, Ny = Ny, Nz = Nz)
z_faces = exponential_z_faces(; Nz = Nz, depth = 20, h = 20)
underlying_grid = LatitudeLongitudeGrid(
    setup_grid.arch;
    size = (setup_grid.Nx, setup_grid.Ny, setup_grid.Nz),
    halo = (7, 7, 7),
    z = z_faces,
    latitude = (43.177, 43.214),
    longitude = (27.640, 27.947),
)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)

## River forcing
forcing = rivers_forcing_with_NO₃(setup_grid.Nz)

## Simulation
Δt = 1seconds
ocean_sim = biogeochemical_simulation(grid; Δt, forcing, coriolis = nothing)
model = ocean_sim.model

## Set initial conditions
Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

# set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2409.0, S = 15, T = 10)
set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, S = 15, T = 10)

## Atmosphere
backend = InMemory()
atmosphere = JRA55_prescribed_atmosphere(setup_grid.arch; backend, grid)
radiation = Radiation(setup_grid.arch)
sea_ice = nothing

## Coupled model
coupled_model = OceanSeaIceModel(ocean_sim, sea_ice; atmosphere, radiation)
coupled_simulation = Simulation(coupled_model; Δt, stop_time = 10days)

## Callbacks
wall_time = [time_ns()]
function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf(
        "Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T): %.2f, min(T): %.2f, wtime: %s \n",
        prettytime(ocean.model.clock.time),
        ocean.model.clock.iteration,
        prettytime(ocean.Δt),
        umax...,
        Tmax,
        Tmin,
        prettytime(step_time)
    )

    wall_time[1] = time_ns()
end
coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))
nothing #hide

## Set up output writers
surface_prefix = joinpath(homedir(), "data_Varna", "surface_snapshots")
ocean_sim.output_writers[:surface] = JLD2OutputWriter(
    model,
    merge(model.tracers, model.velocities);
    schedule = TimeInterval(1hour),
    filename = "$surface_prefix.jld2",
    indices = (:, :, grid.Nz),
    overwrite_existing = true,
    array_type = Array{Float32},
)

profile_prefix = joinpath(homedir(), "data_Varna", "profile_snapshots")
ocean_sim.output_writers[:profile] = JLD2OutputWriter(
    model,
    merge(model.tracers, model.velocities);
    schedule = TimeInterval(1hour),
    filename = "$profile_prefix.jld2",
    indices = (:, 18, :),
    overwrite_existing = true,
    array_type = Array{Float32},
)

## Spinning up the simulation
#
# As an initial condition, we have interpolated ECCO tracer fields onto our custom grid.
# The bathymetry of the original ECCO data may differ from our grid, so the initialization of the velocity
# field might cause shocks if a large time step is used.
#
# Therefore, we spin up the simulation with a small time step to ensure that the interpolated initial
# conditions adapt to the model numerics and parameterization without causing instability. A 10-day
# integration with a maximum time step of 1.5 minutes should be sufficient to dissipate spurious
# initialization shocks.
# We use an adaptive time step that maintains the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) equal to 0.1.
# For this scope, we use the Oceananigans utility `conjure_time_step_wizard!` (see Oceanigans's documentation).

ocean_sim.stop_time = 10days
coupled_simulation.stop_time = 10days
conjure_time_step_wizard!(ocean_sim; cfl = 0.1, max_Δt = 40seconds, max_change = 1.01)
run!(coupled_simulation)
nothing #hide

## Running the simulation
#
# Now that the simulation has spun up, we can run it for the full 100 days.
# We increase the maximum time step size to 10 minutes and let the simulation run for 100 days.
# This time, we set the CFL in the time_step_wizard to be 0.25 as this is the maximum recommended CFL to be
# used in conjunction with Oceananigans' hydrostatic time-stepping algorithm ([two step Adams-Bashfort](https://en.wikipedia.org/wiki/Linear_multistep_method))

ocean_sim.stop_time = 355days
coupled_simulation.stop_time = 355days
conjure_time_step_wizard!(ocean_sim; cfl = 0.2, max_Δt = 40seconds, max_change = 1.01)
run!(coupled_simulation)
nothing #hide
