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

using Oceananigans.Units: second, seconds, minute, minutes, hour, hours, day, days
using Oceananigans.Utils: TimeInterval, IterationInterval
using Oceananigans.Simulations: Callback, conjure_time_step_wizard!, run!
using Oceananigans.OutputWriters: JLD2OutputWriter, NetCDFOutputWriter
using Oceanostics
using FjordsSim: coupled_hydrostatic_simulation

include("setup.jl")

## Model Setup
sim_setup = setup_region_3d()

coupled_simulation = coupled_hydrostatic_simulation(sim_setup)

## Set up output writers
ocean_sim = coupled_simulation.model.ocean
ocean_sim.callbacks[:progress] = Callback(ProgressMessengers.TimedMessenger(), IterationInterval(1000));
ocean_model = ocean_sim.model

prefix = joinpath(sim_setup.results_dir, "snapshots_ocean")
ocean_sim.output_writers[:ocean] = JLD2OutputWriter(
    ocean_model,
    merge(ocean_model.tracers, ocean_model.velocities, coupled_simulation.model.interfaces.atmosphere_ocean_interface.fluxes);
    schedule = TimeInterval(1hours),
    filename = "$prefix.jld2",
    overwrite_existing = true,
    array_type = Array{Float32},
)

atmosphere_fields = coupled_simulation.model.interfaces.exchanger.exchange_atmosphere_state
atmosphere_data = NamedTuple((
    u_atm  = atmosphere_fields.u,
    v_atm  = atmosphere_fields.v,
    T_atm  = atmosphere_fields.T,
    p_atm  = atmosphere_fields.p,
    q_atm  = atmosphere_fields.q,
    Qs_atm = atmosphere_fields.Qs,
    Qℓ_atm = atmosphere_fields.Qℓ,
    Mp_atm = atmosphere_fields.Mp,
))
prefix = joinpath(sim_setup.results_dir, "snapshots_atmosphere")
ocean_sim.output_writers[:atmosphere] = JLD2OutputWriter(
    ocean_model,
    atmosphere_data;
    schedule = TimeInterval(1hours),
    filename = "$prefix.jld2",
    overwrite_existing = true,
    array_type = Array{Float32},
)

## Running the simulation
# We use an adaptive time step that maintains the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)
# We set the CFL in the time_step_wizard to be 0.25 as this is the maximum recommended CFL to be
# used in conjunction with Oceananigans' hydrostatic time-stepping algorithm ([two step Adams-Bashfort](https://en.wikipedia.org/wiki/Linear_multistep_method))
ocean_sim.stop_time = 365days
coupled_simulation.stop_time = 365days

conjure_time_step_wizard!(coupled_simulation; cfl = 0.25, max_Δt = 10minutes, max_change = 1.01)
run!(coupled_simulation)
