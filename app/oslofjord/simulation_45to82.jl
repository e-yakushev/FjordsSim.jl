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
using Oceananigans.OutputWriters: JLD2Writer, NetCDFWriter
using Oceanostics
using FjordsSim: coupled_hydrostatic_simulation, progress
using Printf

include("setup.jl")

## Model Setup
sim_setup = setup_region_3d_OXYDEP_45to82()

coupled_simulation = coupled_hydrostatic_simulation(sim_setup)

coupled_simulation.callbacks[:progress] = Callback(progress, TimeInterval(12hours))

ocean_sim = coupled_simulation.model.ocean
ocean_model = ocean_sim.model

# ocean model output
net_ocean_fluxes = NamedTuple((
    u_atm_ocean_flux = coupled_simulation.model.interfaces.net_fluxes.ocean_surface.u,
    v_atm_ocean_flux = coupled_simulation.model.interfaces.net_fluxes.ocean_surface.v,
))
prefix = joinpath(sim_setup.results_dir, "snapshots_ocean")

# write as JLD2
#ocean_sim.output_writers[:ocean] = JLD2Writer(
#    ocean_model,
#    merge(
#        ocean_model.tracers,
#        ocean_model.velocities,
#        coupled_simulation.model.interfaces.atmosphere_ocean_interface.fluxes,
#        net_ocean_fluxes,
#    );
#    filename = "$prefix",
#    schedule = TimeInterval(3hours),
#    overwrite_existing = true,
#)

# write as NetCDF
ocean_sim.output_writers[:ocean] = NetCDFWriter(
    ocean_model,
    merge(
        ocean_model.tracers,
        ocean_model.velocities,
        coupled_simulation.model.interfaces.atmosphere_ocean_interface.fluxes,
        net_ocean_fluxes,
    );
    filename = "$prefix",
    schedule = TimeInterval(6hours), #(3hours),
    overwrite_existing = true,
)


# prescribed interpolated atmosphere
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
#prefix = joinpath(sim_setup.results_dir, "snapshots_atmosphere")
#ocean_sim.output_writers[:atmosphere] = NetCDFWriter(
#ocean_sim.output_writers[:atmosphere] = JLD2Writer(    
#    ocean_model,
#    atmosphere_data;
#    filename = "$prefix",
#    schedule = TimeInterval(6hours),
#    overwrite_existing = true,
#)

# interface fluxes
#interface_fluxes = coupled_simulation.model.interfaces.atmosphere_ocean_interface.fluxes
#prefix = joinpath(sim_setup.results_dir, "snapshots_interface_fluxes")
#ocean_sim.output_writers[:interface_fluxes] = NetCDFWriter(
#ocean_sim.output_writers[:interface_fluxes] = JLD2Writer(        
#    ocean_model,
#    interface_fluxes;
#    filename = "$prefix",
#    schedule = TimeInterval(6hours),
#    overwrite_existing = true,
#)

# net fluxes
net_ocean_fluxes = coupled_simulation.model.interfaces.net_fluxes.ocean_surface
net_ocean_fluxes = (u = net_ocean_fluxes.u,
                    v = net_ocean_fluxes.v,
                    T = net_ocean_fluxes.T,
                    S = net_ocean_fluxes.S)
#prefix = joinpath(sim_setup.results_dir, "snapshots_net_ocean_fluxes")
#ocean_sim.output_writers[:net_fluxes] = NetCDFWriter(
#ocean_sim.output_writers[:net_fluxes] = JLD2Writer(         
#    ocean_model,
#    net_ocean_fluxes;
#    filename = "$prefix",
#    schedule = TimeInterval(6hours),
#    overwrite_existing = true,
#)

## Spinning up the simulation
# We use an adaptive time step that maintains the [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition) equal to 0.1.
ocean_sim.stop_time = 10days
coupled_simulation.stop_time = 10days

conjure_time_step_wizard!(coupled_simulation; cfl = 0.1, max_Δt = 1minute, max_change = 1.01)
run!(coupled_simulation)

## Running the simulation
# This time, we set the CFL in the time_step_wizard to be 0.25 as this is the maximum recommended CFL to be
# used in conjunction with Oceananigans' hydrostatic time-stepping algorithm ([two step Adams-Bashfort](https://en.wikipedia.org/wiki/Linear_multistep_method))
ocean_sim.stop_time = 365days
coupled_simulation.stop_time = 365days

conjure_time_step_wizard!(coupled_simulation; cfl = 0.25, max_Δt = 3minutes, max_change = 1.01)
run!(coupled_simulation)
