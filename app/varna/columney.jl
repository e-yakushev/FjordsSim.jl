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

## Model setup
using OceanBioME, Oceananigans, Printf
using OceanBioME: GasExchange
using OceanBioME.Sediments: sinking_flux
#using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units
using Oceananigans: Forcing
using Interpolations
using JLD2
using CairoMakie
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

include("setup.jl")

using .FjordsSim: OXYDEP, read_TSU_forcing

const year = 365days
stoptime = 1095days  # Set simulation stoptime here!

## Grid
#depth_extent = 100meters
Nz = 12
grid = RectilinearGrid(
    size = (1, 1, Nz),
    extent = (500meters, 500meters, 67meters),
    topology = (Bounded, Bounded, Bounded),
)

## Model
add_contaminants = false

biogeochemistry = OXYDEP(;
    grid,
    args_oxydep...,
    surface_photosynthetically_active_radiation = PAR⁰,
    TS_forced = true,
    Chemicals = add_contaminants,
    scale_negatives = true,
)
biogeochemistry
## Hydrophysics forcing
filename = "../../data_Varna/Varna_brom.nc"

Tnc, Snc, Unc, Kznc, depth, times = read_TSU_forcing(filename)
Kznc = 8.0 * Kznc #10
Kznc[:, 1] = Kznc[:, 1] ./ 10.0  # we decrease Kz above the bottom

# restore z-faces from nc file, as it provides us only centers of layers. dz=5
# z-faces are needed to construct input_grid
z_faces = depth .+ 2.6
z_faces

times = collect(range(0, stop = 366 * 24 * 3600, step = 3600))[1:8784]
temp_itp = interpolate((times, z_faces), Tnc, Gridded(Linear()))
sal_itp = interpolate((times, z_faces), Snc, Gridded(Linear()))
kz_itp = interpolate((times, z_faces), Kznc, Gridded(Linear()))

# Define a function to perform bilinear interpolation
function bilinear_interpolate(itp, t, z)
    return itp(t, z)
end

T_function(x, y, z, t) = bilinear_interpolate(temp_itp, mod(t, 365days), z)
S_function(x, y, z, t) = bilinear_interpolate(sal_itp, mod(t, 365days), z)
Kz_function(x, y, z, t) = bilinear_interpolate(kz_itp, mod(t, 365days), clamp(z, -67, 0))

clock = Clock(; time = times[1])

T = FunctionField{Center,Center,Center}(T_function, grid; clock)
S = FunctionField{Center,Center,Center}(S_function, grid; clock)

κ = 5.0 * FunctionField{Center,Center,Center}(Kz_function, grid; clock)

#- - - - - - - - - - - - - - - - - - - - - - 
## Boundary conditions for OxyDep
O2_suboxic = 30.0  # OXY threshold for oxic/suboxic switch (mmol/m3)
Trel = 10000.0     # Relaxation time for exchange with the sediments (s/m)
b_ox = 15.0        # difference of OXY in the sediment and water, 
b_NUT = 10.0       # NUT in the sediment, (mmol/m3)  18
b_DOM_ox = 6.0     # OM in the sediment (oxic conditions), (mmol/m3) 
b_DOM_anox = 20.0  # OM in the sediment (anoxic conditions), (mmol/m3)  
bu = 0.85 #0.2 0.8=hyp     # Burial coeficient for lower boundary (0<Bu<1), 1 - for no burying, (nd)

windspeed = 5.0    # m/s windspeed for gases exchange

@inline F_ox(conc, threshold) = (0.5 + 0.5 * tanh(conc - threshold))
@inline F_subox(conc, threshold) = (0.5 - 0.5 * tanh(conc - threshold))

## oxy
include("sea_water_flux.jl")

Oxy_top_cond(i, j, grid, clock, fields) = @inbounds (OxygenSeaWaterFlux(
    fields.T[i, j, Nz],
    fields.S[i, j, Nz],
    0.0,                # sea surface pressure
    fields.O₂[i, j, Nz],
    windspeed,
))
OXY_top = FluxBoundaryCondition(Oxy_top_cond; discrete_form = true)

OXY_bottom_cond(i, j, grid, clock, fields) = @inbounds (
    -(
        F_ox(fields.O₂[i, j, 1], O2_suboxic) * b_ox +
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * (0.0 - fields.O₂[i, j, 1])
    ) / Trel
)
OXY_bottom = FluxBoundaryCondition(OXY_bottom_cond, discrete_form = true)

## nut
#NUT_bottom_cond(i, j, grid, clock, fields) =
NUT_bottom_cond(i, j, grid, clock, fields) = @inbounds (
    (
        F_ox(fields.O₂[i, j, 1], O2_suboxic) * (b_NUT - fields.NUT[i, j, 1]) +
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * (0.0 - fields.NUT[i, j, 1])
    ) / Trel
)
NUT_bottom = FluxBoundaryCondition(NUT_bottom_cond, discrete_form = true)

## phy
w_PHY = biogeochemical_drift_velocity(biogeochemistry, Val(:PHY)).w[1, 1, 1]
PHY_bottom_cond(i, j, grid, clock, fields) = @inbounds (-bu * w_PHY * fields.PHY[i, j, 1])
PHY_bottom = FluxBoundaryCondition(PHY_bottom_cond, discrete_form = true)

## het
w_HET = biogeochemical_drift_velocity(biogeochemistry, Val(:HET)).w[1, 1, 1]
HET_bottom_cond(i, j, grid, clock, fields) = @inbounds (-bu * w_HET * fields.HET[i, j, 1])
HET_bottom = FluxBoundaryCondition(HET_bottom_cond, discrete_form = true)

## pom
w_POM = biogeochemical_drift_velocity(biogeochemistry, Val(:POM)).w[1, 1, 1]
POM_bottom_cond(i, j, grid, clock, fields) = @inbounds (-bu * w_POM * fields.POM[i, j, 1])
POM_bottom = FluxBoundaryCondition(POM_bottom_cond, discrete_form = true)

## dom
DOM_top = ValueBoundaryCondition(0.0)
DOM_bottom_cond(i, j, grid, clock, fields) = @inbounds (
    (
        F_ox(fields.O₂[i, j, 1], O2_suboxic) * (b_DOM_ox - fields.DOM[i, j, 1]) +
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * 2.0 * (b_DOM_anox - fields.DOM[i, j, 1])
    ) / Trel
)
DOM_bottom = FluxBoundaryCondition(DOM_bottom_cond, discrete_form = true)
#- - - - - - - - - - - - - - - - - - - - - - 

## Model instantiation
model = NonhydrostaticModel(;
    grid,
    clock,
    #closure = VerticallyImplicitTimeDiscretization(), #SmagorinskyLilly(), 
    closure = ScalarDiffusivity(ν = κ, κ = κ), #(ν = 1e-4, κ = 1e-4),
    biogeochemistry,
    boundary_conditions = (
        O₂ = FieldBoundaryConditions(top = OXY_top, bottom = OXY_bottom),
        NUT = FieldBoundaryConditions(bottom = NUT_bottom),
        DOM = FieldBoundaryConditions(top = DOM_top, bottom = DOM_bottom),
        POM = FieldBoundaryConditions(bottom = POM_bottom),
        PHY = FieldBoundaryConditions(bottom = PHY_bottom),
        HET = FieldBoundaryConditions(bottom = HET_bottom),
    ),
    auxiliary_fields = (; S, T),
)

model

## Set model
if add_contaminants == false
    set!(model, NUT = 10.0, PHY = 0.01, HET = 0.05, O₂ = 350.0, DOM = 1.0)
else
    set!(model, NUT = 10.0, PHY = 0.01, HET = 0.05, O₂ = 350.0, DOM = 1.0, Ci_free = 0.123)
end

## Simulation
simulation = Simulation(model, Δt = 6minutes, stop_time = stoptime)
progress_message(sim) = @printf(
    "Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
    iteration(sim),
    prettytime(sim),
    prettytime(sim.Δt),
    prettytime(sim.run_wall_time),
)

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))

#add_contaminants ? (Ci_free, NUT, PHY, HET, POM, DOM, O₂ = model.tracers) : (NUT, PHY, HET, POM, DOM, O₂ = model.tracers)
if add_contaminants == false
    NUT, PHY, HET, POM, DOM, O₂ = model.tracers
else
    NUT, PHY, HET, POM, DOM, O₂, Ci_free, Ci_PHY, Ci_HET, Ci_POM, Ci_DOM = model.tracers
end
PAR = model.auxiliary_fields.PAR
T = model.auxiliary_fields.T
S = model.auxiliary_fields.S

output_prefix = joinpath(homedir(), "data_Varna", "columney_snapshots")
# output_prefix = joinpath("out")
if add_contaminants == false
    simulation.output_writers[:profiles] = JLD2OutputWriter(
        model,
        (; NUT, PHY, HET, POM, DOM, O₂, T, S, PAR, κ),
        filename = "$output_prefix.jld2",
        schedule = TimeInterval(1day),
        overwrite_existing = true,
    )
else
    simulation.output_writers[:profiles] = JLD2OutputWriter(
        model,
        (; NUT, PHY, HET, POM, DOM, O₂, T, S, PAR, κ, Ci_free, Ci_PHY, Ci_HET, Ci_POM, Ci_DOM),
        filename = "$output_prefix.jld2",
        schedule = TimeInterval(1day),
        overwrite_existing = true,
    )
end
## Run!
run!(simulation)

## Make plots.
model

include("images.jl")

println(" BOT XPEHb, OCTAHOBKA...")
