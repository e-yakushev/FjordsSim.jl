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

using .FjordsSim: OXYDEP, read_TSU_forcing, OxygenSeaWaterFlux

const year = 365days
stoptime = 1095days  # Set simulation stoptime here!

## Grid
#depth_extent = 100meters
Nz = 12
grid_column = RectilinearGrid(
   size = (1, 1, Nz),
   extent = (500meters, 500meters, 67meters),
   topology = (Bounded, Bounded, Bounded),
)

## Model
add_contaminants = false

biogeochemistry = OXYDEP(;
    grid=grid_column,
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

temp_itp

T = FunctionField{Center,Center,Center}(T_function, grid_column; clock)
S = FunctionField{Center,Center,Center}(S_function, grid_column; clock)

κ = 5.0 * FunctionField{Center,Center,Center}(Kz_function, grid_column; clock)

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
w_PHY = biogeochemical_drift_velocity(biogeochemistry, Val(:P)).w[1, 1, 1]
PHY_bottom_cond(i, j, grid, clock, fields) = @inbounds (-bu * w_PHY * fields.P[i, j, 1])
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
    grid = grid_column,
    clock,
    #closure = VerticallyImplicitTimeDiscretization(), #SmagorinskyLilly(), 
    closure = ScalarDiffusivity(ν = κ, κ = κ), #(ν = 1e-4, κ = 1e-4),
    biogeochemistry,
    boundary_conditions = (
        O₂ = FieldBoundaryConditions(top = OXY_top, bottom = OXY_bottom),
        NUT = FieldBoundaryConditions(bottom = NUT_bottom),
        DOM = FieldBoundaryConditions(top = DOM_top, bottom = DOM_bottom),
        POM = FieldBoundaryConditions(bottom = POM_bottom),
        P = FieldBoundaryConditions(bottom = PHY_bottom),
        HET = FieldBoundaryConditions(bottom = HET_bottom),
    ),
    auxiliary_fields = (; S, T),
)

model

## Set model
if add_contaminants == false
    set!(model, NUT = 10.0, P = 0.01, HET = 0.05, O₂ = 350.0, DOM = 1.0)
else
    set!(model, NUT = 10.0, P = 0.01, HET = 0.05, O₂ = 350.0, DOM = 1.0, Ci_free = 0.123)
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
    NUT, P, HET, POM, DOM, O₂ = model.tracers
else
    NUT, P, HET, POM, DOM, O₂, Ci_free, Ci_PHY, Ci_HET, Ci_POM, Ci_DOM = model.tracers
end
PAR = model.auxiliary_fields.PAR
T = model.auxiliary_fields.T
S = model.auxiliary_fields.S

output_prefix = joinpath(homedir(), "data_Varna", "columney_snapshots")
# output_prefix = joinpath("out")
if add_contaminants == false
    simulation.output_writers[:profiles] = JLD2OutputWriter(
        model,
        (; NUT, P, HET, POM, DOM, O₂, T, S, PAR, κ),
        filename = "$output_prefix.jld2",
        schedule = TimeInterval(1day),
        overwrite_existing = true,
    )
else
    simulation.output_writers[:profiles] = JLD2OutputWriter(
        model,
        (; NUT, P, HET, POM, DOM, O₂, T, S, PAR, κ, Ci_free, Ci_PHY, Ci_HET, Ci_POM, Ci_DOM),
        filename = "$output_prefix.jld2",
        schedule = TimeInterval(1day),
        overwrite_existing = true,
    )
end
## Run!
run!(simulation)

## Make plots.
model

"
Plotting of images
"
filename = joinpath(homedir(), "FjordsSim_results", "columney_snapshots")
 
## Load saved output
@info "Loading saved outputs..."
P = FieldTimeSeries("$filename.jld2", "P")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
HET = FieldTimeSeries("$filename.jld2", "HET")
POM = FieldTimeSeries("$filename.jld2", "POM")
DOM = FieldTimeSeries("$filename.jld2", "DOM")
O₂ =  FieldTimeSeries("$filename.jld2", "O₂")
T =   FieldTimeSeries("$filename.jld2", "T")
S =   FieldTimeSeries("$filename.jld2", "S")
PAR = FieldTimeSeries("$filename.jld2", "PAR")
κ =   FieldTimeSeries("$filename.jld2", "κ")
#Ci_free =   FieldTimeSeries("$filename.jld2", "Ci_free")
@info "Saved outputs loaded..."

z = jldopen("$filename.jld2")["grid"]["zᵃᵃᶜ"]
times = T.times
#times = StepRangeLen(0.0, 1, (160 * 86400))
# 2
#biogeochemistry =
#    OXYDEP(; grid, particles = nothing)
#model = NonhydrostaticModel(; grid,
#       biogeochemistry,
#)
#/2

nitrogen_burying = zeros(length(times))
POM_burying = zeros(length(times))
PHY_burying = zeros(length(times))
HET_burying = zeros(length(times))

for (i, t) in enumerate(times)
    nitrogen_burying[i] = (
        POM[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:POM)).w[1, 1, 1] +
        P[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:P)).w[1, 1, 1] +
        HET[1, 1, 1, i] * 
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:HET)).w[1, 1, 1]
    )
    POM_burying[i] = (
        POM[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:POM)).w[1, 1, 1]
    )
    PHY_burying[i] = (
        P[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:P)).w[1, 1, 1]
    )
    HET_burying[i] = (
        HET[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:HET)).w[1, 1, 1]
    )
end


## Plot

fig = Figure(size = (1500, 1000), fontsize = 20)

axis_kwargs = (
    xlabel = "Time (days)",
    ylabel = "z (m)",
 #   limits = ((0, times[end] / days), (-(depth_extent + 10), 10)),
     xticks = (0:365:stoptime),
     xtickformat = "{:.0f}" #   values -> ["$(value)kg" for value in values]     
 #    xtickformat = x -> @sprintf("%.1f", x)
 #   xticks = collect(0:stoptime, 365),
)

axPHY = Axis(fig[1, 3]; title = "PHY, mmolN/m³", axis_kwargs...)
hmPHY = heatmap!(times / days, z, interior(P, 1, 1, :, :)', colormap = Reverse(:cubehelix)) #(:davos10))
Colorbar(fig[1, 4], hmPHY)

axHET = Axis(fig[2, 3]; title = "HET, mmolN/m³", axis_kwargs...)
hmHET = heatmap!(times / days, z, interior(HET, 1, 1, :, :)', colormap = Reverse(:afmhot))
Colorbar(fig[2, 4], hmHET)

axPOM = Axis(fig[3, 3]; title = "POM, mmolN/m³", axis_kwargs...)
hmPOM =
    heatmap!(times / days, z, interior(POM, 1, 1, :, :)', colormap = Reverse(:greenbrownterrain)) #(:bilbao25))
hmPOM =
    heatmap!(times / days, z, interior(POM, 1, 1, :, :)', colormap = Reverse(:greenbrownterrain)) #(:bilbao25))
Colorbar(fig[3, 4], hmPOM)

axDOM = Axis(fig[3, 1]; title = "DOM, mmolN/m³", axis_kwargs...)
hmDOM = heatmap!(times / days, z, interior(DOM, 1, 1, :, :)', colormap = Reverse(:CMRmap)) #(:devon10))
Colorbar(fig[3, 2], hmDOM)

axNUT = Axis(fig[1, 1]; title = "NUT, mmolN/m³", axis_kwargs...)
hmNUT = heatmap!(times / days, z, interior(NUT, 1, 1, :, :)', colormap = Reverse(:cherry))
hmNUT = heatmap!(times / days, z, interior(NUT, 1, 1, :, :)', colormap = Reverse(:cherry))
Colorbar(fig[1, 2], hmNUT)

axOXY = Axis(fig[2, 1]; title = "OXY, mmol/m³", axis_kwargs...)
hmOXY = heatmap!(times / days, z, interior(O₂, 1, 1, :, :)', colormap = :turbo)
hmOXY = heatmap!(times / days, z, interior(O₂, 1, 1, :, :)', colormap = :turbo)
Colorbar(fig[2, 2], hmOXY)

axκ = Axis(fig[1, 5]; title = "κ  m³/s", axis_kwargs...)
hmκ = heatmap!(times / days, z, interior(κ, 1, 1, :, :)', colormap = Reverse(:RdYlBu)) # :linear_grey_0_100_c0_n256)
Colorbar(fig[1, 6], hmκ)

axT = Axis(fig[2, 5]; title = "T, oC", axis_kwargs...)
hmT = heatmap!(times / days, z, interior(T, 1, 1, :, :)', colormap = Reverse(:RdYlBu))
Colorbar(fig[2, 6], hmT)

axS = Axis(fig[3, 5]; title = "S, psu", axis_kwargs...)
hmS = heatmap!(times / days, z, interior(S, 1, 1, :, :)', colormap = Reverse(:RdYlBu))
Colorbar(fig[3, 6], hmS)

axPAR = Axis(fig[4, 1]; title = "PAR  μE⋅m-2⋅s-1", axis_kwargs...)
hmPAR = heatmap!(times / days, z, interior(PAR, 1, 1, :, :)', colormap = :grayC100) # :linear_grey_0_100_c0_n256)
Colorbar(fig[4, 2], hmPAR)

@info "VARIABLES Z-Time plots made"

save("out_vars.png", fig)

fig2 = Figure(size = (1500, 1000), fontsize = 20)

axNburying = Axis(
    fig2[1, 1],
    xlabel = "Time (days)",
    ylabel = "Flux (mmolN/m²/year)",
    title = "N burying",
    limits = ((0, times[end] / days), nothing),
)
lines!(axNburying, times / days, nitrogen_burying / 1e3 * year, linewidth = 3, label = "N burying")

axPOMburying = Axis(
    fig2[1, 2],
    xlabel = "Time (days)",
    ylabel = "Flux (mmolN/m²/year)",
    title = "POM burying",
    limits = ((0, times[end] / days), nothing),
)
lines!(axPOMburying, times / days, POM_burying / 1e3 * year, linewidth = 3, color = :brown,  label = "POM burying")

axPHYburying = Axis(
    fig2[2, 1],
    xlabel = "Time (days)",
    ylabel = "Flux (mmolN/m²/year)",
    title = "PHY burying",
    limits = ((0, times[end] / days), nothing),
)
lines!(axPHYburying, times / days, PHY_burying / 1e3 * year, linewidth = 3, color = :red,  label = "PHYburying")

axHETburying = Axis(
    fig2[2, 2],
    xlabel = "Time (days)",
    ylabel = "Flux (mmolN/m²/year)",
    title = "HET burying",
    limits = ((0, times[end] / days), nothing),
)
lines!(axHETburying, times / days, HET_burying / 1e3 * year, linewidth = 3, color = :black, label = "HET burying")


@info "FLUXES plots made"

save("out_fluxes.png", fig2)

println(" BOT XPEHb, OCTAHOBKA...")
