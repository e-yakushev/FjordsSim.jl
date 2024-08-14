using Oceananigans
using Oceananigans.Units: hour, days, meters
using JLD2
using CairoMakie

# 1
using OceanBioME, Oceananigans, Printf
using OceanBioME: Boundaries, GasExchange
using OceanBioME.Boundaries.Sediments: sinking_flux
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

include("../../src/Oxydep.jl")
using .OXYDEPModel
#/1

year=365days
stoptime = 365
#depth_extent = 10meters
#grid = RectilinearGrid(size = (1, 1, 10), extent = (500meters, 500meters, 10meters), topology = (Bounded, Bounded, Bounded))
grid = RectilinearGrid(size = (1, 1, 12), extent = (500meters, 500meters, 67meters), topology = (Bounded, Bounded, Bounded))


# filename = joinpath(homedir(), "data_Varna", "columney_snapshots")
filename = "out"
## Load saved output
@info "Loading saved outputs..."
PHY = FieldTimeSeries("$filename.jld2", "PHY")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
HET = FieldTimeSeries("$filename.jld2", "HET")
POM = FieldTimeSeries("$filename.jld2", "POM")
DOM = FieldTimeSeries("$filename.jld2", "DOM")
O₂ = FieldTimeSeries("$filename.jld2", "O₂")
T = FieldTimeSeries("$filename.jld2", "T")
S =   FieldTimeSeries("$filename.jld2", "S")
PAR = FieldTimeSeries("$filename.jld2", "PAR")

@info "Saved outputs loaded..."

z = jldopen("$filename.jld2")["grid"]["zᵃᵃᶜ"]
times = T.times
#times = StepRangeLen(0.0, 1, (160 * 86400))
# 2
biogeochemistry =
    OXYDEP(; grid, particles = nothing)
model = NonhydrostaticModel(; grid,
       biogeochemistry,
)
#/2

nitrogen_burying = zeros(length(times))
POM_burying = zeros(length(times))
PHY_burying = zeros(length(times))
HET_burying = zeros(length(times))

for (i, t) in enumerate(times)
    nitrogen_burying[i] = (
        POM[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:POM)).w[1, 1, 1] +
        PHY[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:PHY)).w[1, 1, 1] +
        HET[1, 1, 1, i] * 
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:HET)).w[1, 1, 1]
    )
    POM_burying[i] = (
        POM[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:POM)).w[1, 1, 1]
    )
    PHY_burying[i] = (
        POM[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:PHY)).w[1, 1, 1]
    )
    HET_burying[i] = (
        POM[1, 1, 1, i] *
        biogeochemical_drift_velocity(model.biogeochemistry, Val(:HET)).w[1, 1, 1]
    )
end


## Plot

fig = Figure(size = (1500, 1000), fontsize = 20)

axis_kwargs = (
    xlabel = "Time (days)",
    ylabel = "z (m)",
 #   limits = ((0, times[end] / days), (-(depth_extent + 10), 10)),
 #   xticks = collect(0:10:stoptime),
)

axPHY = Axis(fig[1, 3]; title = "PHY, mmolN/m³", axis_kwargs...)
hmPHY = heatmap!(times / days, z, interior(PHY, 1, 1, :, :)', colormap = Reverse(:cubehelix)) #(:davos10))
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

axPAR = Axis(fig[1, 5]; title = "PAR  μE⋅m-2⋅s-1", axis_kwargs...)
hmPAR = heatmap!(times / days, z, interior(PAR, 1, 1, :, :)', colormap = :grayC100) # :linear_grey_0_100_c0_n256)
Colorbar(fig[1, 6], hmPAR)

axT = Axis(fig[2, 5]; title = "T, oC", axis_kwargs...)
hmT = heatmap!(times / days, z, interior(T, 1, 1, :, :)', colormap = Reverse(:RdYlBu))
Colorbar(fig[2, 6], hmT)

axS = Axis(fig[3, 5]; title = "S, psu", axis_kwargs...)
hmS = heatmap!(times / days, z, interior(S, 1, 1, :, :)', colormap = Reverse(:RdYlBu))
Colorbar(fig[3, 6], hmS)

@info "Z-Time plots made"

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


@info "Fluxes plots made"

save("out_fluxes.png", fig2)

