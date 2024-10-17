using Oceananigans
using JLD2
using Oceananigans.Units
using CairoMakie

include("../../src/FjordsSim.jl")

using .FjordsSim: plot_1d_phys, extract_z_faces, record_vertical_tracer, record_surface_speed, record_surface_tracer, plot_ztime

Nz = 10

folder = joinpath(homedir(), "data_Varna")
filename = joinpath(folder, "snapshots")
T =   FieldTimeSeries("$filename.jld2", "T")
S =   FieldTimeSeries("$filename.jld2", "S")
u =   FieldTimeSeries("$filename.jld2", "u")
v =   FieldTimeSeries("$filename.jld2", "v")
O₂ =  FieldTimeSeries("$filename.jld2", "O₂")
NUT =  FieldTimeSeries("$filename.jld2", "NUT")
PHY =  FieldTimeSeries("$filename.jld2", "P")
HET =  FieldTimeSeries("$filename.jld2", "HET")
DOM =  FieldTimeSeries("$filename.jld2", "DOM")
POM =  FieldTimeSeries("$filename.jld2", "POM")
C =  FieldTimeSeries("$filename.jld2", "C")       
times = T.times

PHY
grid = jldopen("$filename.jld2")["grid"]

println(keys(grid["underlying_grid"]))
println(grid["underlying_grid"]["Δyᶠᶜᵃ"])

# stupid, but I cannot find a right way with znodes
# znodes(grid["underlying_grid"], with_halos=false)
z = grid["underlying_grid"]["zᵃᵃᶜ"][8:17]

# z = extract_z_faces(grid)

plot_ztime(PHY, HET, POM, DOM, NUT, O₂, T, S, 84, 14, times, z, folder)

# HORIZONTAL
# # plot_1d_phys(T, S, z, times, folder)
record_surface_speed(u, v, 10, times, folder)

record_surface_tracer(
    C, Nz, times, folder, "Contsurf", "Contaminant",
    colorrange=(0, 100), colormap=:matter,
    )

record_surface_tracer(
    T, Nz, times, folder, "Tsurf", "Temperature",
    colorrange=(5, 21), colormap=Reverse(:RdYlBu),
    )

record_surface_tracer(
    S, Nz, times, folder, "Ssurf", "Salinity",
    colorrange=(0, 17), colormap=:viridis,
    )

record_surface_tracer(
    O₂, Nz, times, folder, "O2surf", "Oxy",
    colorrange=(100, 350), colormap=:turbo,
    )

record_surface_tracer(
    NUT, Nz, times, folder, "NUTsurf", "NUT",
    colorrange=(0, 20), colormap=Reverse(:cherry),
    )

record_surface_tracer(
    PHY, Nz, times, folder, "PHYsurf", "PHY",
    colorrange=(0, 2), colormap=Reverse(:cubehelix),
    )

# VERTICAL
record_vertical_tracer(
    T, 18, times, folder, "Tprofile", "temperature",
    colorrange=(5, 21), colormap=Reverse(:RdYlBu),
    )

record_vertical_tracer(
    u, 18, times, folder, "Uprofile", "u component",
    colorrange=(-0.5, 0.5), colormap=:deep,
    )

record_vertical_tracer(
    S, 18, times, folder, "Sprofile", "S",
    colorrange=(0, 17), colormap=:viridis,
        )

record_vertical_tracer(
    PHY, 18, times, folder, "PHYprofile", "PHY",
    colorrange=(0, 2), colormap=Reverse(:cubehelix),
        )

record_vertical_tracer(
    DOM, 18, times, folder, "DOMprofile", "DOM",
    colorrange=(0, 20), colormap=Reverse(:CMRmap),
        )

record_vertical_tracer(
    POM, 18, times, folder, "POMprofile", "POM",
    colorrange=(0, 100), colormap=Reverse(:CMRmap),
        )

record_vertical_tracer(
    NUT, 18, times, folder, "NUTprofile", "NUT",
    colorrange=(0, 10), colormap=Reverse(:cherry),
        )

record_vertical_tracer(
    O₂, 18, times, folder, "O2profile", "O2",
    colorrange=(100, 350), colormap=:turbo,
        )