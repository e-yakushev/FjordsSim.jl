using Oceananigans
using JLD2
using Oceananigans.Units
using CairoMakie

include("../../src/FjordsSim.jl")

using .FjordsSim: plot_1d_phys, extract_z_faces, record_vertical_tracer, record_surface_speed, record_horizontal_tracer, plot_ztime

Nz = 20

folder = joinpath(homedir(), "FjordsSim_results", "sognefjord")
filename = raw"D:\sogn_snapshots"
i = 50
j = 14

T =   FieldTimeSeries("$filename.jld2", "T"; backend = OnDisk())
S =   FieldTimeSeries("$filename.jld2", "S"; backend = OnDisk())
u =   FieldTimeSeries("$filename.jld2", "u"; backend = OnDisk())
v =   FieldTimeSeries("$filename.jld2", "v"; backend = OnDisk())
O₂ =  FieldTimeSeries("$filename.jld2", "O₂"; backend = OnDisk())
NUT =  FieldTimeSeries("$filename.jld2", "NUT"; backend = OnDisk())
PHY =  FieldTimeSeries("$filename.jld2", "P"; backend = OnDisk())
HET =  FieldTimeSeries("$filename.jld2", "HET"; backend = OnDisk())
DOM =  FieldTimeSeries("$filename.jld2", "DOM"; backend = OnDisk())
POM =  FieldTimeSeries("$filename.jld2", "POM"; backend = OnDisk())
C =  FieldTimeSeries("$filename.jld2", "C"; backend = OnDisk())       
times = T.times

grid = jldopen("$filename.jld2")["grid"]

println(keys(grid["underlying_grid"]))
println(grid["underlying_grid"]["Δyᶠᶜᵃ"])

# stupid, but I cannot find a right way with znodes
# znodes(grid["underlying_grid"], with_halos=false)
z = grid["underlying_grid"]["zᵃᵃᶜ"][8:Nz+7]

# z = extract_z_faces(grid)

# doesnt work with OnDisk
# plot_ztime(PHY, HET, POM, DOM, NUT, O₂, T, S, i, j, times, z, folder)

# HORIZONTAL
# plot_1d_phys(T, S, z, times, folder)

record_surface_speed(u, v, Nz, times, folder)

record_horizontal_tracer(
    C, times, folder, "Contsurf", "Contaminant (% of max. concentration)",
    colorrange=(0, 100), colormap=:matter, iz=Nz,
    )

record_horizontal_tracer(
    T, times, folder, "Tsurf", "Temperature (°C)",
    colorrange=(-1, 10), colormap=Reverse(:RdYlBu), iz=Nz,
    )

record_horizontal_tracer(
    S, times, folder, "Ssurf", "Salinity (PSU)", iz=Nz,
    colorrange=(0, 17), colormap=:viridis,
    )

record_horizontal_tracer(
    O₂, times, folder, "O2bottom", "Dissolved oxygen (μM)",
    colorrange=(100, 350), colormap=:turbo, iz=3,
    )

record_horizontal_tracer(
    NUT, times, folder, "NUTsurf", "Nutrients (μM N)",
    colorrange=(0, 20), colormap=Reverse(:cherry), iz=Nz,
    )

record_horizontal_tracer(
    PHY, times, folder, "PHYsurf", "Phytoplankton (μM N)",
    colorrange=(0, 0.5), colormap=Reverse(:cubehelix), iz=Nz,
    )

# VERTICAL
record_vertical_tracer(
    T, z, 18, times, folder, "Tprofile", "Temperature (°C)",
    colorrange=(5, 21), colormap=Reverse(:RdYlBu),
    )

record_vertical_tracer(
    u, z, 18, times, folder, "Uprofile", "u velocity component (ms⁻¹)",
    colorrange=(-0.5, 0.5), colormap=:deep,
    )

record_vertical_tracer(
    S, z, 18, times, folder, "Sprofile", "Salinity (PSU)",
    colorrange=(0, 17), colormap=:viridis,
        )

record_vertical_tracer(
    PHY, z, 18, times, folder, "PHYprofile", "Phytoplankton (μM N)",
    colorrange=(0, 2), colormap=Reverse(:cubehelix),
        )

record_vertical_tracer(
    DOM, z, 18, times, folder, "DOMprofile", "DOM (μM N)",
    colorrange=(0, 20), colormap=Reverse(:CMRmap),
        )

record_vertical_tracer(
    POM, z, 18, times, folder, "POMprofile", "POM (μM N)",
    colorrange=(0, 100), colormap=Reverse(:CMRmap),
        )

record_vertical_tracer(
    NUT, z, 18, times, folder, "NUTprofile", "Nutrients (μM N)",
    colorrange=(0, 10), colormap=Reverse(:cherry),
        )

record_vertical_tracer(
    O₂, z, 18, times, folder, "O2profile", "Dissolved Oxygen (μM)",
    colorrange=(100, 350), colormap=:turbo,
        )
