using Oceananigans
using JLD2
using Oceananigans.Units

include("../../src/FjordsSim.jl")

using .FjordsSim: plot_1d_phys, extract_z_faces, record_vertical_tracer, record_surface_speed, record_surface_tracer

Nz = 10

folder = joinpath(homedir(), "data_Varna")
filename = joinpath(folder, "snapshots26daysGOOD3")
T =   FieldTimeSeries("$filename.jld2", "T")
S =   FieldTimeSeries("$filename.jld2", "S")
u =   FieldTimeSeries("$filename.jld2", "u")
v =   FieldTimeSeries("$filename.jld2", "v")
O₂ =  FieldTimeSeries("$filename.jld2", "O₂")
NUT =  FieldTimeSeries("$filename.jld2", "NUT")
PHY =  FieldTimeSeries("$filename.jld2", "P")
DOM =  FieldTimeSeries("$filename.jld2", "DOM")  
times = T.times

grid = jldopen("$filename.jld2")["grid"]
# z = extract_z_faces(grid)

# HORIZONTAL
# # plot_1d_phys(T, S, z, times, folder)
record_surface_speed(u, v, 10, times, folder)

record_surface_tracer(
    T, Nz, times, folder, "Tsurf", "Temperature",
    colorrange=(8, 11), colormap=:magma,
    )

record_surface_tracer(
    S, Nz, times, folder, "Ssurf", "Salinity",
    colorrange=(0, 18), colormap=:viridis,
    )

record_surface_tracer(
    O₂, Nz, times, folder, "O2surf", "Oxy",
    colorrange=(300, 360), colormap=:viridis,
    )

record_surface_tracer(
    NUT, Nz, times, folder, "NUTsurf", "NUT",
    colorrange=(0, 10), colormap=:viridis,
    )

# VERTICAL
record_vertical_tracer(T, 18, times, folder, "Tprofile", "temperature")
record_vertical_tracer(
    u, 18, times, folder, "Uprofile", "u component",
    colorrange=(8, 11), colormap=:deep,
    )

record_vertical_tracer(
    S, 18, times, folder, "Sprofile", "S",
    colorrange=(0, 18), colormap=:viridis,
        )

record_vertical_tracer(
    PHY, 18, times, folder, "PHYprofile", "PHY",
    colorrange=(0, 0.2), colormap=:viridis,
        )

record_vertical_tracer(
    DOM, 18, times, folder, "DOMprofile", "DOM",
    colorrange=(0, 2), colormap=:viridis,
        )

record_vertical_tracer(
    NUT, 18, times, folder, "NUTprofile", "NUT",
    colorrange=(0, 10), colormap=:viridis,
        )

record_vertical_tracer(
    O₂, 18, times, folder, "O2profile", "O2",
    colorrange=(200, 400), colormap=:viridis,
        )
