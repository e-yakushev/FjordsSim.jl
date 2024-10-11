using Oceananigans
using JLD2
using Oceananigans.Units

include("../../src/FjordsSim.jl")

using .FjordsSim: plot_1d_phys, extract_z_faces, record_vertical_tracer, record_surface_speed

folder = joinpath(homedir(), "data_Varna")
filename = joinpath(folder, "snapshots270days")
T =   FieldTimeSeries("$filename.jld2", "T")
S =   FieldTimeSeries("$filename.jld2", "S")
u =   FieldTimeSeries("$filename.jld2", "u")
v =   FieldTimeSeries("$filename.jld2", "v")
O₂ =  FieldTimeSeries("$filename.jld2", "O₂")
NUT =  FieldTimeSeries("$filename.jld2", "NUT")
PHY =  FieldTimeSeries("$filename.jld2", "P") 
times = T.times

grid = jldopen("$filename.jld2")["grid"]
# z = extract_z_faces(grid)

# # plot_1d_phys(T, S, z, times, folder)
# record_surface_speed(u, v, 10, times, folder)

# record_vertical_tracer(T, 18, times, folder, "Tprofile", "temperature")
# record_vertical_tracer(
#     u, 18, times, folder, "Uprofile", "u component",
#     colorrange=(-0.5, 0.5), colormap=:deep,
#     )

# record_vertical_tracer(
#     S, 18, times, folder, "Sprofile", "S",
#     colorrange=(0, 20), colormap=:viridis,
#         )

record_vertical_tracer(
    PHY, 18, times, folder, "PHYprofile", "PHY",
    colorrange=(0, 10), colormap=:viridis,
        )

# record_vertical_tracer(
#     NUT, 18, times, folder, "NUTprofile", "NUT",
#     colorrange=(0, 20), colormap=:viridis,
#         )

# record_vertical_tracer(
#     O₂, 18, times, folder, "O2profile", "O2",
#     colorrange=(100, 500), colormap=:viridis,
#         )