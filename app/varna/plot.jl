using Oceananigans
using JLD2
using Oceananigans.Units

include("../../src/FjordsSim.jl")

using .FjordsSim: record_surface_speed, plot_1d_phys, extract_z_faces

folder = joinpath(homedir(), "FjordsSim_results")
filename = joinpath(folder, "varna_snapshots")

u = FieldTimeSeries("$filename.jld2", "u")
v = FieldTimeSeries("$filename.jld2", "v")
grid = jldopen("$filename.jld2")["grid"]["underlying_grid"]
z = extract_z_faces(grid)
times = u.times

record_surface_speed(u, v, grid["Nz"], times, folder)

T = FieldTimeSeries("$filename.jld2", "T")
S = FieldTimeSeries("$filename.jld2", "S")
times = T.times

plot_1d_phys(T, S, z, times, folder, 104, 28)
