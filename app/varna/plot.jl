using Oceananigans
using JLD2
using Oceananigans.Units

include("../../src/FjordsSim.jl")

using .FjordsSim: plot_1d_phys, extract_z_faces

folder = joinpath(homedir(), "FjordsSim_results")
filename = joinpath(folder, "snapshots")
T =   FieldTimeSeries("$filename.jld2", "T")
S =   FieldTimeSeries("$filename.jld2", "S")
times = T.times

grid = jldopen("$filename.jld2")["grid"]
z = extract_z_faces(grid)

plot_1d_phys(T, S, z, times, folder)
