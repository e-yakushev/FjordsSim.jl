using GLMakie
using Oceananigans

include("../../src/FjordsSim.jl")
include("setup.jl")

using .FjordsSim

setup_grid = SetupGridRegrid(;args_grid...)
grid = ImmersedBoundaryGrid(setup_grid)
λ, φ, z = nodes(grid.underlying_grid, Center(), Center(), Center())

h = interior(grid.immersed_boundary.bottom_height)
land = h .>= 0
h[land] .= NaN
h = h[:, :, 1]

fig = Figure(; size=(900, 700))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, λ, φ, h, nan_color=:white, colorrange=(-500, 0))
Colorbar(fig[1, 2], hm; label = "m")

display(fig)
