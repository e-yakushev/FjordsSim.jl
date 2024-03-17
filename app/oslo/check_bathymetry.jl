using FjordsSim
using FjordsSim: ImmersedBoundaryGrid
using GLMakie
using Oceananigans

fjords_setup = FjordsSetup(
    CPU(),
    (50, 50, 10),
    (58.8, 59.9),
    (10.1, 11.1),
    (-500, 0),
    (4, 4, 4),
    joinpath(homedir(), "fjords_data"),
    "ETOPO_2022_v1_15s_N60E000_surface.nc",
    0,
    5,
)

grid = ImmersedBoundaryGrid(fjords_setup)
λ, φ, z = nodes(grid.underlying_grid, Center(), Center(), Center())

h = grid.immersed_boundary.bottom_height
land = interior(h) .>= 0
interior(h)[land] .= NaN

fig = Figure(; size=(700, 700))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, λ, φ, interior(h, :, :, 1), nan_color=:white, colorrange=(-1000, 0))
Colorbar(fig[1, 2], hm; label = "m")

display(fig)
