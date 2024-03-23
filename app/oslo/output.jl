using Oceananigans
using GLMakie

filepath = joinpath(homedir(), "fjords_data", "oslo_fjord_snapshots.jld2")

time_series = FieldTimeSeries(filepath, "T")
x, y, z = nodes(time_series)
snap = interior(time_series)[:, :, 10, 4]
h = interior(time_series.grid.immersed_boundary.bottom_height)[:, :, 1]
land = h .>= 0
snap[land] .= NaN

f = Figure(size = (900, 700))
axis_kwargs = (xlabel="x", ylabel="y")
ax = Axis(f[1, 1]; axis_kwargs...)
hm = heatmap!(ax, x, y, snap; nan_color=:white, colormap = :thermal)
Colorbar(f[1, 2], hm; label = "m")

display(f)
