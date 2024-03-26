using Oceananigans
using GLMakie

filepath = joinpath(homedir(), "fjords_data", "oslo_fjord_snapshots.jld2")
filename = joinpath(homedir(), "fjords_data", "oslo_fjord_snapshots")

time_series = FieldTimeSeries(filepath, "T")
x, y, z = nodes(time_series)

# get a land mask
h = interior(time_series.grid.immersed_boundary.bottom_height)[:, :, 1]
land = h .>= -1

times = time_series.times
intro = 1
n = Observable(intro)

Txy = @lift begin
    # time_series[1] takes the first time step (it is the last dimension)
    # 3 dims left: (x, y, z), we take (x, y, end) from interior grid points
    snap = interior(time_series[$n],  :, :, size(time_series[1])[end])
    snap[land] .= NaN
    snap
end

Txz = @lift interior(time_series[$n],  :, 150, :)

fig = Figure(size = (1000, 500))

grid = time_series.grid.underlying_grid
axis_kwargs_xy = (
    xlabel = "x",
    ylabel = "y",
    aspect = AxisAspect(grid.Lx/grid.Ly),
    limits = ((0, grid.Ly), (-grid.Ly, 0))
)
axis_kwargs_xz = (
    xlabel = "x",
    ylabel = "z",
    aspect = AxisAspect(grid.Lx/grid.Lz),
    limits = ((0, grid.Lx), (-grid.Lz, 0))
)

ax_xy = Axis(fig[1, 1])  # ; axis_kwargs_xy...)
ax_xz = Axis(fig[1, 3])  # ; axis_kwargs_xz...)

Tlims = (15, 20)

hm_xy = heatmap!(ax_xy, x, y, Txy; nan_color=:white, colormap = :thermal)  # , colorrange = Tlims)
Colorbar(fig[1, 2], hm_xy; label = "m")

hm_xz = heatmap!(ax_xz, x, z, Txz; nan_color=:white, colormap = :thermal)  # , colorrange = Tlims)
Colorbar(fig[1, 4], hm_xz; label = "m")

# display(fig)

frames = intro:length(times)
record(fig, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end
