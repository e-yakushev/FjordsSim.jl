using Oceananigans
using GLMakie

filepath = joinpath(homedir(), "fjords_data", "oslo_fjord_snapshots.jld2")
filename = joinpath(homedir(), "fjords_data", "oslo_fjord_snapshots")

time_series = FieldTimeSeries(filepath, "T")
x, y, z = nodes(time_series)

h = interior(time_series.grid.immersed_boundary.bottom_height)[:, :, 1]
land = h .>= -1

times = time_series.times
intro = 1

n = Observable(intro)
Tₙ = @lift begin
    snap = interior(time_series[$n],  :, :, size(time_series[1])[end])
    snap[land] .= NaN
    snap
end

f = Figure(size = (900, 700))
axis_kwargs = (xlabel="x", ylabel="y")
ax = Axis(f[1, 1]; axis_kwargs...)
hm = heatmap!(ax, x, y, Tₙ; nan_color=:white, colormap = :thermal)
Colorbar(f[1, 2], hm; label = "m")

frames = intro:length(times)
record(f, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end
