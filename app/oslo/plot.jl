using Oceananigans
using Oceananigans.Units: hour

include("../../src/FjordsSim.jl")

using .FjordsSim: record_surface_speed

folder = joinpath(homedir(), "FjordsSim_results")
filename = joinpath(folder, "oslo_surface_snapshots")
u =   FieldTimeSeries("$filename.jld2", "u")
v =   FieldTimeSeries("$filename.jld2", "v")
times = u.times

record_surface_speed(u, v, times, folder)

x, y, z = nodes(time_series)

# get a land mask
h = interior(time_series.grid.immersed_boundary.bottom_height)[:, :, 1]
land = h .>= -3

times = time_series.times
intro = searchsortedfirst(times, 1hour)
n = Observable(intro)

Txy = @lift begin
    # time_series[1] takes the first time step (it is the last dimension)
    # 3 dims left: (x, y, z), we take (x, y, end) from interior grid points
    snap = interior(time_series[$n],  :, :, size(time_series[1])[end])
    snap[land] .= NaN
    snap
end

Txz = @lift interior(time_series[$n],  :, 50, :)

fig = Figure(size = (1000, 500))

Lxₘᵢₙ, Lxₘₐₓ = x[1], x[end]
Lyₘᵢₙ, Lyₘₐₓ = y[1], y[end]
Lzₘᵢₙ, Lzₘₐₓ = z[1], z[end]

axis_kwargs_xy = (
    xlabel = "x",
    ylabel = "y",
    aspect = AxisAspect(1),
    limits = ((Lxₘᵢₙ, Lxₘₐₓ), (Lyₘᵢₙ, Lyₘₐₓ))
)
axis_kwargs_xz = (
    xlabel = "x",
    ylabel = "z",
    aspect = AxisAspect(1),
    limits = ((Lxₘᵢₙ, Lxₘₐₓ), (Lzₘᵢₙ, Lzₘₐₓ))
)

ax_xy = Axis(fig[1, 1]; axis_kwargs_xy...)
ax_xz = Axis(fig[1, 3]; axis_kwargs_xz...)

Tlims = (3, 20)

hm_xy = heatmap!(ax_xy, x, y, Txy; nan_color=:white, colormap = :thermal, colorrange = Tlims)
Colorbar(fig[1, 2], hm_xy; label = "m")

hm_xz = heatmap!(ax_xz, x, z, Txz; nan_color=:white, colormap = :thermal, colorrange = Tlims)
Colorbar(fig[1, 4], hm_xz; label = "m")

frames = intro:length(times)
record(fig, filename * ".mp4", frames, framerate=1) do i
    n[] = i
end
