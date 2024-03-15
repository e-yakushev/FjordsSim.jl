using Oceananigans
using GLMakie

Nx = 50
Ny = 50
Nz = 10;

grid = LatitudeLongitudeGrid(CPU();
                             size = (Nx, Ny, Nz),
                             latitude = (58.8, 59.9),
                             longitude = (10.1, 11.1),
                             z = (-500, 0),
                             halo = (4, 4, 4))

filepath = joinpath(homedir(), "src", "FjOs.jl", "oslo_fjord_50_50_10_fine_snapshots.jld2")

time_series = (w = FieldTimeSeries(filepath, "w"; grid),
               T = FieldTimeSeries(filepath, "T"; grid))
xw, yw, zw = nodes(time_series.w)
xT, yT, zT = nodes(time_series.T)

snap = interior(time_series.w, :, :, 10, 3)

fig = Figure(size = (500, 400))

axis_kwargs = (xlabel="x",
               ylabel="y",
               aspect = AxisAspect(grid.Lx/grid.Lz),
               limits = ((0, grid.Lx), (-grid.Lz, 0)))

ax_w  = Axis(fig[1, 1]; title = "Vertical velocity", axis_kwargs...)
ax_T  = Axis(fig[2, 1]; title = "Temperature", axis_kwargs...)

level = 10
time_snap = 1
wlims = (-0.05, 0.05)
Tlims = (18.0, 22.00)

wₙ = interior(time_series.w, :, :, level, time_snap)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colormap = :balance, colorrange = wlims)
Colorbar(fig[1, 2], hm_w; label = "m s⁻¹")

Tₙ = interior(time_series.T, :, :, level, time_snap)
hm_T = heatmap!(ax_T, xT, zT, Tₙ; colormap = :thermal, colorrange = Tlims)
Colorbar(fig[2, 2], hm_T; label = "ᵒC")

display(fig)
