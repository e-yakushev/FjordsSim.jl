using Oceananigans
using CairoMakie

u = FieldTimeSeries("surface.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("surface.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("surface.jld2", "T"; backend = OnDisk())
e = FieldTimeSeries("surface.jld2", "e"; backend = OnDisk())

times = u.times
Nt = length(times)

iter = Observable(Nt)

Ti = @lift begin
     Ti = interior(T[$iter], :, :, 1)
     Ti[Ti .== 0] .= NaN
     Ti
end

ei = @lift begin
     ei = interior(e[$iter], :, :, 1)
     ei[ei .== 0] .= NaN
     ei
end

si = @lift begin
     s = Field(sqrt(u[$iter]^2 + v[$iter]^2))
     compute!(s)
     s = interior(s, :, :, 1)
     s[s .== 0] .= NaN
     s
end


fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, si, colorrange = (0, 0.5), colormap = :deep)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Surface speed (ms⁻¹)")
hidedecorations!(ax)

CairoMakie.record(fig, "near_global_ocean_surface_s.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
end
nothing #hide
 
 # ![](near_global_ocean_surface_s.mp4)
 
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Ti, colorrange = (-1, 30), colormap = :magma)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Surface Temperature (Cᵒ)")
hidedecorations!(ax)

CairoMakie.record(fig, "near_global_ocean_surface_T.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
end
nothing #hide
 
# ![](near_global_ocean_surface_T.mp4)

fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, ei, colorrange = (0, 1e-3), colormap = :solar)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Turbulent Kinetic Energy (m²s⁻²)")
hidedecorations!(ax)

CairoMakie.record(fig, "near_global_ocean_surface_e.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
end
nothing #hide
