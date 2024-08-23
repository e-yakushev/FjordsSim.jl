using Oceananigans
using CairoMakie

prefix = "surface"
filepath = joinpath(homedir(), "data_Varna", "$(prefix)_snapshots.jld2")

u = FieldTimeSeries(filepath, "u"; backend = OnDisk())
v = FieldTimeSeries(filepath, "v"; backend = OnDisk())
T = FieldTimeSeries(filepath, "T"; backend = OnDisk())

times = u.times
Nt = length(times)

iter = Observable(Nt)

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
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(prefix) speed (ms⁻¹)")
hidedecorations!(ax)

CairoMakie.record(fig, joinpath(homedir(), "data_Varna", "$(prefix)_s.mp4"), 1:Nt, framerate = 8) do i
    iter[] = i
end

Ti = @lift begin
     Ti = interior(T[$iter], :, :, 1)
     Ti[Ti .== 0] .= NaN
     Ti
end
 
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Ti, colorrange = (-1, 30), colormap = :magma)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(prefix) Temperature (Cᵒ)")
hidedecorations!(ax)

CairoMakie.record(fig, joinpath(homedir(), "data_Varna", "$(prefix)_T.mp4"), 1:Nt, framerate = 8) do i
    iter[] = i
end

prefix = "profile"
filepath = joinpath(homedir(), "data_Varna", "$(prefix)_snapshots.jld2")

u = FieldTimeSeries(filepath, "u"; backend = OnDisk())
v = FieldTimeSeries(filepath, "v"; backend = OnDisk())
T = FieldTimeSeries(filepath, "T"; backend = OnDisk())

times = u.times
Nt = length(times)

iter = Observable(Nt)

Ti = @lift begin
     Ti = interior(T[$iter], :, 1, :)
     Ti[Ti .== 0] .= NaN
     Ti
end
 
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Ti, colorrange = (-1, 30), colormap = :magma)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(prefix) Temperature (Cᵒ)")
hidedecorations!(ax)

CairoMakie.record(fig, joinpath(homedir(), "data_Varna", "$(prefix)_T.mp4"), 1:Nt, framerate = 8) do i
    iter[] = i
end
