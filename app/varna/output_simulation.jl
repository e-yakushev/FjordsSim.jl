using Oceananigans
using CairoMakie

prefix = "surface"
filepath = joinpath(homedir(), "data_Varna", "$(prefix)_snapshots.jld2")

u = FieldTimeSeries(filepath, "u"; backend = OnDisk())
v = FieldTimeSeries(filepath, "v"; backend = OnDisk())
T = FieldTimeSeries(filepath, "T"; backend = OnDisk())
S = FieldTimeSeries(filepath, "S"; backend = OnDisk())
P = FieldTimeSeries(filepath, "P"; backend = OnDisk())
NO₃ = FieldTimeSeries(filepath, "NO₃"; backend = OnDisk()) 

times = u.times
Nt = length(times)

iter = Observable(Nt)

## Speed
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

## Temperature
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

## Salinity
Si = @lift begin
     Si = interior(S[$iter], :, :, 1)
     Si[Si .== 0] .= NaN
     Si
end
 
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Si, colorrange = (0, 15), colormap = :magma)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(prefix) Salinity (psu)")
hidedecorations!(ax)

CairoMakie.record(fig, joinpath(homedir(), "data_Varna", "$(prefix)_S.mp4"), 1:Nt, framerate = 8) do i
    iter[] = i
end

## Phy
Pi = @lift begin
     Pi = interior(P[$iter], :, :, 1)
     Pi[Pi .== 0] .= NaN
     Pi
end
 
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Pi, colorrange = (0, 15), colormap = :magma)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(prefix) Phytoplankton")
hidedecorations!(ax)

CairoMakie.record(fig, joinpath(homedir(), "data_Varna", "$(prefix)_Phy.mp4"), 1:Nt, framerate = 8) do i
    iter[] = i
end

## NO₃
Ni = @lift begin
     Ni = interior(NO₃[$iter], :, :, 1)
     Ni[Ni .== 0] .= NaN
     Ni
end
 
fig = Figure(size = (800, 400))
ax = Axis(fig[1, 1])
hm = heatmap!(ax, Ni, colorrange = (0, 10), colormap = :magma)
cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(prefix) NO₃ (ml/l)")
hidedecorations!(ax)

CairoMakie.record(fig, joinpath(homedir(), "data_Varna", "$(prefix)_NO3.mp4"), 1:Nt, framerate = 8) do i
    iter[] = i
end

## Profile temperature
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
