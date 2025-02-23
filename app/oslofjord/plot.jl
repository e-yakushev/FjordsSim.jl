using Oceananigans
using JLD2
using Oceananigans.Units
using CairoMakie: Axis, Figure, Colorbar, Observable, Reverse, record, heatmap!, @lift
using FjordsSim:
    record_bottom_tracer

function record_variable(
    variable,
    var_name,
    Nz,
    times,
    folder,
    figsize;
    colorrange = (0, 0.5),
    colormap = :deep,
    framerate = 12,
)
    Nt = length(times)
    iter = Observable(Nt)

    f = @lift begin
        x = variable[$iter]
        x = interior(x, :, :, Nz)
        x[x.==0] .= NaN
        x
    end

    fig = Figure(size = figsize)
    title = @lift "$(var_name) at " * prettytime(times[$iter])
    ax = Axis(
        fig[1, 1];
        title = title,
        xlabel = "Grid points, eastward direction",
        ylabel = "Grid points, northward direction",
    )
    hm = heatmap!(ax, f, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(var_name) (ms⁻¹)")

    record(fig, joinpath(folder, "$(var_name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

folder = joinpath(homedir(), "FjordsSim_results", "oslofjord")
filename = joinpath(folder, "snapshots")
T = FieldTimeSeries("$filename.jld2", "T")
S = FieldTimeSeries("$filename.jld2", "S")
u = FieldTimeSeries("$filename.jld2", "u")
v = FieldTimeSeries("$filename.jld2", "v")
O₂ = FieldTimeSeries("$filename.jld2", "O₂")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
PHY = FieldTimeSeries("$filename.jld2", "P")
HET = FieldTimeSeries("$filename.jld2", "HET")
DOM = FieldTimeSeries("$filename.jld2", "DOM")
POM = FieldTimeSeries("$filename.jld2", "POM")
P = FieldTimeSeries("$filename.jld2", "P")
Z = FieldTimeSeries("$filename.jld2", "Z")
NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
NH₄ = FieldTimeSeries("$filename.jld2", "NH₄")
DIC = FieldTimeSeries("$filename.jld2", "DIC")
Alk = FieldTimeSeries("$filename.jld2", "Alk")

grid = jldopen("$filename.jld2")["grid"]
Nz = grid["underlying_grid"]["Nz"]

record_bottom_tracer(O₂, "oxygen", Nz, O₂.times, folder; figsize = (300, 700))

record_variable(T, "temperature surface", Nz, T.times, folder, (300, 450); colorrange = (2, 6))
record_variable(S, "salinity surface", Nz, S.times, folder, (300, 450); colorrange = (20, 37))
record_variable(u, "u velocity surface", Nz, u.times, folder, (300, 450); colorrange = (-1, 1))
record_variable(v, "v velocity surface", Nz, v.times, folder, (300, 450); colorrange = (-1, 1))
record_variable(O₂, "O₂ surface", Nz, O₂.times, folder, (300, 700); colorrange = (0, 400))
record_variable(NUT, "NUT surface", Nz, v.times, folder, (300, 700); colorrange = (0, 5))
record_variable(PHY, "PHY surface", Nz, v.times, folder, (300, 700); colorrange = (0, 5))
record_variable(P, "P surface", Nz, P.times, folder, (300, 700); colorrange = (0, 1))
record_variable(Z, "Z surface", Nz, Z.times, folder, (300, 700); colorrange = (0, 1))
record_variable(NO₃, "NO₃ surface", Nz, NO₃.times, folder, (300, 700); colorrange = (0, 20))
record_variable(DIC, "DIC surface", Nz, DIC.times, folder, (300, 700); colorrange = (2230, 2270))
record_variable(Alk, "Alk surface", Nz, Alk.times, folder, (300, 700); colorrange = (2300, 2500))
