using Oceananigans
using JLD2
using Printf
using Oceananigans.Units
using Oceananigans.Utils: prettytimeunits, maybe_int
using CairoMakie: Auto, Axis, Figure, GridLayout, Colorbar, Observable, Reverse, record, heatmap!, @lift
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
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(var_name)")

    record(fig, joinpath(folder, "$(var_name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

function prettiertime(t, longform=true)
    s = longform ? "seconds" : "s" 
    iszero(t) && return "0 $s"
    t < 1e-9 && return @sprintf("%.3e %s", t, s) # yah that's small

    t = maybe_int(t)
    value, units = prettytimeunits(t, longform)
    return @sprintf("%d %s", Int(trunc(Int, value)), units)
end

function record_variable_multilayer(
    variable,
    var_name,
    Nz_layers,
    times,
    folder;
    colorrange = (0, 0.5),
    colormap = :deep,
    framerate = 30,
)
    Nt = length(times)
    iter = Observable(Nt)
    num_layers = length(Nz_layers)
    figsize = (200 * num_layers, 400)  # Adaptive figure size based on number of layers

    figs = []
    titles = []
    heatmaps = []
    axes = []

    fig = Figure(size = figsize)
    grid = GridLayout(fig[1, 1])  # Stack vertically
    
    for (i, Nz) in enumerate(Nz_layers)
        f = @lift begin
            x = variable[$iter]
            x = interior(x, :, :, Nz)
            x[x .== 0] .= NaN
            x
        end

        title = @lift "Layer $Nz - " * prettiertime(times[$iter])
        push!(titles, title)

        ax = Axis(
            grid[1, i];
            title = title,
            titlealign = :left,
            width = Auto(),  # Adaptive width
            height = Auto()  # Adaptive height
        )
        push!(axes, ax)

        hm = heatmap!(ax, f, colorrange = colorrange, colormap = colormap)
        push!(heatmaps, hm)
    end

    cb = Colorbar(fig[1, 2], heatmaps[1], vertical = true, label = "$(var_name)")

    record(fig, joinpath(folder, "$(var_name)_multi.mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

folder = joinpath(homedir(), "FjordsSim_results", "oslofjord")
filename = joinpath(folder, "snapshots_ocean")
datafile = jldopen("$filename.jld2")
grid = datafile["grid"]
Nz = grid["underlying_grid"]["Nz"]

η = FieldTimeSeries("$filename.jld2", "free_surface", backend=OnDisk())
T = FieldTimeSeries("$filename.jld2", "T", backend=OnDisk())
S = FieldTimeSeries("$filename.jld2", "S", backend=OnDisk())
u = FieldTimeSeries("$filename.jld2", "u", backend=OnDisk())
u_ao_flux = FieldTimeSeries("$filename.jld2", "u_atm_ocean_flux", backend=OnDisk())
x_momentum = FieldTimeSeries("$filename.jld2", "x_momentum", backend=OnDisk())
v = FieldTimeSeries("$filename.jld2", "v", backend=OnDisk())
v_ao_flux = FieldTimeSeries("$filename.jld2", "v_atm_ocean_flux", backend=OnDisk())
y_momentum = FieldTimeSeries("$filename.jld2", "y_momentum", backend=OnDisk())
O₂ = FieldTimeSeries("$filename.jld2", "O₂")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
PHY = FieldTimeSeries("$filename.jld2", "P")
HET = FieldTimeSeries("$filename.jld2", "HET")
DOM = FieldTimeSeries("$filename.jld2", "DOM")
POM = FieldTimeSeries("$filename.jld2", "POM")
#P = FieldTimeSeries("$filename.jld2", "P")
#Z = FieldTimeSeries("$filename.jld2", "Z")
#NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
#NH₄ = FieldTimeSeries("$filename.jld2", "NH₄")
#DIC = FieldTimeSeries("$filename.jld2", "DIC")
#Alk = FieldTimeSeries("$filename.jld2", "Alk")

record_variable_multilayer(T, "temperature", [10, 12, 14, 16, 18], T.times, folder; colorrange = (0, 20), colormap = Reverse(:RdYlBu))
record_variable_multilayer(S, "salinity", [10, 12, 14, 16, 18], S.times, folder; colorrange = (20, 37), colormap = :viridis)
#record_variable_multilayer(DIC, "dissolved inorganic carbon", [10, 12, 14, 16, 18], DIC.times, folder; colorrange = (2200, 2300), colormap = Reverse(:CMRmap))
record_variable_multilayer(P, "phytoplankton", [10, 12, 14, 16, 18], P.times, folder; colorrange = (0, 1), colormap = Reverse(:cubehelix))

record_bottom_tracer(O₂, "oxygen", Nz, O₂.times, folder; figsize = (300, 700))

#record_variable(η, "free surface", 1, η.times, folder, (300, 450); colorrange = (0., 0.1))
record_variable(T, "temperature surface", Nz, T.times, folder, (300, 450); colorrange = (2, 6), colormap = Reverse(:RdYlBu))
record_variable(S, "salinity surface", Nz, S.times, folder, (300, 450); colorrange = (20, 37), colormap = :viridis)
record_variable(u, "u velocity surface", Nz, u.times, folder, (300, 450); colorrange = (-1, 1))
#record_variable(u_ao_flux, "u ao flux", 1, u_ao_flux.times, folder, (300, 450); colorrange = (-1, 1))
#record_variable(x_momentum, "x momentum", 1, x_momentum.times, folder, (300, 450); colorrange = (-1, 1))
record_variable(v, "v velocity surface", Nz, v.times, folder, (300, 450); colorrange = (-1, 1))
#record_variable(v_ao_flux, "v ao flux", 1, v_ao_flux.times, folder, (300, 450); colorrange = (-0.1, 0.1))
#record_variable(y_momentum, "y momentum", 1, y_momentum.times, folder, (300, 450); colorrange = (-0.1, 0.1))
record_variable(O₂, "O₂ surface", Nz, O₂.times, folder, (300, 700); colorrange = (0, 400), colormap = :turbo)
record_variable(NUT, "NUT surface", Nz, v.times, folder, (300, 700); colorrange = (0, 5))
record_variable(PHY, "PHY surface", Nz, v.times, folder, (300, 700); colorrange = (0, 5), colormap = Reverse(:cubehelix))
#record_variable(P, "P surface", Nz, P.times, folder, (300, 700); colorrange = (0, 1))
#record_variable(Z, "Z surface", Nz, Z.times, folder, (300, 700); colorrange = (0, 1))
#record_variable(NO₃, "NO₃ surface", Nz, NO₃.times, folder, (300, 700); colorrange = (0, 20))
#record_variable(DIC, "DIC surface", Nz, DIC.times, folder, (300, 700); colorrange = (2230, 2270))
#record_variable(Alk, "Alk surface", Nz, Alk.times, folder, (300, 700); colorrange = (2300, 2500))

#atmo_filename = joinpath(folder, "snapshots_atmosphere")
#atmo_datafile = jldopen("$atmo_filename.jld2")

#u_atm = FieldTimeSeries("$atmo_filename.jld2", "u_atm")
#v_atm = FieldTimeSeries("$atmo_filename.jld2", "v_atm")

#record_variable(u_atm, "u atm", 1, u_atm.times, folder, (300, 450); colorrange = (-20, 20))
#record_variable(v_atm, "v atm", 1, v_atm.times, folder, (300, 450); colorrange = (-20, 20))

