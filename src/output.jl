using CairoMakie: Axis, Figure, Colorbar, Observable, Reverse, record, heatmap!, @lift, lines!, axislegend
using FileIO: save
using Oceananigans.Fields: Field, interior, compute!

function plot_1d_phys(T, S, z, times, folder, x, y)
    fig = Figure(size = (1000, 400), fontsize = 20)

    axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", xticks = (0:30:times[end]/days), xtickformat = "{:.0f}")

    Axis(fig[1, 1]; title = "T, ⁰C", axis_kwargs...)
    hmT = heatmap!(times / days, z, interior(T, x, y, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[1, 2], hmT)

    Axis(fig[2, 1]; title = "S, psu", axis_kwargs...)
    hmS = heatmap!(times / days, z, interior(S, x, y, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 2], hmS)

    save(joinpath(folder, "1d_phys.png"), fig)
end

map_axis_kwargs = (xlabel = "Grid points, eastward direction", ylabel = "Grid points, northward direction")
transect_axis_kwargs = (xlabel = "Grid points, eastward direction", ylabel = "z (m)")
framerate = 12

function record_surface_speed(u, v, Nz, times, folder; colorrange = (0, 0.5), colormap = :deep)
    Nt = length(times)
    iter = Observable(Nt)

    ## Speed
    si = @lift begin
        s = Field(sqrt(u[$iter]^2 + v[$iter]^2))
        compute!(s)
        s = interior(s, :, :, Nz)
        s[s.==0] .= NaN
        s
    end

    fig = Figure(size = (1000, 400))

    title = @lift "Surface speed at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, map_axis_kwargs...)
    hm = heatmap!(ax, si, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Surface speed (ms⁻¹)")
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "surface_speed.mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

function record_bottom_tracer(
    variable,
    var_name,
    Nz,
    times,
    folder;
    colorrange = (-1, 350),
    colormap = :turbo,
    figsize = (1000, 400),
)

    # bottom_z evaluation
    bottom_z = ones(Int, size(variable, 1), size(variable, 2))
    for i = 1:size(variable, 1)
        for j = 1:size(variable, 2)
            for k = size(variable, 3):-1:1  # Loop backwards to find the latest non-zero
                if variable[i, j, k, 1] == 0
                    bottom_z[i, j] = k
                    if k != Nz
                        bottom_z[i, j] = k + 1
                    end
                    break
                end
            end
        end
    end

    iter = Observable(1)
    f = @lift begin
        x = [variable[i, j, bottom_z[i, j], $iter] for i = 1:size(variable, 1), j = 1:size(variable, 2)]
        x[x.==0] .= NaN
        x
    end
    title = @lift "bottom $(var_name), mmol/m³ at " * prettytime(times[$iter])
    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1]; title = title, map_axis_kwargs...)
    hm = heatmap!(ax, f, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "$(var_name), mmol/m³")

    Nt = length(times)
    record(fig, joinpath(folder, "bottom_$(var_name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

function record_horizontal_tracer(tracer, times, folder, name, label; colorrange = (-1, 30), colormap = :magma, iz = 10)
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        Ti = interior(tracer[$iter], :, :, iz)
        Ti[Ti.==0] .= NaN
        Ti
    end

    title = @lift label * " at " * prettytime(times[$iter])
    fig = Figure(size = (1000, 400))
    ax = Axis(fig[1, 1]; title = title, map_axis_kwargs...)
    hm = heatmap!(ax, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

function record_vertical_tracer(tracer, depth, iy, times, folder, name, label; colorrange = (-1, 30), colormap = :magma)

    xs = 1:size(tracer)[1] # get x-values for x-axis
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        Ti = interior(tracer[$iter], :, iy, :)
        Ti[Ti.==0] .= NaN
        Ti
    end

    fig = Figure(size = (1000, 400))

    title = @lift label * " at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, transect_axis_kwargs...)
    hm = heatmap!(ax, xs, depth, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

function record_vertical_tracer_points(
    tracer,
    depth,
    indices::Vector{Tuple{Int, Int}},  # List of (ix, iy) pairs
    times,
    folder,
    name,
    label;
    colorrange = (-1, 30),
    colormap = :magma,
)

    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        # Collect slices from all given (ix, iy) points
        slices = [interior(tracer[$iter], ix, iy, :) for (ix, iy) in indices]
        Ti = hcat(slices...)  # Stack slices
        Ti[Ti.==0] .= NaN
        Ti = Ti'
        Ti
    end

    xs = 1:size(indices)[1] # get x-values for x-axis
    fig = Figure(size = (1000, 400))

    title = @lift label * " at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, transect_axis_kwargs...)
    hm = heatmap!(ax, xs, depth, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

function record_vertical_diff(
    tracer,
    depth,
    iy,
    times,
    folder,
    name,
    label;
    colorrange = (-1, 30),
    colormap = :magma,
)

    xs = 1:size(tracer)[1] # get x-values for x-axis
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
        Ti = tracer[:, iy, :, $iter]
        Ti[Ti.==0] .= NaN
        Ti
    end

    fig = Figure(size = (1000, 400))

    title = @lift label * " at " * prettytime(times[$iter])
    ax = Axis(fig[1, 1]; title = title, transect_axis_kwargs...)
    hm = heatmap!(ax, xs, depth, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = label)
    # hidedecorations!(ax)

    record(fig, joinpath(folder, "$(name).mp4"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

function plot_ztime(PHY, HET, POM, DOM, NUT, O₂, T, S, i, j, times, z, folder)

    fig = Figure(size = (1500, 1000), fontsize = 20)

    axis_kwargs = (
        xlabel = "Time (days)",
        ylabel = "z (m)",
        xticks = (0:50:times[end]),
        xtickformat = "{:.0f}", #   values -> ["$(value)kg" for value in values]     
    )

    axPHY = Axis(fig[1, 3]; title = "PHY, mmolN/m³", axis_kwargs...)
    hmPHY = heatmap!(times / days, z, interior(PHY, i, j, :, :)', colormap = Reverse(:cubehelix)) #(:davos10))
    Colorbar(fig[1, 4], hmPHY)

    axHET = Axis(fig[2, 3]; title = "HET, mmolN/m³", axis_kwargs...)
    hmHET = heatmap!(times / days, z, interior(HET, i, j, :, :)', colormap = Reverse(:afmhot))
    Colorbar(fig[2, 4], hmHET)

    axPOM = Axis(fig[3, 3]; title = "POM, mmolN/m³", axis_kwargs...)
    hmPOM = heatmap!(times / days, z, interior(POM, i, j, :, :)', colormap = Reverse(:greenbrownterrain)) #(:bilbao25))
    hmPOM = heatmap!(times / days, z, interior(POM, i, j, :, :)', colormap = Reverse(:greenbrownterrain)) #(:bilbao25))
    Colorbar(fig[3, 4], hmPOM)

    axDOM = Axis(fig[3, 1]; title = "DOM, mmolN/m³", axis_kwargs...)
    hmDOM = heatmap!(times / days, z, interior(DOM, i, j, :, :)', colormap = Reverse(:CMRmap)) #(:devon10))
    Colorbar(fig[3, 2], hmDOM)

    axNUT = Axis(fig[1, 1]; title = "NUT, mmolN/m³", axis_kwargs...)
    hmNUT = heatmap!(times / days, z, interior(NUT, i, j, :, :)', colormap = Reverse(:cherry))
    hmNUT = heatmap!(times / days, z, interior(NUT, i, j, :, :)', colormap = Reverse(:cherry))
    Colorbar(fig[1, 2], hmNUT)

    axOXY = Axis(fig[2, 1]; title = "OXY, mmol/m³", axis_kwargs...)
    hmOXY = heatmap!(times / days, z, interior(O₂, i, j, :, :)', colormap = :turbo)
    hmOXY = heatmap!(times / days, z, interior(O₂, i, j, :, :)', colormap = :turbo)
    Colorbar(fig[2, 2], hmOXY)

    axT = Axis(fig[2, 5]; title = "T, oC", axis_kwargs...)
    hmT = heatmap!(times / days, z, interior(T, i, j, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 6], hmT)

    axS = Axis(fig[3, 5]; title = "S, psu", axis_kwargs...)
    hmS = heatmap!(times / days, z, interior(S, i, j, :, :)', colormap = :viridis)
    Colorbar(fig[3, 6], hmS)

    @info "VARIABLES Z-Time plots made"

    save(joinpath(folder, "ztime.png"), fig)
end

function plot_bottom_tracer(tracer, bottom_z, time, folder)

    bottom_tracer = [tracer[i, j, bottom_z[i, j], time] for i = 1:size(tracer, 1), j = 1:size(tracer, 2)]
    fig = Figure(size = (1000, 400), fontsize = 20)

    axis_kwargs = (xlabel = "Grid points, eastward direction", ylabel = "Grid points, northward direction")

    axOXY = Axis(fig[2, 1]; title = "$(tracer), mmol/m³, " * prettytime(time), axis_kwargs...)
    hmOXY = heatmap!([i for i = 1:size(tracer, 1)], [j for j = 1:size(tracer, 2)], bottom_tracer, colormap = :turbo)
    Colorbar(fig[2, 2], hmOXY)

    save(joinpath(folder, "bottom_$(tracer).png"), fig)
end

###########################

function plot_ratio_under_thresh(O₂, Nz, times, folder)

    # bottom_z evaluation
    bottom_z = ones(Int, size(O₂, 1), size(O₂, 2))
    for i = 1:size(O₂, 1)
        for j = 1:size(O₂, 2)
            for k = 1:size(O₂, 3)
                if O₂[i, j, k, 1] .!= 0
                    bottom_z[i, j] = k
                    if k == 20
                        bottom_z[i, j] = 20
                    end
                    break
                end
            end
        end
    end

    # bottom O₂ Matrix evaluation
    O₂_bot = Array{eltype(O₂)}(undef, size(O₂, 1), size(O₂, 2), size(O₂, 4))
    for i in 1:size(O₂, 1)
        for j in 1:size(O₂, 2)
            O₂_bot[i, j, :] = O₂[i, j, bottom_z[i, j], :]
        end
    end

    land_cells = length(O₂_bot[:,:,1][O₂_bot[:,:,1] .< 10])
    water_cells = length(bottom_z) - land_cells
    
    under_thresh(thresh) = [length(O₂_bot[:,:,i][O₂_bot[:,:,i] .< thresh]) - land_cells for i in 1:size(O₂_bot, 3)]
    
    fig = Figure(size = (1000, 400), fontsize = 20)
    
    axis_kwargs = (xlabel = "Time (days)", ylabel = "% of bottom grid cells", xticks = (0:30:times[end]/days), xtickformat = "{:.0f}")
    
    Axis(fig[1, 1]; axis_kwargs...)
    
    lines!(times / days, 100*under_thresh(90)/water_cells, colormap = Reverse(:RdYlBu), label="> 90 μM")
    lines!(times / days, 100*under_thresh(15)/water_cells, colormap = Reverse(:RdYlBu), label="> 15 μM")
    lines!(times / days, 100*under_thresh(1)/water_cells, colormap = Reverse(:RdYlBu), label="> 1 μM")
    axislegend()

    save(joinpath(folder, "ratio_uder_thresh.png"), fig)
end