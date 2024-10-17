using Oceananigans
using CairoMakie
using Oceananigans.Units

function plot_1d_phys(T, S, z, times, folder)
    fig = Figure(size = (1000, 400), fontsize = 20)

    axis_kwargs = (
        xlabel = "Time (days)",
        ylabel = "z (m)",
        xticks = (0:30:times[end]/days),
        xtickformat = "{:.0f}",
    )

    Axis(fig[1, 1]; title = "T, ⁰C", axis_kwargs...)
    hmT = heatmap!(times / days, z, interior(T, 1, 1, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[1, 2], hmT)

    Axis(fig[2, 1]; title = "S, psu", axis_kwargs...)
    hmS = heatmap!(times / days, z, interior(S, 1, 1, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 2], hmS)

    save(joinpath(folder,"1d_phys.png"), fig)
end

function record_surface_speed(
    u, v, Nz, times, folder;
    colorrange = (0, 0.5), colormap = :deep,
    )
    Nt = length(times)
    iter = Observable(Nt)

    ## Speed
    si = @lift begin
         s = Field(sqrt(u[$iter]^2 + v[$iter]^2))
         compute!(s)
         s = interior(s, :, :, Nz)
         s[s .== 0] .= NaN
         s
    end

    fig = Figure(size = (1000, 400))

    title = @lift "Surface speed at " * prettytime(times[$iter]) 
    ax = Axis(fig[1, 1], title = title)
    hm = heatmap!(ax, si, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false, label = "Surface speed (ms⁻¹)")
    hidedecorations!(ax)

    CairoMakie.record(fig, joinpath(folder, "surface_speed.mp4"), 1:Nt, framerate = 8) do i
        iter[] = i
    end
end

function record_surface_tracer(
    tracer, Nz, times, folder, name, label;
    colorrange=(-1, 30), colormap=:magma,
    )
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
         Ti = interior(tracer[$iter], :, :, Nz)
         Ti[Ti .== 0] .= NaN
         Ti
    end

    title = @lift label * " at " * prettytime(times[$iter]) 
    fig = Figure(size = (1000, 400))
    ax = Axis(fig[1, 1], title=title)
    hm = heatmap!(ax, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false)
    hidedecorations!(ax)

    CairoMakie.record(fig, joinpath(folder, "$(name).mp4"), 1:Nt, framerate = 8) do i
        iter[] = i
    end
end

function record_vertical_tracer(
    tracer, iy, times, folder, name, label;
    colorrange=(-1, 30), colormap=:magma,
    )
    Nt = length(times)
    iter = Observable(Nt)

    Ti = @lift begin
         Ti = interior(tracer[$iter], :, iy, :)
         Ti[Ti .== 0] .= NaN
         Ti
    end
     
    fig = Figure(size = (1000, 400))

    title = @lift label * " at " * prettytime(times[$iter]) 
    ax = Axis(fig[1, 1], title = title)
    hm = heatmap!(ax, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false)
    hidedecorations!(ax)

    CairoMakie.record(fig, joinpath(folder, "$(name).mp4"), 1:Nt, framerate = 8) do i
        iter[] = i
    end
end



function plot_ztime(
    PHY, HET, POM, DOM, NUT, O₂, T, S, i, j, times, z, folder,
    )

    fig = Figure(size = (1500, 1000), fontsize = 20)

    axis_kwargs = (
        xlabel = "Time (days)",
        ylabel = "z (m)",
        xticks = (0:30:times[end]),
        xtickformat = "{:.0f}" #   values -> ["$(value)kg" for value in values]     
    )

    axPHY = Axis(fig[1, 3]; title = "PHY, mmolN/m³", axis_kwargs...)
    hmPHY = heatmap!(times / days, z, interior(PHY, i, j, :, :)', colormap = Reverse(:cubehelix)) #(:davos10))
    Colorbar(fig[1, 4], hmPHY)

    axHET = Axis(fig[2, 3]; title = "HET, mmolN/m³", axis_kwargs...)
    hmHET = heatmap!(times / days, z, interior(HET, i, j, :, :)', colormap = Reverse(:afmhot))
    Colorbar(fig[2, 4], hmHET)

    axPOM = Axis(fig[3, 3]; title = "POM, mmolN/m³", axis_kwargs...)
    hmPOM =
        heatmap!(times / days, z, interior(POM, i, j, :, :)', colormap = Reverse(:greenbrownterrain)) #(:bilbao25))
    hmPOM =
        heatmap!(times / days, z, interior(POM, i, j, :, :)', colormap = Reverse(:greenbrownterrain)) #(:bilbao25))
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

    # axκ = Axis(fig[1, 5]; title = "κ  m³/s", axis_kwargs...)
    # hmκ = heatmap!(times / days, z, interior(κ, 1, 1, :, :)', colormap = Reverse(:RdYlBu)) # :linear_grey_0_100_c0_n256)
    # Colorbar(fig[1, 6], hmκ)

    axT = Axis(fig[2, 5]; title = "T, oC", axis_kwargs...)
    hmT = heatmap!(times / days, z, interior(T, i, j,  :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 6], hmT)

    axS = Axis(fig[3, 5]; title = "S, psu", axis_kwargs...)
    hmS = heatmap!(times / days, z, interior(S, i, j, :, :)', colormap = :viridis)
    Colorbar(fig[3, 6], hmS)

    # axPAR = Axis(fig[4, 1]; title = "PAR  μE⋅m-2⋅s-1", axis_kwargs...)
    # hmPAR = heatmap!(times / days, z, interior(PAR, 1, 1, :, :)', colormap = :grayC100) # :linear_grey_0_100_c0_n256)
    # Colorbar(fig[4, 2], hmPAR)

    @info "VARIABLES Z-Time plots made"

    save(joinpath(folder, "ztime.png"), fig)

end