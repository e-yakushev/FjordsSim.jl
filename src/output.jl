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
    ax = Axis(fig[1, 1])
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
     
    fig = Figure(size = (1000, 400))
    ax = Axis(fig[1, 1])
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
    ax = Axis(fig[1, 1])
    hm = heatmap!(ax, Ti, colorrange = colorrange, colormap = colormap)
    cb = Colorbar(fig[0, 1], hm, vertical = false)
    hidedecorations!(ax)

    CairoMakie.record(fig, joinpath(folder, "$(name).mp4"), 1:Nt, framerate = 8) do i
        iter[] = i
    end
end