using Oceananigans
using JLD2
using NCDatasets
using NetCDF
using Printf
using Oceananigans.Units
using Oceananigans.Utils: prettytimeunits, maybe_int
using CairoMakie: 
      Auto, Axis, Figure, GridLayout, Colorbar, 
      rowgap!, colgap!,GridLayout, Relative,
      Observable, Reverse, record, heatmap!, @lift
using FjordsSim:
    plot_1d_phys,
    extract_z_faces,
    record_vertical_tracer,
    record_surface_speed,
#    plot_ztime,   
#    record_bottom_tracer,
    BGCModel,
    oxygen_saturation  
include("/home/eya/src/FjordsSim.jl/src/BGCModels/BGCModels.jl")
using .BGCModels
include("/home/eya/src/FjordsSim.jl/src/BGCModels/boundary_conditions.jl")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

map_axis_kwargs = (xlabel = "Grid points, East", ylabel = "Grid points, North")
transect_axis_kwargs = (xlabel = "Grid points, East", ylabel = "z (m)")
framerate = 12

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make animated gif of changes at selected depth
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make animated gif of changes at the bottom
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function record_bottom_tracer(variable, var_name, Nz, times, folder;
    colorrange = (-1, 350), colormap = :turbo, figsize = (1000, 400),)

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
    record(fig, joinpath(folder, "movie_$(var_name).gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make animated gif of changes at selected depth
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function record_horizontal_tracer(tracer, times, folder, name, label;
                                  colorrange = (-1, 30), colormap = :magma, iz = 10)

    Nt = length(times)
    iter = Observable(1)
    Ti = @lift begin
        if tracer isa AbstractArray
            Ti = tracer[:, :, iz, $iter]
        elseif tracer isa FieldTimeSeries
            Ti = interior(tracer[$iter], :, :, iz)
        else
            error("Unsupported tracer type: $(typeof(tracer))")
        end
        Ti[Ti .== 0] .= NaN
        Ti
    end

    title = @lift "$label at $(prettytime(times[$iter]))"
    fig = Figure(size = (400, 550), fontsize = 20)
    ax = Axis(fig[1, 1]; title = title, map_axis_kwargs...)
    hm = heatmap!(ax, Ti, colorrange = colorrange, colormap = colormap, nan_color = :silver)
    cb = Colorbar(fig[0, 1], hm, vertical = false)

    record(fig, joinpath(folder, "movie_$(name)_iz_$iz.gif"), 1:Nt, framerate = framerate) do i
        iter[] = i
    end

    @info "movie_$(name)_iz_$iz record made"
end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Do something with time (Shamil' function  )
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot tracer distribution at given day (plot_day) at given depth (iz)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_surface_tracer(tracer, name, iz, times, plot_day, folder; 
                             colorrange = (0, 0.5), colormap = :turbo, )
# calculate the index in the time array corresponding to the plot_day                             
   nday = plot_day*round(Int, (length(times)/365))
# extract tracer data at iz at given day   
   z_tracer = [tracer[i, j, iz, nday] == 0 ? NaN : tracer[i, j, iz, nday] for i in 1:size(tracer, 1), j in 1:size(tracer, 2)]
# create the plot
   fig = Figure(size = (400, 550), fontsize = 20)
    axis_kwargs = (xlabel = "Grid points, East ", ylabel = "Grid points, North")
    axTRAC = Axis(fig[1, 1]; title = "$(name), μM, day $(plot_day) ", axis_kwargs...)
    hmTRAS = heatmap!([i for i = 1:size(tracer, 1)], [j for j = 1:size(tracer, 2)], z_tracer, colorrange = colorrange, colormap = colormap, nan_color = :grey) #:gray)
    Colorbar(fig[1, 2], hmTRAS)
    save(joinpath(folder, "surface_$(name)_day_$(plot_day)_iz_$(iz).png"), fig)
    @info "surface_$(name)_day_$(plot_day) plot made"    
end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot tracer distribution at given day (plot_day) at the bottom  (bottom_z)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_bottom_tracer!(tracer, name, bottom_z, times, plot_day, folder; colorrange = (0, 0.5), colormap = :turbo)
# calculate the index in the time array corresponding to the plot_day      
    nday = plot_day*round(Int, (length(times)/365))
# extract tracer data at iz at given day  
    z_tracer = [tracer[i, j, bottom_z[i, j], nday] == 0 ? NaN : tracer[i, j, bottom_z[i, j], nday] for i in 1:size(tracer, 1), j in 1:size(tracer, 2)]
# create the plot
    fig = Figure(size = (400, 550), fontsize = 20)
    axis_kwargs = (xlabel = "Grid points, East ", ylabel = "Grid points, North")
    axTRAC = Axis(fig[1, 1]; title = "$(name), μM, day $(plot_day) ", axis_kwargs...) 
    hmTRAS = heatmap!([i for i = 1:size(tracer, 1)], [j for j = 1:size(tracer, 2)], z_tracer, colorrange = colorrange, colormap = colormap, nan_color = :silver)
    Colorbar(fig[1, 2], hmTRAS)
    save(joinpath(folder, "$(name)_bottom_day_$(plot_day).png"), fig)  
    @info "$(name)_bottom_day_$(plot_day) plot made"    
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# subplot function for tracer plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_tracer_subplot!(fig, pos, 
                        data, title_str; colorrange=(0,1), colormap=:viridis)
    ax = Axis(fig[pos...]; title=title_str, 
            width = 200, #Auto(),  # Adaptive width
            height = 200, #Auto()  # Adaptive height
            xlabel="", ylabel="")
    hm = heatmap!(ax, data; colorrange=colorrange, colormap=colormap, nan_color=:silver)
    Colorbar(fig[pos[1], pos[2]+1], hm, vertical=true)
end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# replace zeros with NaN in 4D array slice
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function replace_zeros_with_NaN!(A, depth_index, day_index)
    slice = Float64.(view(A, :, :, depth_index, day_index))
    @. slice = ifelse(slice == 0, NaN, slice)
    return slice
end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# vertical distribution changes at a point (i,j)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Helper function to safely get the "interior" data
get_interior(A, inds...) = 
    hasmethod(interior, Tuple{typeof(A), Vararg{Any}}) ? interior(A, inds...) : A[inds...]
# Seasonal changes at a point (i,j)
function plot_ztime(NUT, O₂, O₂_relative, PHY, HET, T, DOM, POM, S, i, j, times, z, folder)
    fig = Figure(size = (1500, 1000), fontsize = 20)

    axis_kwargs = (
        xlabel = "Time (days)",
        ylabel = "z (m)",
        xticks = (0:50:times[end]),
        xtickformat = "{:.0f}",
    )

    axNUT = Axis(fig[1, 1]; title = "NUT, mmolN/m³", axis_kwargs...)
    hmNUT = heatmap!(times / days, z, get_interior(NUT, i, j, :, :)', colormap = Reverse(:cherry))
    Colorbar(fig[1, 2], hmNUT)
    
    axOXY = Axis(fig[1, 3]; title = "O₂, mmol/m³", axis_kwargs...)
    hmOXY = heatmap!(times / days, z, get_interior(O₂, i, j, :, :)', colormap = :turbo)
    Colorbar(fig[1, 4], hmOXY)
    
    axOXY_rel = Axis(fig[1, 5]; title = "O₂ saturation, %", axis_kwargs...)
    hmOXY_rel = heatmap!(times / days, z, get_interior(O₂_relative, i, j, :, :)', colormap = :gist_stern)
    Colorbar(fig[1, 6], hmOXY_rel)

    axPHY = Axis(fig[2, 1]; title = "PHY, mmolN/m³", axis_kwargs...)
    hmPHY = heatmap!(times / days, z, get_interior(PHY, i, j, :, :)', colormap = Reverse(:cubehelix))
    Colorbar(fig[2, 2], hmPHY)

    axHET = Axis(fig[2, 3]; title = "HET, mmolN/m³", axis_kwargs...)
    hmHET = heatmap!(times / days, z, get_interior(HET, i, j, :, :)', colormap = Reverse(:afmhot))
    Colorbar(fig[2, 4], hmHET)

    axT = Axis(fig[2, 5]; title = "T, °C", axis_kwargs...)
    hmT = heatmap!(times / days, z, get_interior(T, i, j, :, :)', colormap = Reverse(:RdYlBu))
    Colorbar(fig[2, 6], hmT)

    axDOM = Axis(fig[3, 1]; title = "DOM, mmolN/m³", axis_kwargs...)
    hmDOM = heatmap!(times / days, z, get_interior(DOM, i, j, :, :)', colormap = Reverse(:CMRmap))
    Colorbar(fig[3, 2], hmDOM)

    axPOM = Axis(fig[3, 3]; title = "POM, mmolN/m³", axis_kwargs...)
    hmPOM = heatmap!(times / days, z, get_interior(POM, i, j, :, :)', colormap = Reverse(:greenbrownterrain))
    Colorbar(fig[3, 4], hmPOM)

    axS = Axis(fig[3, 5]; title = "S, psu", axis_kwargs...)
    hmS = heatmap!(times / days, z, get_interior(S, i, j, :, :)', colormap = :viridis)
    Colorbar(fig[3, 6], hmS)

    save(joinpath(folder, "ztime_$(i)_$(j).png"), fig)
    @info "Saved ztime_$(i)_$(j) plot in $folder"
end
# ====================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN CODE STARTS HERE: open file and extract data to "ds"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder = joinpath(homedir(), "FjordsSim_results", "oslofjord")
filename = joinpath(folder, "snapshots_ocean")
ds = NCDataset("$filename.nc", "r")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
println("1D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 1
            println("$varname, Dimensions: ", dimnames(var), ", Sizes: ", size(var))
        end
    end
println("2D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 2
            println("$varname, Dimensions: ", dimnames(var),", Sizes: ", size(var))
        end
    end
println("3D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 3
            println("$varname, Dimensions: ", dimnames(var),", Sizes: ", size(var))
        end
    end
println("4D Variables in $filename.nc:")
    for (varname, var) in ds
        if ndims(var) == 4
            println("$varname, Dimensions: ", dimnames(var), ", Sizes: ", size(var))
        end
    end
        
# Get grid dimensions and properties
grid_group = ds.group["grid_reconstruction"]
Nx = grid_group.attrib["Nx"]
Ny = grid_group.attrib["Ny"]
Nz = grid_group.attrib["Nz"]

println("Grid dimensions: Nx=$Nx,", typeof(Nx),", Ny=$Ny,", typeof(Ny),", Nz=$Nz,", typeof(Nz))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract from NetCDF dataset 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
times = ds["time"][:]         # time in seconds, float
depth = ds["z_aac"][:]        # depth in m, starting from the deepest, negative, float
println("Time stats (in seconds) — min: ", minimum(times), ", max: ", maximum(times))
# --- Extract  T from NetCDF dataset --- 
T = ds["T"][:,:,:,:]               
println("T stats — min: ", minimum(T), ", max: ", maximum(T))
# --- Extract  S from NetCDF dataset --- 
S = ds["S"][:,:,:,:]                
println("S stats — min: ", minimum(S), ", max: ", maximum(S))
# --- Extract  P from NetCDF dataset --- 
P = ds["P"][:,:,:,:]               
println("P stats — min: ", minimum(P),  ", max: ", maximum(P))
# --- Extract  NUT from NetCDF dataset --- 
HET = ds["HET"][:,:,:,:]                 
println("HET stats — min: ", minimum(HET),  ", max: ", maximum(HET))
# --- Extract  NUT from NetCDF dataset --- 
NUT = ds["NUT"][:,:,:,:]                 
println("NUT stats — min: ", minimum(NUT),  ", max: ", maximum(NUT))
# --- Extract  DOM from NetCDF dataset --- 
POM = ds["POM"][:,:,:,:]                 
println("POM stats — min: ", minimum(POM),  ", max: ", maximum(POM))
# --- Extract  DOM from NetCDF dataset --- 
DOM = ds["DOM"][:,:,:,:]                 
println("DOM stats — min: ", minimum(DOM),  ", max: ", maximum(DOM))
# --- Extract  O₂ from NetCDF dataset --- 
O₂ = ds["O₂"][:,:,:,:]                 
println("O₂ stats — min: ", minimum(O₂),  ", max: ", maximum(O₂))

@info "BGH arrays loaded"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate additional fields
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Pressure = similar(T[:,:,:,1])  
# Pressure is needed in atm; we calculate  ro*g*z  in Pa and convert to atm using 1 atm = 101325 Pa
for k in 1:Nz
    Pressure[:,:,k] .= 1. + 0.0992*(-depth[k])
end

# oxygen saturation related
O₂_sat_val = similar(T)  # Oxygen saturation concentrtaion
O₂_sat = similar(T)       # Oxygen %
ϵ = eps(Float32)            # small positive number needed in division to zero

# bottom depths index for the bottom maps
bottom_z = ones(Int, size(O₂, 1), size(O₂, 2))
for i = 1:size(O₂, 1)
    for j = 1:size(O₂, 2)
        for k = 1:size(O₂, 3)
            if O₂[i, j, k, 1] .!= 0
                bottom_z[i, j] = k
                break
                if k == Nz
                    bottom_z[i, j] = Nz
                end
            end
            if k == Nz
                bottom_z[i, j] = Nz
            end
        end
    end
end
# Fill oxygen saturation and percentage arrays
for i = 1:size(O₂, 1)
    for j = 1:size(O₂, 2)
        for k = 1:size(O₂, 3)
            for it = 1:size(O₂, 4)
                O₂_sat_val[i, j, k, it] = oxygen_saturation(
                    Float64(T[i, j, k, it]),
                    Float64(S[i, j, k, it]),
                    Float64(Pressure[i, j, k])
                )
                denom = O₂_sat_val[i, j, k, it] == 0f0 ? ϵ : O₂_sat_val[i, j, k, it]
                O₂_sat[i, j, k, it] = 100f0 * O₂[i, j, k, it] / denom
            end
        end
    end
end
println("O₂_sat_val stats — min: ", minimum(O₂_sat_val),  ", max: ", maximum(O₂_sat_val))
println("O₂_sat % stats — min: ", minimum(O₂_sat),  ", max: ", maximum(O₂_sat))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  vertical distributions changes in a point
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 15, 52, times, depth, folder) # Vestfjorden
plot_ztime(NUT, O₂, O₂_sat, P, HET, T, DOM, POM, S, 35, 50, times, depth, folder) # Bunnefjorden

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot parameters MAPs as subplots at given days(plot_dates), depths(depth_indexes) and bottom
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_dates = [36, 72, 108, 144, 180, 226, 262, 298, 334]
bottom_layer = Nz
depth_indexes =  [12]  # surface slice index
fig_width =  1200         # figure width
fig_height = 1200         # figure height

for plot_day in plot_dates
    day_index = plot_day * round(Int, length(times)/365)  
    for depth_index in depth_indexes
        println("Plotting full map figure for day $plot_day ...")

        # Extract horizontal slices
        T_slice = replace_zeros_with_NaN!(T, depth_index, day_index)
        S_slice = replace_zeros_with_NaN!(S, depth_index, day_index)
        O₂_slice = replace_zeros_with_NaN!(O₂, depth_index, day_index)
        NUT_slice = replace_zeros_with_NaN!(NUT, depth_index, day_index)
        P_slice = replace_zeros_with_NaN!(P, depth_index, day_index)
        HET_slice = replace_zeros_with_NaN!(HET, depth_index, day_index)
        DOM_slice = replace_zeros_with_NaN!(DOM, depth_index, day_index)
        POM_slice  = replace_zeros_with_NaN!(POM, depth_index, day_index)
        O₂_sat_slice = replace_zeros_with_NaN!(O₂_sat, depth_index, day_index)
println("Creating figure for day $plot_day at depth index $depth_index ...")

        # Create 3×3 subplot figure
        fig = Figure(size=(fig_width, fig_height))

        plot_tracer_subplot!(fig, (1, 1), T_slice, "T [°C]";  colorrange=(0, 20), colormap=:turbo)
        plot_tracer_subplot!(fig, (1, 3), S_slice, "S [PSU]"; colorrange=(15, 35))
        plot_tracer_subplot!(fig, (1, 5), O₂_slice,"O₂ [μM]"; colorrange=(0, 350), colormap=:turbo)

        plot_tracer_subplot!(fig, (2, 1), P_slice,     "PHY";   colorrange=(0, 5), colormap=Reverse(:cubehelix))
        plot_tracer_subplot!(fig, (2, 3), HET_slice,   "HET";   colorrange=(0, 5), colormap=Reverse(:afmhot))
        plot_tracer_subplot!(fig, (2, 5), O₂_sat_slice,"O₂ [%]";colorrange=(0, 150), colormap=:gist_stern)

        plot_tracer_subplot!(fig, (3, 1), DOM_slice,   "DOM"; colorrange=(0, 15), colormap=Reverse(:CMRmap))
        plot_tracer_subplot!(fig, (3, 3), POM_slice,   "POM"; colorrange=(0, 5), colormap=Reverse(:greenbrownterrain))
        plot_tracer_subplot!(fig, (3, 5), NUT_slice,   "NUT"; colorrange=(0, 40), colormap=Reverse(:cherry))

        save(joinpath(folder, "map_iz_$(depth_index)_day_$(plot_day).png"), fig)
        @info "Saved: map_iz_$(depth_index)_day_$(plot_day).png"
    end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now plot bottom maps
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create 3×3 subplot figure
    fig_b = Figure(size=(fig_width, fig_height))

    # Preallocate result arrays (same horizontal dimensions as O₂)
    T_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    S_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    O₂_slice_bot = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))

    P_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    HET_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    O₂_sat_slice_bot = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))

    DOM_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    POM_slice_bot  = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    NUT_slice_bot = Array{Float64}(undef, size(O₂, 1), size(O₂, 2))
    # Fill them
    for i in 1:size(O₂, 1)
        for j in 1:size(O₂, 2)
            z = bottom_z[i, j]
            T_slice_bot[i, j]  = Float64(T[i, j, z, day_index])
            S_slice_bot[i, j]  = Float64(S[i, j, z, day_index])
            O₂_slice_bot[i, j] = Float64(O₂[i, j, z, day_index])

            P_slice_bot[i, j]  = Float64(P[i, j, z, day_index])
            HET_slice_bot[i, j] = Float64(HET[i, j, z, day_index])
            O₂_sat_slice_bot[i, j] = Float64(O₂_sat[i, j, z, day_index])
        
            DOM_slice_bot[i, j] = Float64(DOM[i, j, z, day_index])
            POM_slice_bot[i, j] = Float64(POM[i, j, z, day_index])
            NUT_slice_bot[i, j] = Float64(NUT[i, j, z, day_index])
                             
        end
    end
    # Replace zeros with NaN
    @. T_slice_bot  = ifelse(T_slice_bot == 0, NaN, T_slice_bot)
    @. S_slice_bot  = ifelse(S_slice_bot == 0, NaN, S_slice_bot)
    @. O₂_slice_bot = ifelse(O₂_slice_bot == 0, NaN, O₂_slice_bot)
    @. P_slice_bot  = ifelse(P_slice_bot == 0, NaN, P_slice_bot)
    @. HET_slice_bot  = ifelse(HET_slice_bot == 0, NaN, HET_slice_bot)
    @. O₂_sat_slice_bot = ifelse(O₂_sat_slice_bot == 0, NaN, O₂_sat_slice_bot)
    @. DOM_slice_bot  = ifelse(DOM_slice_bot == 0, NaN, DOM_slice_bot)
    @. POM_slice_bot  = ifelse(POM_slice_bot == 0, NaN, POM_slice_bot)
    @. NUT_slice_bot = ifelse(NUT_slice_bot == 0, NaN, NUT_slice_bot)

    # Plot bottom maps
    plot_tracer_subplot!(fig_b, (1, 1), T_slice_bot, "T [°C]";  colorrange=(0, 20), colormap=:turbo)
    plot_tracer_subplot!(fig_b, (1, 3), S_slice_bot, "S [PSU]"; colorrange=(15, 35))
    plot_tracer_subplot!(fig_b, (1, 5), O₂_slice_bot, "O₂ [μM]"; colorrange=(0, 350), colormap=:turbo)

    plot_tracer_subplot!(fig_b, (2, 1), P_slice_bot,     "PHY";   colorrange=(0, 5), colormap=Reverse(:cubehelix))
    plot_tracer_subplot!(fig_b, (2, 3), HET_slice_bot,   "HET";   colorrange=(0, 5), colormap=Reverse(:afmhot))
    plot_tracer_subplot!(fig_b, (2, 5), O₂_sat_slice_bot,"O₂ [%]";colorrange=(0, 150), colormap=:gist_stern)

    plot_tracer_subplot!(fig_b, (3, 1), DOM_slice_bot,   "DOM"; colorrange=(0, 15), colormap=Reverse(:CMRmap))
    plot_tracer_subplot!(fig_b, (3, 3), POM_slice_bot,   "POM"; colorrange=(0, 5), colormap=Reverse(:greenbrownterrain))
    plot_tracer_subplot!(fig_b, (3, 5), NUT_slice_bot,   "NUT"; colorrange=(0, 40), colormap=Reverse(:cherry))
    # Save figure with bottom maps
            save(joinpath(folder, "map_bottom_day_$(plot_day).png"), fig_b)
            @info "Saved: map_bottom_day_$(plot_day).png"        
end 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make movies at the bottom        
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
record_bottom_tracer(O₂, "O2_bottom", Nz, times, folder;
    colorrange = (-1, 350), colormap = :turbo, figsize = (400, 550),)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make movies at the surface (iz = Nz) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
 record_horizontal_tracer(NUT, times, folder, "NUTsurf", "Nutrients (μM N)",
                          colorrange = (0, 40),colormap = Reverse(:cherry),iz = Nz, )
 # ~~~~~~~~~~~~~~~~~~~~~~~~~
 record_horizontal_tracer(O₂, times, folder, "O2surf", "Dissolved oxygen (μM)",
     colorrange = (0, 350), colormap = :turbo, iz = Nz, )
 # ~~~~~~~~~~~~~~~~~~~~~~~~~
 record_horizontal_tracer(P, times, folder, "PHYsurf", "PHY (μM)",
     colorrange = (0, 5), colormap = Reverse(:cubehelix), iz = Nz, )     
# ~~~~~~~~~~~~~~~~~~~~~~~~~
 record_horizontal_tracer(HET, times, folder, "HETsurf", "HET (μM)",
     colorrange = (0, 5), colormap = Reverse(:afmhot), iz = Nz, )
# ~~~~~~~~~~~~~~~~~~~~~~~~~
 record_horizontal_tracer(POM, times, folder, "POMsurf", "POM (μM)",
     colorrange = (0, 5), colormap = Reverse(:greenbrownterrain), iz = Nz, )    
# ~~~~~~~~~~~~~~~~~~~~~~~~~
    record_horizontal_tracer(DOM, times, folder, "DOMsurf", "DOM (μM)",
     colorrange = (0, 20), colormap = Reverse(:CMRmap), iz = Nz, )
#println("⏸ Press Enter to continue...")
#readline()

