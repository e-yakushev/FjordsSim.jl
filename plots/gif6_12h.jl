# Increase the memory limit
ENV["JULIA_GC_SOFT_MEMORY_LIMIT"] = "80%"  # Use up to 80% of memory
using Oceananigans
using JLD2
using NCDatasets
using Printf
using Oceananigans.Units
using Oceananigans.Utils: prettytime
using CairoMakie

using FjordsSim:
    plot_1d_phys,
    extract_z_faces,
    record_vertical_tracer,
    record_surface_speed,
    BGCModel,
    oxygen_saturation  

include("/home/ocean12400/src/FjordsSim.jl/src/BGCModels/BGCModels.jl")
using .BGCModels
include("/home/ocean12400/src/FjordsSim.jl/src/BGCModels/boundary_conditions.jl")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAIN CODE STARTS HERE: open file and extract data to "ds"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
folder = joinpath(homedir(), "FjordsSim_results", "oslofjord")
filename = joinpath(folder, "snapshots_ocean")

# Example for the Oslofjord (approximate boundaries)
oslo_fjord_lon_min = 10.45   # East longitude
oslo_fjord_lon_max = 10.8
oslo_fjord_lat_min = 59.6  # Northern latitude  
oslo_fjord_lat_max = 59.95

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to get time indices every 12 hours
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function get_12hour_indices(times_seconds)
    times_hours = times_seconds ./ 3600
    max_hours = maximum(times_hours)
    target_hours = 0:12:max_hours
    
    indices = Int[]
    sizehint!(indices, length(target_hours))
    
    for target_hour in target_hours
        idx = argmin(abs.(times_hours .- target_hour))
        push!(indices, idx)
    end
    
    unique_indices = unique(indices)
    
    total_days = round(max_hours / 24, digits=1)
    @info "Selected $(length(unique_indices)) frames every 12 hours"
    @info "Time range: 0 to $(round(max_hours, digits=1)) hours ($total_days days)"
    
    return unique_indices
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to format time for display
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function format_time_decimal_days(total_seconds)
    total_days = total_seconds / (24 * 3600)
    return @sprintf("Day %.1f", total_days)
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Optimized function to plot six animations every 12 hours
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_six_animations_12hour_optimized(filename, folder, labels, lon_range, lat_range;
                                   colorranges, colormaps, framerate=1, 
                                   figsize=(1500, 1000), fontsize=22)
    
    # Open the file only to read the required data
    ds = NCDataset(filename * ".nc", "r")
    
    try
        # Receive only the necessary metadata
        grid_group = ds.group["grid_reconstruction"]
        Nx = Int(grid_group.attrib["Nx"])
        Ny = Int(grid_group.attrib["Ny"])
        Nz = Int(grid_group.attrib["Nz"])
        
        println("Grid dimensions: Nx=$Nx, Ny=$Ny, Nz=$Nz")
        
        # Load only the necessary data
        times = ds["time"][:]
        depth = ds["z_aac"][:]
        
        # Get indexes for 12-hour intervals BEFORE loading large arrays
        hour_indices = get_12hour_indices(times)
        Nt = length(hour_indices)
        
        # Create realistic coordinates
        real_lon = range(lon_range[1], lon_range[2], length=Nx)
        real_lat = range(lat_range[1], lat_range[2], length=Ny)

        # Load only selected time steps for each track
        tracer_names = ["NUT", "O₂", "P", "HET", "POM", "DOM"]
        tracers = Vector{Array{Float32, 4}}(undef, 6)
        
        for (idx, name) in enumerate(tracer_names)
            # Load all spatial data, but only selected time steps
            full_data = ds[name][:, :, :, hour_indices]
            tracers[idx] = full_data
            println("$name stats — min: ", minimum(full_data), ", max: ", maximum(full_data))
        end

        times_selected = times[hour_indices]
        times_hours = times_selected ./ 3600
        
        output_file = joinpath(folder, "six_animations_12hour_iz_$Nz.gif")
        
        @info "Starting optimized 12-hour interval animation recording with $Nt frames..."
        @info "Output file: $output_file"
        @info "Frame rate: $framerate fps"
        
        # Optimized animation recording
        record(Figure(size=figsize, fontsize=fontsize), output_file, 1:Nt, framerate=framerate) do i
            fig = current_figure()
            current_time_seconds = times_selected[i]
            current_hours = times_hours[i]
            current_day = floor(Int, current_time_seconds/(24*3600))
    
            # Pre-allocate arrays for data
            surface_data = Matrix{Float64}(undef, Nx, Ny)
            
            # Create 6 subplots in 2x3 grid
            for j in 1:6
                row = ((j-1) ÷ 3) + 1
                col = ((j-1) % 3) + 1
                
                # Create main plot axis
                ax = Axis(fig[row, col*2-1], 
                          xlabel=(row == 2 ? "Longitude (°E)" : ""),
                          ylabel=(col == 1 ? "Latitude (°N)" : ""),
                          title=labels[j])
                
                # Copy surface data directly, avoiding intermediate arrays
                surface_view = @view tracers[j][:, :, Nz, i]
                copyto!(surface_data, surface_view)
                
                # Replacing zeros with NaNs more efficiently
                for idx in eachindex(surface_data)
                    if surface_data[idx] == 0
                        surface_data[idx] = NaN
                    end
                end
                
                # Create heatmap
                hm = heatmap!(ax, real_lon, real_lat, surface_data, 
                             colorrange=colorranges[j], 
                             colormap=colormaps[j], 
                             nan_color=:lightgray)
                
                # Create colorbar next to the plot
                Colorbar(fig[row, col*2], hm, width=20, label="")
            end
            
            # Add super title with time information
            if !isdefined(Main, :title_label) || i == 1
            global title_label = Label(fig[0, :], "Day $current_day - Surface", 
                              fontsize=22, font=:bold)
    
            else
                 title_label.text = "Day $current_day - Surface"
            end

            # Adjust layout
            rowgap!(fig.layout, 15)
            colgap!(fig.layout, 10)
            
            if i % 8 == 0 || i <= 12 || i >= Nt-12
                current_days = current_hours / 24
                @info "Processing: frame $i/$Nt (day $(round(current_days, digits=1)))"
            end
        end
        
        @info "Optimized 12-hour interval animation completed: $output_file"
        
    finally
        close(ds)
    end
    
    return nothing
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare and call the optimized function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

labels = ["NUT [μM N]", "O₂ [μM N]", "PHY [μM N]", "HET [μM N]", "POM [μM N]", "DOM [μM N]"]

colorranges = [
    (0.0, 40.0),    # NUT
    (0.0, 350.0),   # O₂
    (0.0, 5.0),     # PHY
    (0.0, 5.0),     # HET  
    (0.0, 5.0),     # POM
    (0.0, 20.0)     # DOM
]

colormaps = [
    :viridis,        # NUT
    :turbo,          # O₂
    :plasma,         # PHY
    :hot,            # HET
    :rainbow,        # POM
    :jet             # DOM
]

lon_range = (oslo_fjord_lon_min, oslo_fjord_lon_max)
lat_range = (oslo_fjord_lat_min, oslo_fjord_lat_max)

# Call the optimized function
@info "Starting optimized animation..."

plot_six_animations_12hour_optimized(filename, folder, labels, lon_range, lat_range,
                                   colorranges=colorranges, colormaps=colormaps, 
                                   framerate=1)

@info "Script completed!"