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
ds = NCDataset("$filename.nc", "r")

# Get grid dimensions and properties
grid_group = ds.group["grid_reconstruction"]
Nx = Int(grid_group.attrib["Nx"])
Ny = Int(grid_group.attrib["Ny"])
Nz = Int(grid_group.attrib["Nz"])

println("Grid dimensions: Nx=$Nx, Ny=$Ny, Nz=$Nz")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract from NetCDF dataset 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
times = ds["time"][:]         # time in seconds, float
depth = ds["z_aac"][:]        # depth in m, starting from the deepest, negative, float
println("Time stats (in seconds) — min: ", minimum(times), ", max: ", maximum(times))

# --- Extract variables from NetCDF dataset --- 
T = ds["T"][:,:,:,:]               
println("T stats — min: ", minimum(T), ", max: ", maximum(T))

S = ds["S"][:,:,:,:]                
println("S stats — min: ", minimum(S), ", max: ", maximum(S))

P = ds["P"][:,:,:,:]               
println("P stats — min: ", minimum(P),  ", max: ", maximum(P))

HET = ds["HET"][:,:,:,:]                 
println("HET stats — min: ", minimum(HET),  ", max: ", maximum(HET))

NUT = ds["NUT"][:,:,:,:]                 
println("NUT stats — min: ", minimum(NUT),  ", max: ", maximum(NUT))

POM = ds["POM"][:,:,:,:]                 
println("POM stats — min: ", minimum(POM),  ", max: ", maximum(POM))

DOM = ds["DOM"][:,:,:,:]                 
println("DOM stats — min: ", minimum(DOM),  ", max: ", maximum(DOM))

O₂ = ds["O₂"][:,:,:,:]                 
println("O₂ stats — min: ", minimum(O₂),  ", max: ", maximum(O₂))

@info "BGC arrays loaded"

# Example for the Oslofjord (approximate boundaries)
oslo_fjord_lon_min = 10.45   # East longitude
oslo_fjord_lon_max = 10.8
oslo_fjord_lat_min = 59.6  # Northern latitude  
oslo_fjord_lat_max = 59.95

# Create realistic coordinates
real_lon = range(oslo_fjord_lon_min, oslo_fjord_lon_max, length=Nx)
real_lat = range(oslo_fjord_lat_min, oslo_fjord_lat_max, length=Ny)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to convert seconds to days and find indices for 365 days
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function get_365_days_indices(times_seconds, total_days=365)
    # Converting seconds to days
    times_days = times_seconds ./ (24 * 3600)
    
    # Find the maximum time in days
    max_days = maximum(times_days)
    
    if max_days < total_days
        @warn "Available data only for $max_days days, but requested $total_days days"
        total_days = Int(floor(max_days))
    end
    
    # Find indexes for 365 days
    target_times = range(0, total_days, length=total_days)
    indices = Int[]
    
    for target_day in target_times
        # Find the nearest available time step
        idx = argmin(abs.(times_days .- target_day))
        push!(indices, idx)
    end
    
    # Remove duplicates
    unique_indices = unique(indices)
    
    @info "Selected $(length(unique_indices)) frames for $total_days days"
    @info "Time range: day 0 to day $total_days"
    
    return unique_indices
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to plot six animations on one sheet for 365 days
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function plot_six_animations_365_days(tracers, times_seconds, folder, labels, longitudes, latitudes;
                                      colorranges, colormaps, iz=Nz, framerate=0.1, 
                                      figsize=(1800, 1200), fontsize=22)
    
    # Receive indexes for 365 days
    day_indices = get_365_days_indices(times_seconds, 365)
    Nt = length(day_indices)
    
    # Create an Observable for the current day
    day_iter = Observable(1)
    current_frame = Observable(1)
    
    # Convert time to days for display
    times_days = times_seconds ./ (24 * 3600)
    
    # Create the figure
    fig = Figure(size=figsize, fontsize=fontsize)
    
    # Create arrays to store axes and heatmaps
    axes_vec = []
    heatmaps_vec = []
    
    # Create 6 subplots in 2x3 grid
    for i in 1:6
        row = ((i-1) ÷ 3) + 1
        col = ((i-1) % 3) + 1
        
        # Create main plot axis
        ax = Axis(fig[row, col*2-1], 
                  xlabel=(row == 2 ? "Longitude (°E)" : ""),
                  ylabel=(col == 1 ? "Latitude (°N)" : ""),
                  title=labels[i])
        
        # Create observable for current tracer и текущего дня
        Ti = @lift begin
            tracer_data = tracers[i]
            frame_idx = day_indices[$day_iter]
            # Check the dimensions of the data
            if ndims(tracer_data) == 4
                Ti_val = tracer_data[:, :, iz, frame_idx]
            else
                error("Expected 4D array, got $(ndims(tracer_data))D")
            end
            # Convert to Float64 and replace zeros with NaN
            Ti_val = Float64.(Ti_val)
            Ti_val[Ti_val .== 0] .= NaN
            Ti_val
        end
        
        # Create heatmap
        hm = heatmap!(ax, longitudes, latitudes, Ti, 
                     colorrange=colorranges[i], 
                     colormap=colormaps[i], 
                     nan_color=:lightgray)
        
        # Create colorbar next to the plot
        Colorbar(fig[row, col*2], hm, width=20, label="")
        
        push!(axes_vec, ax)
        push!(heatmaps_vec, hm)
    end
    
    # Add super title with time in days
    super_title = @lift begin
        day_idx = day_indices[$day_iter]
        current_day = round(Int, times_days[day_idx])
        "Day $current_day/365 - Surface"
    end
    Label(fig[0, :], super_title, fontsize=22, font=:bold)
    
    # Add progress information
    #progress_info = @lift "Frame $($day_iter)/$Nt"
    #Label(fig[3, :], progress_info, fontsize=22, color=:gray)
    
    # Adjust layout
    rowgap!(fig.layout, 15)
    colgap!(fig.layout, 10)
    
    output_file = joinpath(folder, "six_animations_365days_iz_$iz.gif")
    
    @info "Starting 365-day animation recording with $Nt frames..."
    @info "Output file: $output_file"
    @info "Frame rate: $framerate fps"
    
    # Record animation for selected days
    record(fig, output_file, 1:Nt, framerate=framerate) do i
        day_iter[] = i
        current_frame[] = day_indices[i]
        
        if i % 30 == 0 || i <= 5 || i >= Nt-5
            day_idx = day_indices[i]
            current_day = round(Int, times_days[day_idx])
            @info "Processing: frame $i/$Nt (day $current_day)"
        end
    end
    
    @info "365-day animation completed: $output_file"
    
    return fig
end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare and call the function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare data for the six animations
tracers = [NUT, O₂, P, HET, POM, DOM]

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

# Check the available data
times_days = times ./ (24 * 3600)
max_days = maximum(times_days)
@info "Available data: $(length(times)) time steps, up to day $(round(max_days, digits=1))"

# Select the appropriate function depending on the data
if max_days >= 365
    @info "Enough data for 365-day animation"
    fig = plot_six_animations_365_days(tracers, times, folder, labels, real_lon, real_lat,
                                      colorranges=colorranges, colormaps=colormaps, 
                                      iz=Nz, framerate=24)
else
    @info "Using all available data ($(round(max_days, digits=1)) days)"
    fig = plot_six_animations_daily(tracers, times, folder, labels, real_lon, real_lat,
                                   colorranges=colorranges, colormaps=colormaps, 
                                   iz=Nz, framerate=0.1)
end

# Close the dataset
close(ds)

@info "Script completed!"