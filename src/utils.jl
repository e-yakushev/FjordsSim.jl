using Oceananigans.Fields: interior
using Oceananigans.OutputReaders: FieldTimeSeries, OnDisk
using Oceananigans.Utils: prettytime, pretty_filesize
using NCDatasets: Dataset
using JLD2: @save
using Printf: @sprintf

"""
# Example usage
filename = "./app/varna/Varna_brom.nc"
T, S, U, depth, time = read_TSU_forcing(filename)
"""
function read_TSU_forcing(filename::String)
    # reads temperature, salinity and velocities from netcdf.
    # currently for _brom.nc files prepared for 2DBP

    # Open the NetCDF file
    ds = Dataset(filename, "r")

    # Read the variables
    temperature = reverse(ds["temp"][:, :], dims = 2)
    salinity = reverse(ds["salt"][:, :], dims = 2)
    velocity = reverse(ds["u"][:, :], dims = 2)
    Kz = reverse(ds["Kz"][:, :], dims = 2)

    depth = reverse(ds["depth"][:])
    time = ds["oc_time"][:]

    # Close the dataset
    close(ds)

    return temperature, salinity, velocity, Kz, depth, time
end

"""
download_progress
"""
function download_progress(total, now; filename = "")
    messages = 10

    if total > 0
        fraction = now / total

        if fraction < 1 / messages && next_fraction[] == 0
            @info @sprintf("Downloading %s (size: %s)...", filename, pretty_filesize(total))
            next_fraction[] = 1 / messages
            download_start_time[] = time_ns()
        end

        if fraction > next_fraction[]
            elapsed = time_ns() - download_start_time[]
            msg = @sprintf(
                " ... downloaded %s (%d%% complete, %s)",
                pretty_filesize(now),
                100fraction,
                prettytime(elapsed)
            )
            @info msg
            next_fraction[] = next_fraction[] + 1 / messages
        end
    else
        if now > 0 && next_fraction[] == 0
            @info "Downloading $filename..."
            next_fraction[] = 1 / messages
            download_start_time[] = time_ns()
        end
    end

    return nothing
end

wall_time = [time_ns()]
function progress(sim)
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf(
        "Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), max(T): %.2f, min(T): %.2f, wtime: %s \n",
        prettytime(ocean.model.clock.time),
        ocean.model.clock.iteration,
        prettytime(ocean.Δt),
        umax...,
        Tmax,
        Tmin,
        prettytime(step_time)
    )

    wall_time[1] = time_ns()
end

function safe_execute(callable)
    return function (args...)
        if callable === nothing || args === nothing
            return nothing
        elseif isa(callable, Function)
            return callable(args...)
        else
            return nothing
        end
    end
end

function extract_z_faces(grid)
    bar = grid["zᵃᵃᶜ"]
    zero_index = findfirst(x -> x > 0.0, bar)
    n = grid["Nz"] + 1
    if zero_index > 1
        start_index = max(1, zero_index - n)
        z = bar[start_index:zero_index-1]
    else
        z = []
    end
    return z
end


function netcdf_to_jld2(netcdf_file::String, jld2_file::String)
    # Open the NetCDF file in read-only mode
    ds = Dataset(netcdf_file, "r")

    # Create a dictionary to store the variables from the NetCDF file
    data_dict = Dict()

    # Loop through all variables in the NetCDF file and store them in the dictionary
    for varname in keys(ds)
        data_dict[varname] = convert(Array, ds[varname])  #syntax to keep the dimensions
        print(size(convert(Array, ds[varname])))
    end

    # Save the data to a JLD2 file
    @save jld2_file data_dict

    # Close the NetCDF dataset
    close(ds)

    println("Conversion completed: NetCDF to JLD2")
end

function save_fts(; jld2_filepath, fts_name, fts, grid, times, boundary_conditions)
    isfile(jld2_filepath) && rm(jld2_filepath)
    on_disk_fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        boundary_conditions,
        backend = OnDisk(),
        path = jld2_filepath,
        name = fts_name,
    )
    for i = 1:size(fts)[end]
        set!(on_disk_fts, fts[i], i, times[i])
    end
end

