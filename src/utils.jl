using NCDatasets
using Printf

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

    depth = reverse(ds["depth"][:])
    time = ds["oc_time"][:]

    # Close the dataset
    close(ds)

    return temperature, salinity, velocity, depth, time
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
