using Oceananigans.OutputReaders: FieldTimeSeries, Cyclical, AbstractInMemoryBackend, FlavorOfFTS, time_indices
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: interior
using Oceananigans.Forcings: Forcing
using Oceananigans.Grids: Center
using Oceananigans.Units: hours
using ClimaOcean.OceanSimulations: u_immersed_bottom_drag, v_immersed_bottom_drag
using Dates: DateTime, Year, Second
using Statistics: mean
using NCDatasets: Dataset

import Oceananigans.Fields: set!
import Oceananigans.OutputReaders: new_backend

struct NetCDFBackend <: AbstractInMemoryBackend{Int}
    start::Int
    length::Int
end

NetCDFBackend(length) = NetCDFBackend(1, length)
new_backend(::NetCDFBackend, start, length) = NetCDFBackend(start, length)

Base.length(backend::NetCDFBackend) = backend.length
Base.summary(backend::NetCDFBackend) = string("JLD2Backend(", backend.start, ", ", backend.length, ")")

const NetCDFFTS = FlavorOfFTS{<:Any,<:Any,<:Any,<:Any,<:NetCDFBackend}
LX, LY, LZ = Center, Center, Center

function native_times_to_seconds(native_times, start_time = native_times[1])
    times = []
    for native_time in native_times
        time = native_time - start_time
        time = Second(time).value
        push!(times, time)
    end
    return times
end

function load_from_netcdf(; path::String, var_name::String, grid_size::Tuple, time_indices_in_memory::Tuple)
    ds = Dataset(path)
    var = ds[var_name]
    native_times = ds["time"]

    data = fill(-999.0, (grid_size[1:end]..., length(time_indices_in_memory)))
    j = 1
    for i in time_indices_in_memory
        data[:, :, :, j] .= var[:, :, :, i]
        j += 1
    end
    times = convert.(eltype(data), native_times_to_seconds(native_times))

    close(ds)
    return data, times
end

function set!(fts::NetCDFFTS, path::String = fts.path, name::String = fts.name)
    ti = time_indices(fts)
    data, _ = load_from_netcdf(; path, var_name = name, grid_size = size(fts)[1:end-1], time_indices_in_memory = ti)

    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    return nothing
end

function fts_tracer_forcing_func(i, j, k, grid, clock, fields, parameters)
    tr = @inbounds parameters.fts[i, j, k, Time(clock.time)]
    λopen = @inbounds parameters.ftsλ[i, j, k, Time(clock.time)]
    condition = !(tr < -990.0)
    radiation_term = -λopen * (fields.T[i, j, k] - tr)
    return @inbounds ifelse(condition, radiation_term, 0)
end

function forcing_get_tuple(filepath, var_name, grid, time_indices_in_memory, backend)
    grid_size = size(grid)
    data, times = load_from_netcdf(; path = filepath, var_name, grid_size, time_indices_in_memory)
    dataλ, times = load_from_netcdf(; path = filepath, var_name = var_name * "_lambda", grid_size, time_indices_in_memory)
    boundary_conditions = FieldBoundaryConditions(grid, (LX, LY, LZ))

    fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
        path = filepath,
        name = var_name,
    )
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    ftsλ = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        backend,
        time_indexing = Cyclical(),
        boundary_conditions,
        path = filepath,
        name = var_name * "_lambda",
    )
    copyto!(interior(ftsλ, :, :, :, :), dataλ)
    fill_halo_regions!(ftsλ)

    _forcing = Forcing(fts_tracer_forcing_func; discrete_form = true, parameters = (fts = fts, ftsλ = ftsλ))

    result = NamedTuple{(Symbol(var_name),)}((_forcing,))
    return result
end

function forcing_from_file(grid_ref, filepath, tracers)
    grid = grid_ref[]
    ds = Dataset(filepath)
    grid.underlying_grid.Nx == ds.dim["Nx"] &&
        grid.underlying_grid.Ny == ds.dim["Ny"] &&
        grid.underlying_grid.Nz == ds.dim["Nz"] ||
        throw(DimensionMismatch("forcing file dimensions not equal to grid dimensions"))
    tracer_names = map(String, tracers) ∩ keys(ds)
    close(ds)

    backend = NetCDFBackend(2)
    time_indices_in_memory = (1, length(backend))
    result = mapreduce(
        var_name -> forcing_get_tuple(filepath, var_name, grid, time_indices_in_memory, backend),
        merge,
        tracer_names,
    )
    return result
end

function forcing_bottom_drag(bottom_drag_coefficient)
    Fu = Forcing(u_immersed_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    Fv = Forcing(v_immersed_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    return (u = Fu, v = Fv)
end
