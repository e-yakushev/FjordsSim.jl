using Oceananigans.Fields: interior
using Oceananigans.OutputReaders: FieldTimeSeries, OnDisk
using Oceananigans.Utils: prettytime, pretty_filesize
using NCDatasets: Dataset
using JLD2: @save
using Printf: @sprintf

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
    ds = Dataset(netcdf_file, "r")
    data_dict = Dict()
    for varname in keys(ds)
        data_dict[varname] = convert(Array, ds[varname])
        print(size(convert(Array, ds[varname])))
    end

    @save jld2_file data_dict
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

