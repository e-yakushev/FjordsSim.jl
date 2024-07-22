# src/ForcingInput.jl
module ForcingInput

using NCDatasets
using Oceananigans.Fields: interpolate

function read_TSU_forcing(filename::String)
    """
    # Example usage
    filename = "./app/varna/Varna_brom.nc"
    T, S, U = read_TSU_forcing(filename)
    """
    # reads temperature, salinity and velocities from netcdf.
    # currently for _brom.nc files prepared for 2DBP

    # Open the NetCDF file
    ds = Dataset(filename, "r")
    
    # Ensure that the variables "temp", "salt" and "u" exist
    if !haskey(ds, "temp") || !haskey(ds, "salt") || !haskey(ds, "u")
        error("The NetCDF file does not contain the required variables 'temp', 'salt' and 'u'.")
    end
    
    # Read the variables
    temperature = ds["temp"][:, :]
    salinity = ds["salt"][:, :]
    velocity = ds["u"][:, :]
    
    # Close the dataset
    close(ds)
    
    return temperature, salinity, velocity
end

end