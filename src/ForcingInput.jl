# src/ForcingInput.jl
module ForcingInput

using NCDatasets
using Oceananigans.Fields
using Oceananigans

export read_TSU_forcing, make_input_grid, regrid_tracer_depth

function read_TSU_forcing(filename::String)
    """
    # Example usage
    filename = "./app/varna/Varna_brom.nc"
    T, S, U, depth, time = read_TSU_forcing(filename)
    """
    # reads temperature, salinity and velocities from netcdf.
    # currently for _brom.nc files prepared for 2DBP

    # Open the NetCDF file
    ds = Dataset(filename, "r")
    
    # Read the variables
    temperature = reverse(ds["temp"][:,:], dims=2)
    salinity = reverse(ds["salt"][:, :], dims=2)
    velocity = reverse(ds["u"][:, :], dims=2)

    depth = reverse(ds["depth"][:])
    time = ds["oc_time"][:]
    
    # Close the dataset
    close(ds)
    
    return temperature, salinity, velocity, depth, time
end


function make_input_grid(z_faces)
    # Determine the number of elements in z_faces
    Nz = length(z_faces) - 1
    
    # Define the input grid
    input_grid = RectilinearGrid(size = (1, 1, Nz), 
                                 x = (0, 500), 
                                 y = (0, 500), 
                                 z = z_faces, 
                                 topology = (Bounded, Bounded, Bounded))
    return input_grid
end

function regrid_tracer(tracer, input_grid, output_grid)

    Nz = input_grid.Nz
    # Create the input field on the input grid
    input_field = CenterField(input_grid)
    
    # Assign the tracer values on the current timestep to the input field
    input_field[1, 1, 1:Nz] = tracer
    
    output_field = CenterField(output_grid)
    regrid!(output_field, input_field)
    
    return output_field
end


function regrid_tracer_over_timesteps(tracer, input_grid, output_grid)
    # Determine the number of timesteps
    Nt = size(tracer, 1)  # assumes time is the first dimension!
    Nz_out = output_grid.Nz
    # Preallocate an array to hold the reshaped data
    result_array = Array{Float64}(undef, Nt, Nz_out)
    
    # Loop over each timestep and apply the regridding function
    for t in 1:Nt
        output_field = regrid_tracer(tracer[t,:], input_grid, output_grid)
        
        # Extract and reshape data from the output field
        result_array[t, :] = reshape(output_field.data[1, 1, 1:Nz_out], Nz_out)
    end
    
    return result_array

    # # Preallocate an array to hold the output fields for each timestep
    # output_fields = Vector{Any}()
    
    # # Loop over each timestep and apply the regridding function
    # for t in 1:Nt
    #     push!(output_fields, regrid_tracer(tracer[t,:], input_grid, output_grid))
    # end
    
    # return output_fields
end


end
