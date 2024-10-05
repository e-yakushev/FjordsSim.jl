using Oceananigans
using Oceananigans.Architectures
using ClimaOcean

using FileIO
using JLD2

function grid_bathymetry_from_lat_lon(
    arch,
    Nx,
    Ny,
    halo,
    latitude,
    longitude,
    datadir,
    filename,
    depth,
    surface_layer_Δz,
    stretching_factor,
    surface_layer_height,
    height_above_water,
    minimum_depth,
)

    z_faces = stretched_vertical_faces(;
        depth,
        surface_layer_Δz,
        stretching = PowerLawStretching(stretching_factor),
        surface_layer_height,
    )

    Nz = length(z_faces) - 1

    underlying_grid = LatitudeLongitudeGrid(
        arch;
        size = (Nx, Ny, Nz),
        latitude = latitude,
        longitude = longitude,
        z = z_faces,
        halo = halo,
    )

    h = regrid_bathymetry_regional(
        underlying_grid;
        height_above_water = height_above_water,
        minimum_depth = minimum_depth,
        dir = datadir,
        filename = filename,
        interpolation_passes = 1,
        major_basins = 1,
    )
    return underlying_grid, h
end

function grid_from_lat_lon!(sim_setup)
    arch,
    Nx,
    Ny,
    halo,
    latitude,
    longitude,
    datadir,
    filename,
    depth,
    surface_layer_Δz,
    stretching_factor,
    surface_layer_height,
    height_above_water,
    minimum_depth = sim_setup.grid_parameters

    underlying_grid, h = grid_bathymetry_from_lat_lon(
        arch,
        Nx,
        Ny,
        halo,
        latitude,
        longitude,
        datadir,
        filename,
        depth,
        surface_layer_Δz,
        stretching_factor,
        surface_layer_height,
        height_above_water,
        minimum_depth,
    )

    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(h))
    sim_setup.grid[] = grid
    return grid
end

function grid_from_bathymetry_file!(sim_setup)
    arch, Nz, halo, datadir, filename, latitude, longitude = sim_setup.grid_parameters

    filepath_topo = joinpath(datadir, filename)
    @load filepath_topo depth
    Nx, Ny = size(depth)
    depth_max = abs(minimum(depth))
    z_faces = exponential_z_faces(; Nz = Nz, depth = depth_max, h = depth_max)
    underlying_grid = LatitudeLongitudeGrid(
        arch;
        size = (Nx, Ny, Nz),
        halo = halo,
        z = z_faces,
        latitude,
        longitude,
    )
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)
    sim_setup.grid[] = grid
    return grid
end

function grid_latitude_flat!(sim_setup)
    arch, Nx, Ny, Nz, halo, latitude, longitude, depth = sim_setup.grid_parameters

    z_faces = exponential_z_faces(; Nz = Nz, depth = depth, h = depth)
    grid = LatitudeLongitudeGrid(
        arch;
        size = (Nx, Ny, Nz),
        halo = halo,
        z = z_faces,
        latitude,
        longitude,
    )
    sim_setup.grid[] = grid
    return grid
end

function grid_column!(sim_setup)
    arch, Nz, halo, latitude, longitude, depth, h = sim_setup.grid_parameters
    longitude = longitude .+ (-0.03, 0.03)
    latitude = latitude .+ (-0.03, 0.03)

    z_faces = exponential_z_faces(; Nz, depth, h)
    grid = LatitudeLongitudeGrid(
        arch;
        size = (3, 3, Nz),
        halo = halo,
        z = z_faces,
        latitude,
        longitude,
    )
    sim_setup.grid[] = grid
    return grid
end
