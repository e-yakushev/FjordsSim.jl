using FileIO
using JLD2
using ClimaOcean
using Oceananigans
using Oceananigans.Architectures

# import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

struct SetupGridRegrid
    arch::AbstractSerialArchitecture
    size::Tuple
    latitude::Tuple
    longitude::Tuple
    z::Any
    halo::Tuple
    datadir::String
    filename::String
    height_above_water::Any
    minimum_depth::Number
end

function SetupGridRegrid(;
    Nx::Integer,
    Ny::Integer,
    latitude::Tuple,
    longitude::Tuple,
    datadir::String = joinpath(homedir(), "data_fjords"),
    filename::String,
    depth,
    surface_layer_Δz,
    stretching_factor,
    surface_layer_height,
    arch::AbstractSerialArchitecture = CPU(),
    halo::Tuple = (4, 4, 4),
    height_above_water = nothing,
    minimum_depth::Number = 0,
)

    z_faces = stretched_vertical_faces(;
        depth,
        surface_layer_Δz,
        stretching = PowerLawStretching(stretching_factor),
        surface_layer_height,
    )

    Nz = length(z_faces) - 1

    return SetupGridRegrid(
        arch,
        (Nx, Ny, Nz),
        latitude,
        longitude,
        z_faces,
        halo,
        datadir,
        filename,
        height_above_water,
        minimum_depth,
    )
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
