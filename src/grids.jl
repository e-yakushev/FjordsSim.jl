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

# function ImmersedBoundaryGrid(setup_grid::SetupGridRegrid)
#     grid = LatitudeLongitudeGrid(
#         setup_grid.arch;
#         size = setup_grid.size,
#         latitude = setup_grid.latitude,
#         longitude = setup_grid.longitude,
#         z = setup_grid.z,
#         halo = setup_grid.halo,
#     )
#     h = regrid_bathymetry(
#         grid;
#         height_above_water = setup_grid.height_above_water,
#         minimum_depth = setup_grid.minimum_depth,
#         dir = setup_grid.datadir,
#         filename = setup_grid.filename,
#     )
#     return ImmersedBoundaryGrid(grid, GridFittedBottom(h))
# end

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
        halo = (7, 7, 7),
        z = z_faces,
        latitude,
        longitude,
    )
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)
    sim_setup.grid[] = grid
    return grid
end
