module FjordsSim

export SetupGridRegrid, SetupGridPredefinedFromFile, ImmersedBoundaryGrid, OXYDEP

import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

using FileIO
using Oceananigans
using Oceananigans.Architectures
using JLD2

include("DataWrangling.jl")
include("Bathymetry.jl")
include("VerticalGrids.jl")
include("OxydepModel.jl")

using .Bathymetry: regrid_bathymetry
using .VerticalGrids
using .OxydepModel

struct SetupGridRegrid
    arch::AbstractSerialArchitecture
    size::Tuple
    latitude::Tuple
    longitude::Tuple
    z
    halo::Tuple
    datadir::String
    filename::String
    height_above_water
    minimum_depth::Number
end

function SetupGridRegrid(;
    Nx::Integer,
    Ny::Integer,
    latitude::Tuple,
    longitude::Tuple,
    datadir:: String = joinpath(homedir(), "data_fjords"),
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
        surface_layer_height
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

struct SetupGridPredefinedFromFile
    arch::AbstractSerialArchitecture
    Nx::Integer
    Ny::Integer
    dx::Integer
    dy::Integer
    z_levels::Vector
    z_middle::Vector
    halo::Tuple
    datadir::String
    filename::String
end

function ImmersedBoundaryGrid(setup_grid::SetupGridRegrid)
    grid = LatitudeLongitudeGrid(setup_grid.arch;
                                 size = setup_grid.size,
                                 latitude = setup_grid.latitude,
                                 longitude = setup_grid.longitude,
                                 z = setup_grid.z,
                                 halo = setup_grid.halo)
    h = regrid_bathymetry(
        grid;
        height_above_water=setup_grid.height_above_water,
        minimum_depth=setup_grid.minimum_depth,
        dir=setup_grid.datadir,
        filename = setup_grid.filename,
    )

    return ImmersedBoundaryGrid(grid, GridFittedBottom(h))
end

function ImmersedBoundaryGrid(setup_grid::SetupGridPredefinedFromFile)
    Nz = length(setup_grid.z_levels) - 1
    underlying_grid = RectilinearGrid(
        setup_grid.arch,
        topology = (Bounded, Bounded, Bounded),
        size = (setup_grid.Nx, setup_grid.Ny, Nz),
        x = (0, setup_grid.dx * setup_grid.Nx),
        y = (0, setup_grid.dy * setup_grid.Ny),
        z = setup_grid.z_levels,
        halo = setup_grid.halo,
    )
    filepath_topo = joinpath(setup_grid.datadir, setup_grid.filename)
    @load filepath_topo depth
    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)
end

end # module FjordsSim
