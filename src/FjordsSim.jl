module FjordsSim

export SetupGridStretchedVericalFaces, ImmersedBoundaryGrid

import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

using Oceananigans
using Oceananigans.Architectures

include("DataWrangling.jl")
include("Bathymetry.jl")
include("VerticalGrids.jl")

using .Bathymetry: regrid_bathymetry
using .VerticalGrids

struct SetupGridStretchedVericalFaces
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

function SetupGridStretchedVericalFaces(;
    Nx::Integer,
    Ny::Integer,
    latitude::Tuple,
    longitude::Tuple,
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

    return SetupGridStretchedVericalFaces(
        arch,
        (Nx, Ny, Nz),
        latitude,
        longitude,
        z_faces,
        halo,
        joinpath(homedir(), "data_fjords"),
        filename,
        height_above_water,
        minimum_depth,
    )
end

function ImmersedBoundaryGrid(setup_grid::SetupGridStretchedVericalFaces)
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

end # module FjordsSim
