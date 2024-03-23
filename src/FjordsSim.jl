module FjordsSim

export FjordsSetup, ImmersedBoundaryGrid

import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

using Oceananigans
using Oceananigans.Architectures

include("DataWrangling.jl")
include("Bathymetry.jl")
include("VerticalGrids.jl")

using .Bathymetry: regrid_bathymetry
using .VerticalGrids

struct FjordsSetup
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

function FjordsSetup(;
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

    return FjordsSetup(
        arch,
        (Nx, Ny, Nz),
        latitude,
        longitude,
        z_faces,
        halo,
        joinpath(homedir(), "fjords_data"),
        filename,
        height_above_water,
        minimum_depth,
    )
end

function ImmersedBoundaryGrid(setup::FjordsSetup)
    grid = LatitudeLongitudeGrid(setup.arch;
                                 size = setup.size,
                                 latitude = setup.latitude,
                                 longitude = setup.longitude,
                                 z = setup.z,
                                 halo = setup.halo)
    h = regrid_bathymetry(
        grid;
        height_above_water=setup.height_above_water,
        minimum_depth=setup.minimum_depth,
        dir=setup.datadir,
        filename = setup.filename,
    )

    return ImmersedBoundaryGrid(grid, GridFittedBottom(h))
end

end # module FjordsSim
