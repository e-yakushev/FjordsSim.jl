module FjordsSim

export FjordsSetup, ImmersedBoundaryGrid

import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

using Oceananigans
using Oceananigans.Architectures

include("DataWrangling.jl")
include("Bathymetry.jl")

using .Bathymetry: regrid_bathymetry

struct FjordsSetup
    arch::AbstractSerialArchitecture
    size::Tuple
    latitude::Tuple
    longitude::Tuple
    z::Tuple
    halo::Tuple
    datadir::String
    filename::String
    height_above_water
    minimum_depth::Number
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
