module FjordsSim

export SetupGridRegrid, SetupGridPredefinedFromFile, ImmersedBoundaryGrid, OXYDEP

import Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid

include("utils.jl")
include("bathymetry.jl")
include("grids.jl")

include("BGCModels/BGCModels.jl")

using .BGCModels: OXYDEP

end # module FjordsSim
