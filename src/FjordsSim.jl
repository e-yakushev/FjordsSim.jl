module FjordsSim

include("utils.jl")
include("bathymetry.jl")
include("grids.jl")
include("initial_conditions.jl")
include("turbulence_closures.jl")
include("boundary_conditions.jl")
include("forcing.jl")

# include("BGCModels/BGCModels.jl")

# using .BGCModels: OXYDEP

end # module FjordsSim
