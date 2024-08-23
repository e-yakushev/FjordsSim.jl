module FjordsSim

using Oceananigans.Architectures

import Oceananigans.Architectures: on_architecture

function on_architecture(::GPU, a::StepRangeLen)
    on_architecture(GPU(), collect(a))
end

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
