module FjordsSim

using Oceananigans.Architectures
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using OceanBioME
using OceanBioME.Light:
    update_TwoBandPhotosyntheticallyActiveRadiation!,
    TwoBandPhotosyntheticallyActiveRadiation
using ClimaOcean.OceanSimulations:
    default_free_surface,
    default_ocean_closure,
    default_momentum_advection,
    default_tracer_advection,
    u_quadratic_bottom_drag,
    v_quadratic_bottom_drag,
    u_immersed_bottom_drag,
    v_immersed_bottom_drag
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

import Oceananigans.Architectures: on_architecture
import Oceananigans.Biogeochemistry: update_biogeochemical_state!

include("utils.jl")
include("bathymetry.jl")
include("grids.jl")
include("initial_conditions.jl")
include("turbulence_closures.jl")
include("boundary_conditions.jl")
include("forcing.jl")
include("simulations.jl")
include("radiation.jl")

include("BGCModels/BGCModels.jl")

using .BGCModels: OXYDEP

function update_biogeochemical_state!(model, PAR::TwoBandPhotosyntheticallyActiveRadiation)
    arch = architecture(model.grid)
    launch!(
        arch,
        model.grid,
        :xy,
        update_TwoBandPhotosyntheticallyActiveRadiation!,
        PAR.field,
        model.grid,
        model.tracers.PHY,
        PAR.surface_PAR,
        model.clock.time,
        PAR,
    )
end

# there is no a steprangelen method in oceananigans
# but adding it here is type piracy
# we need this when loading atmospheric forcing to the video memory
function on_architecture(::GPU, a::StepRangeLen)
    on_architecture(GPU(), collect(a))
end

end # module FjordsSim
