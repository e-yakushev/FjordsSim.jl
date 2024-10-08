module FjordsSim

using Oceananigans
using Oceananigans.Architectures
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity
using ClimaOcean.OceanSimulations:
    default_free_surface,
    default_ocean_closure,
    default_momentum_advection,
    default_tracer_advection,
    u_quadratic_bottom_drag,
    v_quadratic_bottom_drag,
    u_immersed_bottom_drag,
    v_immersed_bottom_drag
using OceanBioME
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

import Oceananigans.Architectures: on_architecture

include("utils.jl")
include("bathymetry.jl")
include("grids.jl")
include("initial_conditions.jl")
include("turbulence_closures.jl")
include("boundary_conditions.jl")
include("forcing.jl")
include("radiation.jl")
include("output.jl")

include("BGCModels/BGCModels.jl")

using .BGCModels: OXYDEP

# there is no a steprangelen method in oceananigans
# but adding it here is type piracy
# we need this when loading atmospheric forcing to the video memory
function on_architecture(::GPU, a::StepRangeLen)
    on_architecture(GPU(), collect(a))
end

free_surface_default(grid) = SplitExplicitFreeSurface(grid[]; cfl = 0.7)
atmosphere_JRA55(arch, backend, grid) = JRA55_prescribed_atmosphere(arch; backend, grid = grid[])
biogeochemistry_LOBSTER(grid) = LOBSTER(; grid = grid[], carbonates = false, open_bottom = false)
biogeochemistry_OXYDEP(grid, args_oxydep) = OXYDEP(;
    grid = grid[],
    args_oxydep...,
    surface_photosynthetically_active_radiation = PAR⁰,
    TS_forced = false,
    Chemicals = false,
    scale_negatives = false,
)

# Grid
grid = Ref{Any}(nothing)

mutable struct SetupModel
    grid_callable!::Function
    grid_parameters::NamedTuple
    grid::Ref
    buoyancy::Any
    closure::Any
    tracer_advection::Any
    momentum_advection::Any
    tracers::Tuple
    initial_conditions::NamedTuple
    free_surface_callable::Function
    free_surface_args::Tuple
    coriolis::Any
    forcing_callable::Any
    forcing_args::Any
    bc_callable::Any
    bc_args::Any
    atmosphere_callable::Any
    atmosphere_args::Any
    radiation::Any
    biogeochemistry_callable::Any
    biogeochemistry_args::Any
end

function coupled_hydrostatic_simulation(sim_setup::SetupModel)
    grid = sim_setup.grid_callable!(sim_setup)
    buoyancy = sim_setup.buoyancy
    closure = sim_setup.closure
    tracer_advection = sim_setup.tracer_advection
    momentum_advection = sim_setup.momentum_advection
    tracers = sim_setup.tracers
    free_surface = sim_setup.free_surface_callable(sim_setup.free_surface_args...)
    coriolis = sim_setup.coriolis
    forcing = sim_setup.forcing_callable(sim_setup.forcing_args...)
    boundary_conditions = sim_setup.bc_callable(sim_setup.bc_args...)
    biogeochemistry = safe_execute(sim_setup.biogeochemistry_callable)(sim_setup.biogeochemistry_args...)

    ## Model
    ocean_model = HydrostaticFreeSurfaceModel(;
        grid,
        buoyancy,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        free_surface,
        coriolis,
        forcing,
        boundary_conditions,
        biogeochemistry
    )
    println("Done compiling HydrostaticFreeSurfaceModel")

    ## Simulation
    Δt = 1seconds
    ocean_sim = Simulation(ocean_model; Δt)
    println("Initialized simulation")

    ## Set initial conditions
    set!(ocean_model; sim_setup.initial_conditions...)
    println("Set initial conditions")

    ## Coupled model / simulation
    sea_ice = nothing
    atmosphere = safe_execute(sim_setup.atmosphere_callable)(sim_setup.atmosphere_args...)
    println("Initialized atmosphere")
    radiation = sim_setup.radiation
    coupled_model = OceanSeaIceModel(ocean_sim, sea_ice; atmosphere, radiation)
    println("Initialized coupled model")
    coupled_simulation = Simulation(coupled_model; Δt)
    println("Initialized coupled simulation")
    return coupled_simulation
end

end # module FjordsSim
