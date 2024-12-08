module FjordsSim

using Oceananigans
using Oceananigans.Architectures
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities: CATKEVerticalDiffusivity
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
using Oceananigans.Models: update_model_field_time_series!
using Oceananigans.TimeSteppers: tick!
using ClimaOcean
using ClimaOcean.OceanSeaIceModels: MinimumTemperatureSeaIce
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: compute_atmosphere_ocean_fluxes!
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
import Oceananigans.TimeSteppers: time_step!, update_state!
import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: compute_sea_ice_ocean_fluxes!
import ClimaOcean: SimilarityTheoryTurbulentFluxes

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

const OceanOnlyModel = OceanSeaIceModel{Nothing}
const OceanSimplifiedSeaIceModel = OceanSeaIceModel{<:MinimumTemperatureSeaIce}

const NoSeaIceModel = Union{OceanOnlyModel,OceanSimplifiedSeaIceModel}

# there is no a steprangelen method in oceananigans
# but adding it here is type piracy
# we need this when loading atmospheric forcing to the video memory
function on_architecture(::GPU, a::StepRangeLen)
    on_architecture(GPU(), collect(a))
end

free_surface_default(grid_ref) = SplitExplicitFreeSurface(grid_ref[]; cfl = 0.7)

atmosphere_JRA55(arch, backend, grid_ref, start, stop) =
    JRA55_prescribed_atmosphere(arch, start:stop; backend, grid = grid_ref[])
biogeochemistry_LOBSTER(grid_ref) = LOBSTER(; grid = grid_ref[], carbonates = false, open_bottom = false)
biogeochemistry_OXYDEP(grid_ref, args_oxydep) = OXYDEP(;
    grid = grid_ref[],
    args_oxydep...,
    surface_photosynthetically_active_radiation = PAR⁰,
    TS_forced = false,
    Chemicals = false,
    scale_negatives = false,
)
SimilarityTheoryTurbulentFluxes(; grid_ref::Ref, kw...) = SimilarityTheoryTurbulentFluxes(grid_ref[]; kw...)

# Grid
grid_ref = Ref{Any}(nothing)

mutable struct SetupModel
    grid_callable::Function
    grid_args::NamedTuple
    grid_ref::Ref
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
    similarity_theory_callable::Any
    similarity_theory_args::Any
    biogeochemistry_callable::Any
    biogeochemistry_args::Any
    results_dir::String
end

function SetupModel(
    grid_callable,
    grid_args,
    grid_ref,
    buoyance,
    closure,
    tracer_advection,
    momentum_advection,
    tracers,
    initial_conditions,
    free_surface_callable,
    free_surface_args,
    coriolis,
    forcing_callable,
    forcing_args,
    bc_callable,
    bc_args,
    atmosphere_callable,
    atmosphere_args,
    radiation,
    similarity_theory_callable,
    similarity_theory_args,
    biogeochemistry_callable,
    biogeochemistry_args;
    results_dir = joinpath(homedir(), "FjordsSim_results"),
)

    !isdir(results_dir) && mkdir(results_dir)

    SetupModel(
        grid_callable,
        grid_args,
        grid_ref,
        buoyance,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        initial_conditions,
        free_surface_callable,
        free_surface_args,
        coriolis,
        forcing_callable,
        forcing_args,
        bc_callable,
        bc_args,
        atmosphere_callable,
        atmosphere_args,
        radiation,
        similarity_theory_callable,
        similarity_theory_args,
        biogeochemistry_callable,
        biogeochemistry_args,
        results_dir,
    )
end


function coupled_hydrostatic_simulation(sim_setup::SetupModel)
    grid = sim_setup.grid_callable(sim_setup.grid_args...)
    sim_setup.grid_ref[] = grid
    buoyancy = sim_setup.buoyancy
    closure = sim_setup.closure
    tracer_advection = sim_setup.tracer_advection
    momentum_advection = sim_setup.momentum_advection
    tracers = sim_setup.tracers
    free_surface = sim_setup.free_surface_callable(sim_setup.free_surface_args...)
    coriolis = sim_setup.coriolis
    forcing = sim_setup.forcing_callable(sim_setup.forcing_args...)
    boundary_conditions = sim_setup.bc_callable(sim_setup.bc_args...)
    biogeochemistry =
        safe_execute(sim_setup.biogeochemistry_callable)(sim_setup.biogeochemistry_args...)

    println("Start compiling HydrostaticFreeSurfaceModel")
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
        biogeochemistry,
    )
    println("Done compiling HydrostaticFreeSurfaceModel")

    ## Simulation
    Δt = 1second
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
    similarity_theory = sim_setup.similarity_theory_callable(; sim_setup.similarity_theory_args...)
    coupled_model = OceanSeaIceModel(ocean_sim, sea_ice; atmosphere, radiation, similarity_theory)
    println("Initialized coupled model")

    coupled_simulation = Simulation(coupled_model; Δt)
    # function atm_progress(sim)
    #     atmosphere_message = @sprintf("max(Ta) = %.1f", maximum(sim.model.atmosphere.tracers.T))
    #     @info atmosphere_message
    #     return nothing
    # end
    # add_callback!(coupled_simulation, atm_progress, IterationInterval(100))
    println("Initialized coupled simulation")
    return coupled_simulation
end

compute_sea_ice_ocean_fluxes!(::NoSeaIceModel) = nothing

function time_step!(coupled_model::NoSeaIceModel, Δt; callbacks = [], compute_tendencies = true)
    ocean = coupled_model.ocean

    # Be paranoid and update state at iteration 0
    coupled_model.clock.iteration == 0 && update_state!(coupled_model, callbacks)

    time_step!(ocean)

    tick!(coupled_model.clock, ocean.Δt) # An Ocean-only model advances with the ocean time-step!
    update_state!(coupled_model, callbacks; compute_tendencies)

    return nothing
end

function update_state!(coupled_model::NoSeaIceModel, callbacks = []; compute_tendencies = false)
    time = Time(coupled_model.clock.time)
    update_model_field_time_series!(coupled_model.atmosphere, time)
    compute_atmosphere_ocean_fluxes!(coupled_model)
    return nothing
end

end # module FjordsSim
