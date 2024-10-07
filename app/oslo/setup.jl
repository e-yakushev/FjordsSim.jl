using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Units
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.OutputReaders: InMemory
using ClimaOcean
using OceanBioME
using ClimaOcean.OceanSimulations:
    default_ocean_closure, default_momentum_advection, default_tracer_advection
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

include("../../src/FjordsSim.jl")

using .FjordsSim:
    SetupHydrostaticFreeSurface,
    grid_from_lat_lon!,
    grid,
    forcing_bottom_drag,
    bc_ocean,
    PAR⁰,
    free_surface_default,
    atmosphere_JRA55,
    biogeochemistry_LOBSTER,
    biogeochemistry_OXYDEP

function setup_oslo(;
    bottom_drag_coefficient = 0.003,
    reference_density = 1020,
    # Grid
    grid_callable! = grid_from_lat_lon!,
    grid_parameters = (
        arch = GPU(),
        Nx = 50,
        Ny = 50,
        halo = (7, 7, 7),
        latitude = (58.8, 59.8),
        longitude = (10, 11.25),
        datadir = joinpath(homedir(), "FjordsSim_data"),
        filename = "ETOPO_2022_v1_15s_N60E000_surface.nc",
        depth = 500,
        surface_layer_Δz = 3,  # it shouldn't be 1
        stretching_factor = 1.1,
        surface_layer_height = 15,
        height_above_water = 10,
        minimum_depth = 3,
    ),
    # Buoyancy
    buoyancy = SeawaterBuoyancy(;
        gravitational_acceleration = g_Earth,
        equation_of_state = TEOS10EquationOfState(; reference_density),
    ),
    # Closure
    closure = default_ocean_closure(),
    # Tracer advection
    tracer_advection = (
        T = default_tracer_advection(),
        S = default_tracer_advection(),
        e = nothing,
    ),
    # Momentum advection
    momentum_advection = default_momentum_advection(),
    # Tracers
    tracers = (:T, :S, :e),
    initial_conditions = (T = 10, S = 35),
    # Free surface
    free_surface_callable = free_surface_default,
    free_surface_args = (grid,),
    # Coriolis
    coriolis = HydrostaticSphericalCoriolis(rotation_rate = Ω_Earth),
    # Forcing
    forcing_callable = forcing_bottom_drag,
    forcing_args = (bottom_drag_coefficient),
    # Boundary conditions
    bc_callable = bc_ocean,
    bc_args = (grid, bottom_drag_coefficient),
    ## Atmosphere
    atmosphere_callable = atmosphere_JRA55,
    atmosphere_args = (arch = grid_parameters.arch, backend = InMemory(), grid = grid),
    radiation = Radiation(grid_parameters.arch),
    ## Biogeochemistry
    biogeochemistry_callable = nothing,
    biogeochemistry_args = (nothing,),
)

    return SetupHydrostaticFreeSurface(
        grid_callable!,
        grid_parameters,
        grid,
        buoyancy,
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
        biogeochemistry_callable,
        biogeochemistry_args,
    )
end
