using Oceananigans.Architectures
using Oceananigans.Units
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.OutputReaders: InMemory
using ClimaOcean
using ClimaOcean.OceanSimulations:
    default_ocean_closure, default_momentum_advection, default_tracer_advection
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using OceanBioME

include("../../src/FjordsSim.jl")

using .FjordsSim: grid_from_bathymetry_file!, forcing_varna, bc_varna

free_surface_default(grid) = SplitExplicitFreeSurface(grid[]; cfl = 0.7)
atmosphere_JRA55(arch, backend, grid) = JRA55_prescribed_atmosphere(arch; backend, grid = grid[])
biogeochemistry_LOBSTER(grid) = LOBSTER(; grid = grid[], carbonates = false, open_bottom = false)

mutable struct SetupVarna
    grid_callable!::Function
    grid_parameters::NamedTuple
    grid::Ref
    buoyancy::Any
    closure::Any
    tracer_advection::Any
    momentum_advection::Any
    tracers::Tuple
    free_surface_callable::Function
    free_surface_args::Tuple
    coriolis::Any
    forcing_callable::Function
    forcing_args::Tuple
    bc_callable::Function
    bc_args::Tuple
    atmosphere_callable::Function
    atmosphere_args::NamedTuple
    radiation::Any
    biogeochemistry_callable::Function
    biogeochemistry_args::Tuple

    function SetupVarna(Nz = 10, bottom_drag_coefficient = 0.003, reference_density = 1020)
        # Grid
        grid_callable! = grid_from_bathymetry_file!
        grid_parameters = (
            arch = CPU(),
            Nz = Nz,
            halo = (7, 7, 7),
            datadir = joinpath(homedir(), "BadgerArtifacts"),
            filename = "Varna_topo_channels.jld2",
            latitude = (43.177, 43.214),
            longitude = (27.640, 27.947),
        )
        grid = Ref{Any}(nothing)
        # Buoyancy
        buoyancy = SeawaterBuoyancy(;
            gravitational_acceleration = g_Earth,
            equation_of_state = TEOS10EquationOfState(; reference_density),
        )
        # Closure
        closure = default_ocean_closure()
        # Tracer advection
        tracer_advection = default_tracer_advection()
        tracer_advection = (T = tracer_advection, S = tracer_advection, e = nothing)
        # Momentum advection
        momentum_advection = default_momentum_advection()
        # Tracers
        tracers = (:T, :S)
        tracers = tuple(tracers..., :e)
        # Free surface
        free_surface_callable = free_surface_default
        free_surface_args = (grid,)
        # Coriolis
        coriolis = HydrostaticSphericalCoriolis(rotation_rate = Ω_Earth)
        # Forcing
        forcing_callable = forcing_varna
        forcing_args = (bottom_drag_coefficient, Nz)
        # Boundary conditions
        bc_callable = bc_varna
        bc_args = (grid, bottom_drag_coefficient)
        ## Atmosphere
        atmosphere_callable = atmosphere_JRA55
        atmosphere_args = (arch = grid_parameters.arch, backend = InMemory(), grid = grid)
        radiation = Radiation(grid_parameters.arch)
        ## Biogeochemistry
        biogeochemistry_callable = biogeochemistry_LOBSTER
        biogeochemistry_args = (grid,)

        return new(
            grid_callable!,
            grid_parameters,
            grid,
            buoyancy,
            closure,
            tracer_advection,
            momentum_advection,
            tracers,
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
end
