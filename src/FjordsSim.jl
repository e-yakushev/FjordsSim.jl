module FjordsSim

using Oceananigans.Architectures
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using OceanBioME
using ClimaOcean.OceanSimulations:
    default_free_surface,
    default_ocean_closure,
    default_momentum_advection,
    default_tracer_advection,
    u_quadratic_bottom_drag,
    v_quadratic_bottom_drag,
    u_immersed_bottom_drag,
    v_immersed_bottom_drag
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState

import Oceananigans.Architectures: on_architecture

include("utils.jl")
include("bathymetry.jl")
include("grids.jl")
include("initial_conditions.jl")
include("turbulence_closures.jl")
include("boundary_conditions.jl")
include("forcing.jl")

# include("BGCModels/BGCModels.jl")

# using .BGCModels: OXYDEP

# there is no a steprangelen method in oceananigans
# but adding it here is type piracy
# we need this when loading atmospheric forcing to the video memory
function on_architecture(::GPU, a::StepRangeLen)
    on_architecture(GPU(), collect(a))
end

function get_lobster(; grid)
    biogeochemistry = LOBSTER(; grid, carbonates = true, open_bottom = true)
    DIC_bcs = FieldBoundaryConditions(top = CarbonDioxideGasExchangeBoundaryCondition())
    return biogeochemistry, (; DIC = DIC_bcs)
end

function biogeochemical_simulation(
    grid;
    Δt = 5minutes,
    closure = default_ocean_closure(),
    free_surface = default_free_surface(grid),
    reference_density = 1020,
    rotation_rate = Ω_Earth,
    gravitational_acceleration = g_Earth,
    bottom_drag_coefficient = 0.003,
    forcing = NamedTuple(),
    coriolis = HydrostaticSphericalCoriolis(; rotation_rate),
    momentum_advection = default_momentum_advection(),
    tracer_advection = default_tracer_advection(),
    verbose = false,
)

    # Set up boundary conditions using Field
    top_zonal_momentum_flux = Jᵘ = Field{Face,Center,Nothing}(grid)
    top_meridional_momentum_flux = Jᵛ = Field{Center,Face,Nothing}(grid)
    top_ocean_heat_flux = Jᵀ = Field{Center,Center,Nothing}(grid)
    top_salt_flux = Jˢ = Field{Center,Center,Nothing}(grid)

    u_bot_bc = FluxBoundaryCondition(
        u_quadratic_bottom_drag,
        discrete_form = true,
        parameters = bottom_drag_coefficient,
    )
    v_bot_bc = FluxBoundaryCondition(
        v_quadratic_bottom_drag,
        discrete_form = true,
        parameters = bottom_drag_coefficient,
    )

    ocean_boundary_conditions = (
        u = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵘ), bottom = u_bot_bc),
        v = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵛ), bottom = v_bot_bc),
        T = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ)),
        S = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ)),
    )

    if grid isa ImmersedBoundaryGrid
        Fu = Forcing(
            u_immersed_bottom_drag,
            discrete_form = true,
            parameters = bottom_drag_coefficient,
        )
        Fv = Forcing(
            v_immersed_bottom_drag,
            discrete_form = true,
            parameters = bottom_drag_coefficient,
        )
        forcing = merge(forcing, (; u = Fu, v = Fv))
    end

    # Use the TEOS10 equation of state
    teos10 = TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(; gravitational_acceleration, equation_of_state = teos10)

    # Minor simplifications for single column grids
    Nx, Ny, _ = size(grid)
    if Nx == Ny == 1 # single column grid
        tracer_advection = nothing
        momentum_advection = nothing
    end

    tracers = (:T, :S)
    if closure isa CATKEVerticalDiffusivity
        tracers = tuple(tracers..., :e)
        tracer_advection = (; T = tracer_advection, S = tracer_advection, e = nothing)
    end

    biogeochemistry, bgh_boundary_conditions = get_lobster(; grid)
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
        boundary_conditions = merge(ocean_boundary_conditions, bgh_boundary_conditions),
        biogeochemistry,
    )

    ocean = Simulation(ocean_model; Δt, verbose)

    return ocean
end

end # module FjordsSim
