using Oceananigans
using Oceananigans.Units: minute, minutes, days, hour
using Printf

include("../../src/FjordsSim.jl")
include("setup.jl")

using .FjordsSim

setup = FjordsSetup(;oslo_fjord_setup...)
grid = ImmersedBoundaryGrid(setup)

# const surface_νz = 1e-2
# const background_νz = 1e-4
# const background_κz = 1e-5
#
# @inline νz(x, y, z, t) = ifelse(z > -49, surface_νz, background_νz)
#
# horizontal_viscosity = HorizontalScalarDiffusivity(ν=1e4)
# vertical_mixing      = RiBasedVerticalDiffusivity()
# vertical_viscosity   = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(),
#                                                  ν=νz, κ=background_κz)
#
# κ_skew = 9.0      # [m² s⁻¹] skew diffusivity
# κ_symmetric = 9.0 # [m² s⁻¹] symmetric diffusivity
#
# gent_mcwilliams_diffusivity = IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric,
#                                                                 slope_limiter = FluxTapering(1e-2))

# closures = (
#     vertical_viscosity,
#     horizontal_viscosity,
#     vertical_mixing,
#     gent_mcwilliams_diffusivity,
# )
closure = ScalarDiffusivity(ν=1e-4, κ=1e-4)

# free_surface = ImplicitFreeSurface()
# buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = 2e-4,
#                                                                     haline_contraction = 8e-4))
coriolis = HydrostaticSphericalCoriolis()

# model = HydrostaticFreeSurfaceModel(; grid = underlying_grid, free_surface, buoyancy, coriolis,
#                                     momentum_advection = VectorInvariant(),
#                                     tracer_advection = WENO(underlying_grid),
#                                     closure = closures,
#                                     tracers = (:T, :S))

dTdz = 1e-3 # K m⁻¹, temperature gradient

u₁₀ = 3 # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 1e-4 # dimensionless drag coefficient
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀); # m² s⁻²

Qʰ = 200.0  # W m⁻², surface _heat_ flux
cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater

Qᵀ = Qʰ / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux

T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
                                bottom = GradientBoundaryCondition(dTdz))
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

model = HydrostaticFreeSurfaceModel(; coriolis, closure, grid,
                                    momentum_advection = VectorInvariant(),
                                    tracer_advection = WENO(grid.underlying_grid),
                                    tracers = (:T, :S),
                                    boundary_conditions = (u=u_bcs, T=T_bcs)
                                    )

# Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise

# Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, y, z) = 20 + dTdz * z + dTdz * model.grid.Lz * 1e-3 * Ξ(z)

# Velocity initial condition: random noise scaled by the friction velocity.
uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z);

# `set!` the `model` fields using functions or constants:
set!(model, u=uᵢ, v=uᵢ, w=uᵢ, T=Tᵢ, S=35)

# set!(model, u=0.01, w=0.01, T=7, S=35)

Δt = 1minute
simulation = Simulation(model; Δt, stop_iteration=1440) #stop_time=Nyears * years)

start_time = [time_ns()]

function progress(sim)
    wall_time = (time_ns() - start_time[1]) * 1e-9

    u = sim.model.velocities.u
    w = sim.model.velocities.w

    intw  = Array(interior(w))
    max_w = findmax(intw)

    mw = max_w[1]
    iw = max_w[2]

    @info @sprintf("Time: % 12s, iteration: %d, max(|u|): %.2e ms⁻¹, wmax: %.2e , loc: (%d, %d, %d), wall time: %s",
                    prettytime(sim.model.clock.time),
                    sim.model.clock.iteration, maximum(abs, u), mw, iw[1], iw[2], iw[3],
                    prettytime(wall_time))

    start_time[1] = time_ns()

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S

output_prefix = joinpath(homedir(), "fjords_data", "oslo_fjord")
pickup = false
save_interval = 20minutes;

simulation.output_writers[:surface_fields] =
    JLD2OutputWriter(model, (; u, v, w, T, S),
                     schedule = TimeInterval(save_interval),
                     filename = output_prefix * "_snapshots",
                     # with_halos = true,
                     overwrite_existing = true)

# Let's goo!
@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

run!(simulation; pickup)

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""
