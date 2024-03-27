## Packages and modules
using Oceananigans
using Oceananigans.Units: minute, minutes, days, hour
using Printf

include("../../src/FjordsSim.jl")
include("setup.jl")

using .FjordsSim

## Setup
setup = FjordsSetup(;oslo_fjord_setup...)

coriolis = HydrostaticSphericalCoriolis()
closure = ScalarDiffusivity(ν=1e-4, κ=1e-4)
grid = ImmersedBoundaryGrid(setup)
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(thermal_expansion = 2e-4,
                                                                    haline_contraction = 8e-4))

## Boundary conditions
Qʰ = 200.0  # W m⁻², surface _heat_ flux
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater
Qᵀ = Qʰ / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux

dTdz = 1e-2 # K m⁻¹, temperature gradient
T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
                                bottom = GradientBoundaryCondition(dTdz))

u₁₀ = 10 # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3 # dimensionless drag coefficient
ρₐ = 1.225  # kg m⁻³, average density of air at sea-level
Qᵘ = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀); # m² s⁻²
# u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))
vₛ = 1e-1
v_bcs = FieldBoundaryConditions(south = OpenBoundaryCondition(vₛ))

## Model
model = HydrostaticFreeSurfaceModel(; coriolis, closure, grid, buoyancy,
                                    momentum_advection = VectorInvariant(),
                                    tracer_advection = WENO(grid.underlying_grid),
                                    tracers = (:T, :S),
                                    # boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs)
                                    boundary_conditions = (v=v_bcs, T=T_bcs)
                                    )

## Set model, initial conditions
# Random noise damped at top and bottom
Ξ(z) = randn() * z / model.grid.Lz * (1 + z / model.grid.Lz) # noise
# Temperature initial condition: a stable density gradient with random noise superposed.
Tᵢ(x, y, z) = 20 + dTdz * z + dTdz * model.grid.Lz * 1e-3 * Ξ(z)
# Velocity initial condition: random noise scaled by the friction velocity.
uᵢ(x, y, z) = sqrt(abs(Qᵘ)) * 1e-3 * Ξ(z);
# `set!` the `model` fields using functions or constants:
set!(model, u=uᵢ, v=uᵢ, w=uᵢ, T=Tᵢ, S=35)

## Setup a simulation
Δt = 1minute
simulation = Simulation(model; Δt, stop_iteration=1440) #stop_time=Nyears * years)

## Setup output, callbacks
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
save_interval = 1hour;

simulation.output_writers[:surface_fields] =
    JLD2OutputWriter(model, (; u, v, w, T, S),
                     schedule = TimeInterval(save_interval),
                     filename = output_prefix * "_snapshots",
                     # with_halos = true,
                     overwrite_existing = true)

## Run simulation
@info "Running a simulation with Δt = $(prettytime(simulation.Δt))"

run!(simulation; pickup)

@info """
    Simulation took $(prettytime(simulation.run_wall_time))
    Free surface: $(typeof(model.free_surface).name.wrapper)
    Time step: $(prettytime(Δt))
"""
