function get_lobster(; grid)
    biogeochemistry = LOBSTER(; grid, carbonates = false, open_bottom = false)
    # DIC_bcs = FieldBoundaryConditions(top = CarbonDioxideGasExchangeBoundaryCondition())
    return biogeochemistry  # , (; DIC = DIC_bcs)
end

function hydrophysics_simulation(
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
    boundary_conditions = ocean_boundary_conditions(grid, bottom_drag_coefficient),
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
    )

    ocean = Simulation(ocean_model; Δt, verbose)

    return ocean
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
    hydrophysics_boundary_conditions = ocean_boundary_conditions(grid, bottom_drag_coefficient)

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

    # biogeochemistry, bgh_boundary_conditions = get_lobster(; grid)
    biogeochemistry = get_lobster(; grid)
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
        # boundary_conditions = merge(hydrophysics_boundary_conditions, bgh_boundary_conditions),
        boundary_conditions = hydrophysics_boundary_conditions,
        biogeochemistry,
    )

    ocean = Simulation(ocean_model; Δt, verbose)

    return ocean
end
