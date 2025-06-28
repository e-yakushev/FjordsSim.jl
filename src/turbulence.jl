using Oceananigans.Utils
using Oceananigans.TimeSteppers: implicit_step!
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    TKEDissipationDiffusivityFields,
    FlavorOfTD,
    compute_tke_dissipation_diffusivities!,
    compute_TKEDissipation_diffusivities!,
    substep_tke_dissipation!,
    get_time_step

import Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    top_dissipation_flux, compute_diffusivities!, time_step_tke_dissipation_equations!

@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy)

# is this needed?
@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any,<:Any}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy) +
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy)

# is this needed?
@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any,<:Any,<:Any}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy) +
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy) +
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[3], buoyancy)

# TKE dissipation vertical diffusivity with another closure 
@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:FlavorOfTD, <:Any}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[1], buoyancy)

@inline top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple::Tuple{<:Any, <:FlavorOfTD}, buoyancy) =
    top_dissipation_flux(i, j, grid, clock, fields, parameters, closure_tuple[2], buoyancy)

# overload to call the proper time_step_tke_dissipation_equations
function compute_diffusivities!(
    diffusivities::TKEDissipationDiffusivityFields,
    closure::FlavorOfTD,
    model;
    parameters = :xyz,
)

    arch = model.architecture
    grid = model.grid
    velocities = model.velocities
    tracers = model.tracers
    buoyancy = model.buoyancy
    clock = model.clock
    top_tracer_bcs = NamedTuple(c => tracers[c].boundary_conditions.top for c in propertynames(tracers))

    if isfinite(model.clock.last_Δt) # Check that we have taken a valid time-step first.
        # Compute e at the current time:
        #   * update tendency Gⁿ using current and previous velocity field
        #   * use tridiagonal solve to take an implicit step
        time_step_tke_dissipation_equations!(model, diffusivities, closure)
    end

    # Update "previous velocities"
    u, v, w = model.velocities
    u⁻, v⁻ = diffusivities.previous_velocities
    parent(u⁻) .= parent(u)
    parent(v⁻) .= parent(v)

    launch!(
        arch,
        grid,
        parameters,
        compute_TKEDissipation_diffusivities!,
        diffusivities,
        grid,
        closure,
        velocities,
        tracers,
        buoyancy,
    )

    return nothing
end

# pass proper diffusivities and closure
function time_step_tke_dissipation_equations!(model, diffusivities::TKEDissipationDiffusivityFields, closure::FlavorOfTD)

    e = model.tracers.e
    ϵ = model.tracers.ϵ
    arch = model.architecture
    grid = model.grid
    Gⁿe = model.timestepper.Gⁿ.e
    G⁻e = model.timestepper.G⁻.e
    Gⁿϵ = model.timestepper.Gⁿ.ϵ
    G⁻ϵ = model.timestepper.G⁻.ϵ

    κe = diffusivities.κe
    κϵ = diffusivities.κϵ
    Le = diffusivities.Le
    Lϵ = diffusivities.Lϵ
    previous_velocities = diffusivities.previous_velocities
    e_index = findfirst(k -> k == :e, keys(model.tracers))
    ϵ_index = findfirst(k -> k == :ϵ, keys(model.tracers))
    implicit_solver = model.timestepper.implicit_solver

    Δt = model.clock.last_Δt
    Δτ = get_time_step(closure)

    if isnothing(Δτ)
        Δτ = Δt
        M = 1
    else
        M = ceil(Int, Δt / Δτ) # number of substeps
        Δτ = Δt / M
    end

    FT = eltype(grid)

    for m = 1:M # substep
        if m == 1 && M != 1
            χ = convert(FT, -0.5) # Euler step for the first substep
        else
            χ = model.timestepper.χ
        end

        launch!(
            arch,
            grid,
            :xyz,
            compute_tke_dissipation_diffusivities!,
            κe,
            κϵ,
            grid,
            closure,
            model.velocities,
            model.tracers,
            model.buoyancy,
        )

        # Compute the linear implicit component of the RHS (diffusivities, L)
        # and step forward
        launch!(
            arch,
            grid,
            :xyz,
            substep_tke_dissipation!,
            Le,
            Lϵ,
            grid,
            closure,
            model.velocities,
            previous_velocities, # try this soon: model.velocities, model.velocities,
            model.tracers,
            model.buoyancy,
            diffusivities,
            Δτ,
            χ,
            Gⁿe,
            G⁻e,
            Gⁿϵ,
            G⁻ϵ,
        )

        implicit_step!(e, implicit_solver, closure, diffusivities, Val(e_index), model.clock, Δτ)

        implicit_step!(ϵ, implicit_solver, closure, diffusivities, Val(ϵ_index), model.clock, Δτ)
    end

    return nothing
end
