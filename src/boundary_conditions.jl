using CUDA
using Oceananigans.BoundaryConditions:
    FluxBoundaryCondition, ValueBoundaryCondition, FieldBoundaryConditions
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryCondition
using Oceananigans.Operators: Δzᵃᵃᶜ, ℑxyᶜᶠᵃ, ℑxyᶠᶜᵃ
using Oceananigans.Units
using Oceananigans.Architectures
using OceanBioME: GasExchange

const twelve_months = 12
const thirty_days = 30days

function wind_data_hardcoded(Nx, Ny)
    reference_density = 1000.0
    Cd = 0.0025
    ρₐᵢᵣ = 1.225
    Ntimes = 12
    uwind = [9, 3, 2, 1, 3, 1, 1, 2, 1, 1, 1, 3]
    vwind = [6, 9, 2, 3, 1, 1, 2, 2, 1, 2, 1, 2]
    τˣ = Array{Float64}(undef, Nx, Ny, Ntimes)
    τʸ = Array{Float64}(undef, Nx, Ny, Ntimes)
    for i = 1:Nx
        for j = 1:Ny
            τˣ[i, j, :] = -ρₐᵢᵣ * Cd .* (uwind .^ 2) ./ reference_density
            τʸ[i, j, :] = ρₐᵢᵣ * Cd .* (vwind .^ 2) ./ reference_density
        end
    end
    return τˣ, τʸ
end

function wind_data_hardcoded(::CPU, Nx, Ny)
    return wind_data_hardcoded(Nx, Ny)
end

function wind_data_hardcoded(::GPU, Nx, Ny)
    τˣ, τʸ = wind_data_hardcoded(Nx, Ny)
    return cu(τˣ), cu(τʸ)
end

current_time_index(time, tot_months) =
    mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
next_time_index(time, tot_months) =
    mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
cyclic_interpolate(u₁::Number, u₂, time) = u₁ + mod(time / thirty_days, 1) * (u₂ - u₁)

function surface_wind_stress(i, j, grid, clock, fields, τ)
    time = clock.time
    n₁ = current_time_index(time, twelve_months)
    n₂ = next_time_index(time, twelve_months)

    @inbounds begin
        τ₁ = τ[i, j, n₁]
        τ₂ = τ[i, j, n₂]
    end

    return cyclic_interpolate(τ₁, τ₂, time)
end

# Wind stress https://en.wikipedia.org/wiki/Wind_stress
function wind_stress(arch, Nx, Ny)
    τˣ, τʸ = wind_data_hardcoded(arch, Nx, Ny)

    u_wind_stress_bc =
        FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τˣ)
    v_wind_stress_bc =
        FluxBoundaryCondition(surface_wind_stress, discrete_form = true, parameters = τʸ)
    return u_wind_stress_bc, v_wind_stress_bc
end

# Bottom drag
function bottom_drag()
    # Linear bottom drag:
    μ = 0.003 # Non dimensional

    @inline speedᶠᶜᶜ(i, j, k, grid, fields) =
        @inbounds sqrt(fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)^2)
    @inline speedᶜᶠᶜ(i, j, k, grid, fields) =
        @inbounds sqrt(fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)^2)

    @inline u_bottom_drag(i, j, grid, clock, fields, μ) =
        @inbounds -μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
    @inline v_bottom_drag(i, j, grid, clock, fields, μ) =
        @inbounds -μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

    @inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) =
        @inbounds -μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields)
    @inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) =
        @inbounds -μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields)

    drag_u = FluxBoundaryCondition(u_immersed_bottom_drag, discrete_form = true, parameters = μ)
    drag_v = FluxBoundaryCondition(v_immersed_bottom_drag, discrete_form = true, parameters = μ)
    no_slip_bc = ValueBoundaryCondition(0)

    u_immersed_bc = ImmersedBoundaryCondition(
        bottom = drag_u,
        west = no_slip_bc,
        east = no_slip_bc,
        south = no_slip_bc,
        north = no_slip_bc,
    )

    v_immersed_bc = ImmersedBoundaryCondition(
        bottom = drag_v,
        west = no_slip_bc,
        east = no_slip_bc,
        south = no_slip_bc,
        north = no_slip_bc,
    )
    u_bottom_drag_bc = FluxBoundaryCondition(u_bottom_drag, discrete_form = true, parameters = μ)
    v_bottom_drag_bc = FluxBoundaryCondition(v_bottom_drag, discrete_form = true, parameters = μ)

    return u_immersed_bc, v_immersed_bc, u_bottom_drag_bc, v_bottom_drag_bc
end

function physics_boundary_conditions(arch, Nx, Ny)
    # u_wind_stress_bc, v_wind_stress_bc = wind_stress(arch, Nx, Ny)
    u_immersed_bc, v_immersed_bc, u_bottom_drag_bc, v_bottom_drag_bc = bottom_drag()

    u_bcs = FieldBoundaryConditions(
        # top = u_wind_stress_bc,
        bottom = u_bottom_drag_bc,
        immersed = u_immersed_bc,
    )

    v_bcs = FieldBoundaryConditions(
        # top = v_wind_stress_bc,
        bottom = v_bottom_drag_bc,
        immersed = v_immersed_bc,
    )
    return u_bcs, v_bcs
end

# BGC boundary conditions
function bgh_oxydep_boundary_conditions()
    O2_suboxic = 30.0  # OXY threshold for oxic/suboxic switch (mmol/m3)
    Trel = 10000.0     # Relaxation time for exchange with the sediments (s/m)
    b_ox = 15.0        # difference of OXY in the sediment and water, 
    b_NUT = 15.0       # NUT in the sediment, (mmol/m3)  
    b_DOM_ox = 10.0    # OM in the sediment (oxic conditions), (mmol/m3) 
    b_DOM_anox = 20.0  # OM in the sediment (anoxic conditions), (mmol/m3)  
    bu = 0.7           # Burial coeficient for lower boundary (0<Bu<1), 1 - for no burying, (nd)

    @inline F_ox(conc, threshold) = (0.5 + 0.5 * tanh(conc - threshold))
    @inline F_subox(conc, threshold) = (0.5 - 0.5 * tanh(conc - threshold))

    OXY_top = GasExchange(; gas = :O₂)
    @inline OXY_bottom_cond(i, j, grid, clock, fields) = @inbounds -(
        F_ox(fields.O₂[i, j, 1], O2_suboxic) * b_ox +
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * (0.0 - fields.O₂[i, j, 1])
    ) / Trel
    OXY_bottom = FluxBoundaryCondition(OXY_bottom_cond, discrete_form = true)

    @inline NUT_bottom_cond(i, j, grid, clock, fields) = @inbounds (
        F_ox(fields.O₂[i, j, 1], O2_suboxic) * (b_NUT - fields.NUT[i, j, 1]) +
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * (0.0 - fields.NUT[i, j, 1])
    ) / Trel
    NUT_bottom = FluxBoundaryCondition(NUT_bottom_cond, discrete_form = true) #ValueBoundaryCondition(10.0)

    w_PHY = biogeochemical_drift_velocity(biogeochemistry, Val(:PHY)).w[1, 1, 1]
    @inline PHY_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_PHY * fields.PHY[i, j, 1]
    PHY_bottom = FluxBoundaryCondition(PHY_bottom_cond, discrete_form = true)

    w_HET = biogeochemical_drift_velocity(biogeochemistry, Val(:HET)).w[1, 1, 1]
    @inline HET_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_HET * fields.HET[i, j, 1]
    HET_bottom = FluxBoundaryCondition(HET_bottom_cond, discrete_form = true)

    w_POM = biogeochemical_drift_velocity(biogeochemistry, Val(:POM)).w[1, 1, 1]
    @inline POM_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_POM * fields.POM[i, j, 1]
    POM_bottom = FluxBoundaryCondition(POM_bottom_cond, discrete_form = true)

    DOM_top = ValueBoundaryCondition(0.0)
    @inline DOM_bottom_cond(i, j, grid, clock, fields) = @inbounds (
        F_ox(fields.O₂[i, j, 1], O2_suboxic) * (b_DOM_ox - fields.DOM[i, j, 1]) +
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * 2.0 * (b_DOM_anox - fields.DOM[i, j, 1])
    ) / Trel
    DOM_bottom = FluxBoundaryCondition(DOM_bottom_cond, discrete_form = true) #, parameters = (; O2_suboxic, b_DOM_ox, Trel),)

    oxy_bcs = FieldBoundaryConditions(top = OXY_top, bottom = OXY_bottom)
    nut_bcs = FieldBoundaryConditions(bottom = NUT_bottom)
    dom_bcs = FieldBoundaryConditions(top = DOM_top, bottom = DOM_bottom)
    pom_bcs = FieldBoundaryConditions(bottom = POM_bottom)
    phy_bcs = FieldBoundaryConditions(bottom = PHY_bottom)
    het_bcs = FieldBoundaryConditions(bottom = HET_bottom)

    return oxy_bcs, nut_bcs, dom_bcs, pom_bcs, phy_bcs, het_bcs
end
