using Oceananigans.BoundaryConditions:
    FluxBoundaryCondition, ValueBoundaryCondition, FieldBoundaryConditions
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryCondition
using Oceananigans.Operators: ℑxyᶜᶠᵃ, ℑxyᶠᶜᵃ
using Oceananigans.Units: days
using Oceananigans.Architectures: GPU, CPU
using Oceananigans.Grids: Center, Face
using ClimaOcean.OceanSimulations: u_quadratic_bottom_drag, v_quadratic_bottom_drag

import Oceananigans.Biogeochemistry: biogeochemical_drift_velocity

const twelve_months = 12
const thirty_days = 30days

# Coefficients from Garcia and Gordon (1992)
const A1 = -173.4292
const A2 = 249.6339
const A3 = 143.3483
const A4 = -21.8492
const A5 = -0.033096
const A6 = 0.014259
const B1 = -0.035274
const B2 = 0.001429
const B3 = -0.00007292
const C1 = 0.0000826

# Function to calculate oxygen saturation in seawater
function oxygen_saturation(T::Float64, S::Float64, P::Float64)::Float64

    T_kelvin = T + 273.15  # Convert temperature to Kelvin

    # Calculate the natural logarithm of oxygen saturation concentration
    ln_O2_sat =
        A1 +
        A2 * (100 / T_kelvin) +
        A3 * log(T_kelvin / 100) +
        A4 * T_kelvin / 100 +
        A5 * (T_kelvin / 100)^2 +
        A6 * (T_kelvin / 100)^3 +
        S * (B1 + B2 * (T_kelvin / 100) + B3 * (T_kelvin / 100)^2) +
        C1 * S^2

    # Oxygen saturation concentration in µmol/kg
    O2_sat = exp(ln_O2_sat) * 44.66

    # Pressure correction factor (Weiss, 1970)
    P_corr = 1.0 + P * (5.6e-6 + 2.0e-11 * P)

    # Adjusted oxygen saturation with pressure correction
    return (O2_sat * P_corr)
end

# Sc, Schmidt number for O2  following Wanninkhof 2014
@inline function OxygenSchmidtNumber(T::Float64)::Float64
    return ((1920.4 - 135.6 * T + 5.2122 * T^2 - 0.10939 * T^3 + 0.00093777 * T^4))
    # can be replaced by PolynomialParameterisation{4}((a, b, c, d, e)) i.e.:
    #    a = 1953.4, b = - 128.0, c = 3.9918, d = -0.050091, e = 0.00093777  
    # Sc = PolynomialParameterisation{4}((a, b, c, d, e))
end

# WindDependence, [mmol m-2s-1], Oxygen Sea Water Flux 
function WindDependence(windspeed::Float64)::Float64
    return (0.251 * windspeed^2.0) #ko2o=0.251*windspeed^2*(Sc/660)^(-0.5)  Wanninkhof 2014
end

# OxygenSeaWaterFlux, [mmol m-2s-1], Oxygen Sea Water Flux 
function OxygenSeaWaterFlux(
    T::Float64,
    S::Float64,
    P::Float64,
    O₂::Float64,
    windspeed::Float64,
)::Float64
    return (
        WindDependence(windspeed) *
        (OxygenSchmidtNumber(T) / 660.0)^(-0.5) *
        (O₂ - oxygen_saturation(T, S, P)) *
        0.24 / 86400.0        # 0.24 is to convert from [cm/h] to [m/day]  * 0.24  / 86400.0
    )
end

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

current_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days), tot_months) + 1
next_time_index(time, tot_months) = mod(unsafe_trunc(Int32, time / thirty_days) + 1, tot_months) + 1
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

# OXYDEP constants
const O2_suboxic = 30.0  # OXY threshold for oxic/suboxic switch (mmol/m3)
const Trel = 10000.0     # Relaxation time for exchange with the sediments (s/m)
const b_ox = 15.0        # difference of OXY in the sediment and water, 
const b_NUT = 10.0       # NUT in the sediment, (mmol/m3)  
const b_DOM_ox = 6.0     # OM in the sediment (oxic conditions), (mmol/m3) 
const b_DOM_anox = 20.0  # OM in the sediment (anoxic conditions), (mmol/m3)  
const bu = 0.85          # Burial coeficient for lower boundary (0<Bu<1), 1 - for no burying, (nd)
const windspeed = 5.0    # wind speed 10 m, (m/s)

# BGC boundary conditions
function bgh_oxydep_boundary_conditions(biogeochemistry, Nz)

    Oxy_top_cond(i, j, grid, clock, fields) = @inbounds (OxygenSeaWaterFlux(
        fields.T[i, j, Nz],
        fields.S[i, j, Nz],
        0.0,                # sea surface pressure
        fields.O₂[i, j, Nz],
        windspeed,
    ))

    OXY_top = FluxBoundaryCondition(Oxy_top_cond; discrete_form = true)

    # oxic - suboxic switches
    @inline F_ox(conc, threshold) = (0.5 + 0.5 * tanh(conc - threshold))
    @inline F_subox(conc, threshold) = (0.5 - 0.5 * tanh(conc - threshold))

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

    w_P(i, j) = biogeochemical_drift_velocity(biogeochemistry, Val(:P)).w[i, j, 1]
    @inline P_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_P(i, j) * fields.P[i, j, 1]
    P_bottom = FluxBoundaryCondition(P_bottom_cond, discrete_form = true)

    w_HET(i, j) = biogeochemical_drift_velocity(biogeochemistry, Val(:HET)).w[i, j, 1]
    @inline HET_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_HET(i, j) * fields.HET[i, j, 1]
    HET_bottom = FluxBoundaryCondition(HET_bottom_cond, discrete_form = true)

    w_POM(i, j) = biogeochemical_drift_velocity(biogeochemistry, Val(:POM)).w[i, j, 1]
    @inline POM_bottom_cond(i, j, grid, clock, fields) = @inbounds -bu * w_POM(i, j) * fields.POM[i, j, 1]
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
    phy_bcs = FieldBoundaryConditions(bottom = P_bottom)
    het_bcs = FieldBoundaryConditions(bottom = HET_bottom)

    bc_oxydep =
        (O₂ = oxy_bcs, NUT = nut_bcs, DOM = dom_bcs, POM = pom_bcs, P = phy_bcs, HET = het_bcs)

    return bc_oxydep
end

function bc_ocean(grid_ref, bottom_drag_coefficient)
    grid = grid_ref[]
    # Set up boundary conditions using Field
    top_zonal_momentum_flux = τx = Field{Face,Center,Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center,Face,Nothing}(grid)
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

    bc_ocean = (
        u = FieldBoundaryConditions(top = FluxBoundaryCondition(τx), bottom = u_bot_bc),
        v = FieldBoundaryConditions(top = FluxBoundaryCondition(τy), bottom = v_bot_bc),
        T = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ)),
        S = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ)),
    )
    return bc_ocean
end

function bc_varna(grid_ref, bottom_drag_coefficient)
    grid = grid_ref[]
    top_zonal_momentum_flux = τx = Field{Face,Center,Nothing}(grid)
    top_meridional_momentum_flux = τy = Field{Center,Face,Nothing}(grid)
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

    # 42 is length of y ax
    # sin_flux(y, z, t) = @inbounds ifelse(z == grid.Nz, 0.05 * sin((2 * 3.14 / 42) * y) , 0)
    # u_bc = BoundaryCondition(Open, sin_flux)

    bc = (
        u = FieldBoundaryConditions(top = FluxBoundaryCondition(τx), bottom = u_bot_bc),
        # u = FieldBoundaryConditions(top = FluxBoundaryCondition(τx), bottom = u_bot_bc, east = u_bc),
        v = FieldBoundaryConditions(top = FluxBoundaryCondition(τy), bottom = v_bot_bc),
        T = FieldBoundaryConditions(top = FluxBoundaryCondition(Jᵀ)),
        S = FieldBoundaryConditions(top = FluxBoundaryCondition(Jˢ)),
    )
    return bc
end

function bc_varna_bgh_oxydep(grid_ref, bottom_drag_coefficient, biogeochemistry_ref)
    Nz = grid_ref[].Nz
    bgc_model = biogeochemistry_ref[]
    bc_varna_tuple = bc_varna(grid_ref, bottom_drag_coefficient)
    bc_bgh_oxydep_tuple = bgh_oxydep_boundary_conditions(bgc_model, Nz)

    return merge(bc_varna_tuple, bc_bgh_oxydep_tuple)
end
