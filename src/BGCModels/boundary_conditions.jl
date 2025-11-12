using Oceananigans.BoundaryConditions: FluxBoundaryCondition, ValueBoundaryCondition, FieldBoundaryConditions

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

""" Function to calculate oxygen saturation in seawater """
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

    # Pressure correction factor (Weiss, 1970) for pressure in atm
    P_corr = 1.0 + P * (5.6e-6 + 2.0e-11 * P)

    # Adjusted oxygen saturation with pressure correction
    return (O2_sat * P_corr)
end

""" Sc, Schmidt number for O2  following Wanninkhof 2014 """
@inline function OxygenSchmidtNumber(T::Float64)::Float64
    return ((1920.4 - 135.6 * T + 5.2122 * T^2 - 0.10939 * T^3 + 0.00093777 * T^4))
    # can be replaced by PolynomialParameterisation{4}((a, b, c, d, e)) i.e.:
    #    a = 1953.4, b = - 128.0, c = 3.9918, d = -0.050091, e = 0.00093777  
    # Sc = PolynomialParameterisation{4}((a, b, c, d, e))
end

""" WindDependence, [mmol m-2s-1], Oxygen Sea Water Flux """
function WindDependence(windspeed::Float64)::Float64
    return (0.251 * windspeed^2.0) #ko2o=0.251*windspeed^2*(Sc/660)^(-0.5)  Wanninkhof 2014
end

""" OxygenSeaWaterFlux, [mmol m-2s-1], Oxygen Sea Water Flux """
function OxygenSeaWaterFlux(T::Float64, S::Float64, P::Float64, O₂::Float64, windspeed::Float64)::Float64
    return (
        WindDependence(windspeed) * (OxygenSchmidtNumber(T) / 660.0)^(-0.5) * (O₂ - oxygen_saturation(T, S, P)) * 0.24 /
        86400.0        # 0.24 is to convert from [cm/h] to [m/day]  * 0.24  / 86400.0
    )
end

# OXYDEP constants
const O2_suboxic = 10.0  # OXY threshold for oxic/suboxic switch (mmol/m3)
const Trel = 5000000.0 # 250000.0 OK    #25000.0 #10000.0     # Relaxation time for exchange with the sediments (s/m)
const b_ox = 15.0        # difference of OXY in the sediment and water, 
const b_NUT = 10.0       # NUT in the sediment, (mmol/m3)  
const b_DOM_ox = 6.0     # OM in the sediment (oxic conditions), (mmol/m3) 
const b_DOM_anox = 10.0   # OM in the sediment (anoxic conditions), (mmol/m3)  
const bu = 0.8  #0.85 0.6          # Burial coeficient for lower boundary (0<Bu<1), 1 - for no burying, (nd)
const windspeed = 5.0    # wind speed 10 m, (m/s)

""" BGC boundary conditions """
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
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * (fields.O₂[i, j, 1]- 0.0)
    ) / Trel
    OXY_bottom = FluxBoundaryCondition(OXY_bottom_cond, discrete_form = true)

    @inline NUT_bottom_cond(i, j, grid, clock, fields) = @inbounds (
        F_ox(fields.O₂[i, j, 1], O2_suboxic) * (b_NUT - fields.NUT[i, j, 1]) +
        F_subox(fields.O₂[i, j, 1], O2_suboxic) * (0.0 - fields.NUT[i, j, 1])
    ) / Trel
    NUT_bottom = FluxBoundaryCondition(NUT_bottom_cond, discrete_form = true)

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
    DOM_bottom = FluxBoundaryCondition(DOM_bottom_cond, discrete_form = true)

    oxy_bcs = FieldBoundaryConditions(top = OXY_top, bottom = OXY_bottom)
    nut_bcs = FieldBoundaryConditions(bottom = NUT_bottom)
    dom_bcs = FieldBoundaryConditions(top = DOM_top, bottom = DOM_bottom)
    pom_bcs = FieldBoundaryConditions(bottom = POM_bottom)
    phy_bcs = FieldBoundaryConditions(bottom = P_bottom)
    het_bcs = FieldBoundaryConditions(bottom = HET_bottom)

    bc_oxydep = (O₂ = oxy_bcs, NUT = nut_bcs, DOM = dom_bcs, POM = pom_bcs, P = phy_bcs, HET = het_bcs)

    return bc_oxydep
end