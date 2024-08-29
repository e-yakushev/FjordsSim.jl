"""
OXYgen DEPletion model, OXYDEP targests on the silmplest possible way of parameterization of the oxygen  (DO) fate in changeable redox conditions.
It has a simplified ecosystem, and simulates production of DO due to photosynthesis and consumation of DO for biota respiraion,
OM mineralization, nitrification, and oxidation of reduced specied of S, Mn, Fe, present in suboxic conditions.
For the details of  OxyDEP  implemented here see (Berezina et al, 2022)
Tracers
=======
OXYDEP consists of 6 state variables ( in N-units):
    Phy - all the phototrophic organisms (phytoplankton and bacteria).
    Phy grows due to photosynthesis, loses inorganic matter
    due to respiraion, and loses organic matter in dissolved (DOM) and particulate (POM)
    forms due to metabolism and mortality. Phy growth is limited by irradiance, temperature and NUT availability.
    Het - heterotrophs, can consume Phy and POM,  produce DOM and POM and respirate NUT.
    NUT - represents oxydized forms of nutrients (i.e. NO3 and NO2 for N),
    that doesn't need additional  oxygen for nitrification.
    DOM - is dissolved organic matter. DOM  includes all kinds of labile dissolved organic matter
    and reduced forms of inorganic nutrients (i.e. NH4 and Urea for N).
    POM - is particular organic matter (less labile than DOM). Temperature affects DOM and POM mineralization.
    Oxy - is dissolved oxygen.

Required submodels
==================
* Photosynthetically available radiation: PAR (W/m²)
"""

using Oceananigans: fields
using Oceananigans.Units
using Oceananigans.Fields: Field, TracerFields, CenterField, ZeroField
using Oceananigans.BoundaryConditions:
    fill_halo_regions!,
    ValueBoundaryCondition,
    FieldBoundaryConditions,
    regularize_field_boundary_conditions
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Architectures: architecture
using Oceananigans.Utils: launch!
using OceanBioME:
    setup_velocity_fields, show_sinking_velocities, Biogeochemistry, ScaleNegativeTracers
using OceanBioME.Light:
    update_TwoBandPhotosyntheticallyActiveRadiation!,
    default_surface_PAR,
    TwoBandPhotosyntheticallyActiveRadiation
using OceanBioME.Sediments: sinking_flux 
# using OceanBioME.Boundaries.Sediments: sinking_flux

import Adapt: adapt_structure, adapt
import Base: show, summary
import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity,
    update_biogeochemical_state!
import OceanBioME: redfield, conserved_tracers
import OceanBioME: maximum_sinking_velocity
import OceanBioME.Sediments: #import OceanBioME.Boundaries.Sediments:
    nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers

function update_biogeochemical_state!(model, PAR::TwoBandPhotosyntheticallyActiveRadiation)
    arch = architecture(model.grid)
    launch!(
        arch,
        model.grid,
        :xy,
        update_TwoBandPhotosyntheticallyActiveRadiation!,
        PAR.field,
        model.grid,
        model.tracers.PHY,
        PAR.surface_PAR,
        model.clock.time,
        PAR,
    )

    fill_halo_regions!(PAR.field, model.clock, fields(model))
end

struct OXYDEP{FT,B,W} <: AbstractContinuousFormBiogeochemistry
    # PHY
    initial_photosynthetic_slope::FT # α, 1/(W/m²)/s
    Iopt::FT   # Optimal irradiance (W/m2) =50 (Savchuk, 2002)
    alphaI::FT # initial slope of PI-curve [d-1/(W/m2)] (Wallhead?)
    betaI::FT  # photoinhibition parameter [d-1/(W/m2)] (Wallhead?)                
    gammaD::FT # adaptation to daylength parameter (-)    
    Max_uptake::FT # Maximum nutrient uptake rate d-1
    Knut::FT # Half-saturation constant for an uptake of NUT by PHY for the NUT/PHY ratio (nd) 
    r_phy_nut::FT # Specific respiration rate, (1/d)
    r_phy_pom::FT # Specific rate of Phy mortality, (1/d)
    r_phy_dom::FT # Specific rate of Phy excretion, (1/d)
    # HET
    r_phy_het::FT # Max.spec. rate of grazing of HET on PHY, (1/d)
    Kphy::FT # Half-sat.const.for grazing of HET on PHY for PHY/HET ratio (nd)
    r_pom_het::FT # Max.spec. rate of grazing of HET on POM, (1/d)
    Kpom::FT # Half-sat.const.for grazing of HET on POM for POM/HET ratio (nd)
    Uz::FT # Food absorbency for HET (nd)
    Hz::FT # Ratio between diss. and part. excretes of HET (nd)
    r_het_nut::FT # Specific HET respiration rate (1/d)
    r_het_pom::FT # Specific HET mortality rate (1/d)
    # POM
    r_pom_nut_oxy::FT # Specific rate of POM oxic decay, (1/d)
    r_pom_dom::FT # Specific rate of POM decomposition, (1/d)
    # DOM
    r_dom_nut_oxy::FT # Specific rate of DOM oxic decay, (1/d)
    # O₂
    O2_suboxic::FT    # O2 threshold for oxic/suboxic switch (mmol/m3)
    r_pom_nut_nut::FT # Specific rate of POM denitrification, (1/d)
    r_dom_nut_nut::FT # Specific rate of DOM denitrification, (1/d)
    OtoN::FT # Redfield (138/16) to NO3, (uM(O)/uM(N))
    CtoN::FT # Redfield (106/16) to NO3, (uM(C)/uM(N)) 
    NtoN::FT # Richards denitrification (84.8/16.), (uM(N)/uM(N))
    NtoB::FT # N[uM]/BIOMASS [mg/m3], (uM(N) / mgWW/m3)
    optionals::B
    sinking_velocities::W

    function OXYDEP(
        initial_photosynthetic_slope::FT,
        Iopt::FT,
        alphaI::FT,
        betaI::FT,
        gammaD::FT,
        Max_uptake::FT,
        Knut::FT,
        r_phy_nut::FT,
        r_phy_pom::FT,
        r_phy_dom::FT,
        r_phy_het::FT,
        Kphy::FT,
        r_pom_het::FT,
        Kpom::FT,
        Uz::FT,
        Hz::FT,
        r_het_nut::FT,
        r_het_pom::FT,
        r_pom_nut_oxy::FT,
        r_pom_dom::FT,
        r_dom_nut_oxy::FT,
        O2_suboxic::FT,
        r_pom_nut_nut::FT,
        r_dom_nut_nut::FT,
        OtoN::FT,
        CtoN::FT,
        NtoN::FT,
        NtoB::FT,
        optionals::B,
        sinking_velocities::W,
    ) where {FT,B,W}

        return new{FT,B,W}(
            initial_photosynthetic_slope,
            Iopt,
            alphaI,
            betaI,
            gammaD,
            Max_uptake,
            Knut,
            r_phy_nut,
            r_phy_pom,
            r_phy_dom,
            r_phy_het,
            Kphy,
            r_pom_het,
            Kpom,
            Uz,
            Hz,
            r_het_nut,
            r_het_pom,
            r_pom_nut_oxy,
            r_pom_dom,
            r_dom_nut_oxy,
            O2_suboxic,
            r_pom_nut_nut,
            r_dom_nut_nut,
            OtoN,
            CtoN,
            NtoN,
            NtoB,
            optionals,
            sinking_velocities,
        )
    end
end

function OXYDEP(;
    grid,
    initial_photosynthetic_slope::FT = 0.1953 / day, # 1/(W/m²)/s
    Iopt::FT = 50.0,     # (W/m2)
    alphaI::FT = 1.8,   # [d-1/(W/m2)]
    betaI::FT = 5.2e-4, # [d-1/(W/m2)]
    gammaD::FT = 0.71,  # (-)
    Max_uptake::FT = 2.5 / day,  # 1/d 2.0 4 5
    Knut::FT = 2.0,            # (nd)
    r_phy_nut::FT = 0.10 / day, # 1/d
    r_phy_pom::FT = 0.15 / day, # 1/d
    r_phy_dom::FT = 0.17 / day, # 1/d
    r_phy_het::FT = 2.0 / day,  # 1/d 0.4
    Kphy::FT = 0.1,             # (nd) 0.7
    r_pom_het::FT = 0.7 / day,  # 1/d 0.7
    Kpom::FT = 2.0,     # (nd)
    Uz::FT = 0.6,       # (nd)
    Hz::FT = 0.5,       # (nd)
    r_het_nut::FT = 0.15 / day,      # 1/d 0.05
    r_het_pom::FT = 0.10 / day,      # 1/d 0.02
    r_pom_nut_oxy::FT = 0.006 / day, # 1/d
    r_pom_dom::FT = 0.01 / day,      # 1/d
    r_dom_nut_oxy::FT = 0.050 / day,  # 1/d
    O2_suboxic::FT = 30.0,    # mmol/m3
    r_pom_nut_nut::FT = 0.010 / day, # 1/d
    r_dom_nut_nut::FT = 0.003 / day, # 1/d
    OtoN::FT = 8.625, # (nd)
    CtoN::FT = 6.625, # (nd)
    NtoN::FT = 5.3,   # (nd)
    NtoB::FT = 0.016, # (nd)
    surface_photosynthetically_active_radiation = default_surface_PAR,
    light_attenuation_model::LA = TwoBandPhotosyntheticallyActiveRadiation(;
        grid,
        surface_PAR = surface_photosynthetically_active_radiation,
    ),
    sediment_model::S = nothing,
    TS_forced::Bool = false,
    sinking_speeds = (PHY = 0.15 / day, HET = 0.4 / day, POM = 10.0 / day),
    open_bottom::Bool = true,
    scale_negatives = false,
    particles::P = nothing,
    modifiers::M = nothing,
) where {FT,LA,S,P,M}

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)
    optionals = Val(TS_forced)

    underlying_biogeochemistry = OXYDEP(
        initial_photosynthetic_slope,
        Iopt,
        alphaI,
        betaI,
        gammaD,
        Max_uptake,
        Knut,
        r_phy_nut,
        r_phy_pom,
        r_phy_dom,
        r_phy_het,
        Kphy,
        r_pom_het,
        Kpom,
        Uz,
        Hz,
        r_het_nut,
        r_het_pom,
        r_pom_nut_oxy,
        r_pom_dom,
        r_dom_nut_oxy,
        O2_suboxic,
        r_pom_nut_nut,
        r_dom_nut_nut,
        OtoN,
        CtoN,
        NtoN,
        NtoB,
        optionals,
        sinking_velocities,
    )

    if scale_negatives
        scaler = ScaleNegativeTracers(underlying_biogeochemistry, grid)
        modifiers = isnothing(modifiers) ? scaler : (modifiers..., scaler)
    end

    return Biogeochemistry(
        underlying_biogeochemistry;
        light_attenuation = light_attenuation_model,
        sediment = sediment_model,
        particles,
        modifiers,
    )
end

required_biogeochemical_tracers(::OXYDEP{<:Any,<:Val{false},<:Any}) =
    (:NUT, :PHY, :HET, :POM, :DOM, :O₂, :T)
required_biogeochemical_tracers(::OXYDEP{<:Any,<:Val{true},<:Any}) =
    (:NUT, :PHY, :HET, :POM, :DOM, :O₂)
required_biogeochemical_auxiliary_fields(::OXYDEP{<:Any,<:Val{false},<:Any}) = (:PAR,)
required_biogeochemical_auxiliary_fields(::OXYDEP{<:Any,<:Val{true},<:Any}) = (:PAR, :T)

@inline (bgc::OXYDEP)(tracer::Val{:NUT}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args...) = @inbounds bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args[end])
@inline (bgc::OXYDEP)(tracer::Val{:PHY}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args...) = @inbounds bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args[end])
@inline (bgc::OXYDEP)(tracer::Val{:HET}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args...) = @inbounds bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args[end])
@inline (bgc::OXYDEP)(tracer::Val{:POM}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args...) = @inbounds bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args[end])
@inline (bgc::OXYDEP)(tracer::Val{:DOM}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args...) = @inbounds bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args[end])
@inline (bgc::OXYDEP)(tracer::Val{:O₂}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args...)  = @inbounds bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, args[end])

#@inline (bgc::OXYDEP{<:Any, <:Val{(true, true, false)}, <:Any})(tracer::Val{:DIC}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, PAR) = bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, DIC, Alk, PAR)
@inline (bgc::OXYDEP{<:Any, <:Val{(true, true, false)}, <:Any})(tracer::Val{:Ch_free}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, PAR) = bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, Ch_free, Ch_PHY, Ch_HET, Ch_POM, Ch_DOM, PAR)
@inline (bgc::OXYDEP{<:Any, <:Val{(true, true, false)}, <:Any})(tracer::Val{:Ch_PHY}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, PAR)  = bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, Ch_free, Ch_PHY, Ch_HET, Ch_POM, Ch_DOM, PAR)
@inline (bgc::OXYDEP{<:Any, <:Val{(true, true, false)}, <:Any})(tracer::Val{:Ch_HET}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, PAR)  = bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, Ch_free, Ch_PHY, Ch_HET, Ch_POM, Ch_DOM, PAR)
@inline (bgc::OXYDEP{<:Any, <:Val{(true, true, false)}, <:Any})(tracer::Val{:Ch_POM}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, PAR)  = bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, Ch_free, Ch_PHY, Ch_HET, Ch_POM, Ch_DOM, PAR)
@inline (bgc::OXYDEP{<:Any, <:Val{(true, true, false)}, <:Any})(tracer::Val{:Ch_DOM}, x, y, z, t, NUT, PHY, HET, POM, DOM, O₂, PAR)  = bgc(tracer, x, y, z, t, NUT, PHY, HET, POM, DOM, Ch_free, Ch_PHY, Ch_HET, Ch_POM, Ch_DOM, PAR)

include("core.jl")

@inline function biogeochemical_drift_velocity(bgc::OXYDEP, ::Val{tracer_name}) where {tracer_name}
    if tracer_name in keys(bgc.sinking_velocities)
        return (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities[tracer_name])
    else
        return (u = ZeroField(), v = ZeroField(), w = ZeroField())
    end
end

summary(::OXYDEP{FT,NamedTuple{K,V}}) where {FT,K,V} = string("OXYDEP{$FT} model, with $K sinking")
show(io::IO, model::OXYDEP{FT}) where {FT} = print(
    io,
    string(
        "OXYDEP{$FT} model \n",
        "└── Sinking Velocities:",
        "\n",
        show_sinking_velocities(model.sinking_velocities),
    ),
)

@inline maximum_sinking_velocity(bgc::OXYDEP) = maximum(abs, bgc.sinking_velocities.POM.w)

adapt_structure(to, oxydep::OXYDEP) = OXYDEP(
    adapt(to, oxydep.initial_photosynthetic_slope),
    adapt(to, oxydep.Iopt),
    adapt(to, oxydep.alphaI),
    adapt(to, oxydep.betaI),
    adapt(to, oxydep.gammaD),
    adapt(to, oxydep.Max_uptake),
    adapt(to, oxydep.Knut),
    adapt(to, oxydep.r_phy_nut),
    adapt(to, oxydep.r_phy_pom),
    adapt(to, oxydep.r_phy_dom),
    adapt(to, oxydep.r_phy_het),
    adapt(to, oxydep.Kphy),
    adapt(to, oxydep.r_pom_het),
    adapt(to, oxydep.Kpom),
    adapt(to, oxydep.Uz),
    adapt(to, oxydep.Hz),
    adapt(to, oxydep.r_het_nut),
    adapt(to, oxydep.r_het_pom),
    adapt(to, oxydep.r_pom_nut_oxy),
    adapt(to, oxydep.r_pom_dom),
    adapt(to, oxydep.r_dom_nut_oxy),
    adapt(to, oxydep.O2_suboxic),
    adapt(to, oxydep.r_pom_nut_nut),
    adapt(to, oxydep.r_dom_nut_nut),
    adapt(to, oxydep.OtoN),
    adapt(to, oxydep.CtoN),
    adapt(to, oxydep.NtoN),
    adapt(to, oxydep.NtoB),
    adapt(to, oxydep.sinking_velocities),
)

@inline redfield(i, j, k, val_tracer_name, bgc::OXYDEP, tracers) = redfield(val_tracer_name, bgc)
@inline redfield(::Union{Val{:NUT}}, bgc::OXYDEP) = 0
@inline redfield(::Union{Val{:PHY},Val{:HET},Val{:POM},Val{:DOM}}, bgc::OXYDEP) = 6.56

@inline nitrogen_flux(i, j, k, grid, advection, bgc::OXYDEP, tracers) =
    sinking_flux(i, j, k, grid, advection, Val(:POM), bgc, tracers) +
    sinking_flux(i, j, k, grid, advection, Val(:PHY), bgc, tracers)

@inline carbon_flux(i, j, k, grid, advection, bgc::OXYDEP, tracers) =
    nitrogen_flux(i, j, k, grid, advection, bgc, tracers) * redfield(Val(:PHY), bgc)

@inline remineralisation_receiver(::OXYDEP) = :NUT

@inline conserved_tracers(::OXYDEP) = (:NUT, :PHY, :HET, :POM, :DOM, :O₂)
@inline sinking_tracers(bgc::OXYDEP) = keys(bgc.sinking_velocities)
