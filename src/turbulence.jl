using Adapt
using Oceananigans.TurbulenceClosures.TKEBasedVerticalDiffusivities:
    CATKEVerticalDiffusivity, CATKEMixingLength, CATKEEquation, TKEDissipationDiffusivityFields

import Oceananigans.BoundaryConditions: fill_halo_regions!

Adapt.adapt_structure(to, tke_dissipation_diffusivity_fields::TKEDissipationDiffusivityFields) =
    TKEDissipationDiffusivityFields(adapt(to, tke_dissipation_diffusivity_fields.κu),
                                    adapt(to, tke_dissipation_diffusivity_fields.κc),
                                    adapt(to, tke_dissipation_diffusivity_fields.κe),
                                    adapt(to, tke_dissipation_diffusivity_fields.κϵ),
                                    adapt(to, tke_dissipation_diffusivity_fields.Le),
                                    adapt(to, tke_dissipation_diffusivity_fields.Lϵ),
                                    adapt(to, tke_dissipation_diffusivity_fields.previous_velocities),
                                    adapt(to, tke_dissipation_diffusivity_fields._tupled_tracer_diffusivities),
                                    adapt(to, tke_dissipation_diffusivity_fields._tupled_implicit_linear_coefficients))

function fill_halo_regions!(tke_dissipation_diffusivity_fields::TKEDissipationDiffusivityFields, args...; kw...)
    fields_with_halos_to_fill = (tke_dissipation_diffusivity_fields.κu,
                                 tke_dissipation_diffusivity_fields.κc,
                                 tke_dissipation_diffusivity_fields.κe,
                                 tke_dissipation_diffusivity_fields.κϵ)

    return fill_halo_regions!(fields_with_halos_to_fill, args...; kw...)
end

function regional_ocean_closure(FT = Oceananigans.defaults.FloatType)
    mixing_length = CATKEMixingLength(
        Cˢ   = 1.131,  # Surface distance coefficient for shear length scale
        Cᵇ   = 0.28,   # Bottom distance coefficient for shear length scale
        Cˢᵖ  = 0.505,  # Sheared convective plume coefficient
        CRiᵟ  = 1.02,   # Stability function width
        CRi⁰ = 0.254,  # Stability function lower Ri
        Cʰⁱu = 0.242,  # Shear mixing length coefficient for momentum at high Ri
        Cˡᵒu = 0.361,  # Shear mixing length coefficient for momentum at low Ri
        Cᵘⁿu = 0.370,  # Shear mixing length coefficient for momentum at negative Ri
        Cᶜu  = 3.705,  # Convective mixing length coefficient for tracers
        Cᵉu  = 0.0,    # Convective penetration mixing length coefficient for tracers
        Cʰⁱc = 0.098,  # Shear mixing length coefficient for tracers at high Ri
        Cˡᵒc = 0.369,  # Shear mixing length coefficient for tracers at low Ri
        Cᵘⁿc = 0.572,  # Shear mixing length coefficient for tracers at negative Ri
        Cᶜc  = 4.793,  # Convective mixing length coefficient for tracers
        Cᵉc  = 0.112,  # Convective penetration mixing length coefficient for tracers
        Cʰⁱe = 0.548,  # Shear mixing length coefficient for TKE at high Ri
        Cˡᵒe = 7.863,  # Shear mixing length coefficient for TKE at low Ri
        Cᵘⁿe = 1.447,  # Shear mixing length coefficient for TKE at negative Ri
        Cᶜe  = 3.642,  # Convective mixing length coefficient for TKE
        Cᵉe  = 0.0,    # Convective penetration mixing length coefficient for TKE
    )
    turbulent_kinetic_energy_equation = CATKEEquation(
        CʰⁱD = 0.579, # Dissipation length scale shear coefficient for high Ri
        CˡᵒD = 1.604, # Dissipation length scale shear coefficient for low Ri
        CᵘⁿD = 0.923, # Dissipation length scale shear coefficient for high Ri
        CᶜD  = 3.254, # Dissipation length scale convecting layer coefficient
        CᵉD  = 0.0,   # Dissipation length scale penetration layer coefficient
        Cᵂu★ = 3.179, # Surface shear-driven TKE flux coefficient
        CᵂwΔ = 0.383, # Surface convective TKE flux coefficient
        Cᵂϵ  = 1.0,   # Dissipative near-bottom TKE flux coefficient
    )

    return CATKEVerticalDiffusivity(
        VerticallyImplicitTimeDiscretization(),
        FT;
        mixing_length,
        turbulent_kinetic_energy_equation,
    )
end
