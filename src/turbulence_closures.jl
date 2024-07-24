using Oceananigans.TurbulenceClosures:
    HorizontalScalarDiffusivity,
    VerticalScalarDiffusivity,
    VerticallyImplicitTimeDiscretization,
    ConvectiveAdjustmentVerticalDiffusivity,
    IsopycnalSkewSymmetricDiffusivity

function turbulence_closures_a()
    surface_νz = 1e-2
    background_νz = 1e-4
    background_κz = 1e-5

    @inline νz(x, y, z, t) = ifelse(z > -15, surface_νz, background_νz)

    horizontal_viscosity = HorizontalScalarDiffusivity(ν = 1e3)
    vertical_viscosity =
        VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = νz, κ = background_κz)
    convective_adjustment = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 0.1)

    κ_skew = 900.0      # [m² s⁻¹] skew diffusivity
    κ_symmetric = 900.0 # [m² s⁻¹] symmetric diffusivity

    gent_mcwilliams_diffusivity =
        IsopycnalSkewSymmetricDiffusivity(; κ_skew, κ_symmetric, slope_limiter = FluxTapering(1e-2))

    return (
        vertical_viscosity,
        horizontal_viscosity,
        convective_adjustment,
        gent_mcwilliams_diffusivity,
    )
end
