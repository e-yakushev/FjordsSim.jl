using Oceananigans
using Oceananigans.Architectures: GPU, CPU
using Oceananigans.Advection: WENO
using Oceananigans.BuoyancyFormulations: SeawaterBuoyancy, g_Earth
using Oceananigans.Coriolis: HydrostaticSphericalCoriolis, BetaPlane, Ω_Earth
using Oceananigans.TurbulenceClosures:
    TKEDissipationVerticalDiffusivity, ConvectiveAdjustmentVerticalDiffusivity, ScalarDiffusivity
using Oceananigans.OutputReaders: InMemory
using Oceananigans.Units: day
using ClimaOcean
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using FjordsSim:
    SetupModel,
    grid_from_nc,
    grid_ref,
    forcing_from_file,
    # regional_ocean_closure,
    bc_varna_bgh_oxydep,
    bgh_oxydep_boundary_conditions,
    bc_ocean,
    PAR⁰,
    # regional_roughness_lengths,
    free_surface_default,
    JRA55PrescribedAtmosphere,
    ComponentInterfaces,
    biogeochemistry_LOBSTER,
    biogeochemistry_OXYDEP,
    biogeochemistry_ref

const bottom_drag_coefficient = 0.003
const reference_density = 1020

args_oxydep = (
    initial_photosynthetic_slope = 0.1953 / day, # 1/(W/m²)/s
    Iopt = 80.0, # 50.0,     # (W/m2)
    alphaI = 1.8,   # [d-1/(W/m2)]
    betaI = 5.2e-4, # [d-1/(W/m2)]
    gammaD = 0.71,  # (-)
    Max_uptake = 1.7 / day,  # 1/d 2.0 4 5
    Knut = 1.5,            # (nd) 2.0
    r_phy_nut = 0.10 / day, # 1/d
    r_phy_pom = 0.15 / day, # 1/d
    r_phy_dom = 0.17 / day, # 1/d
    r_phy_het = 0.5 / day,  # 1/d 0.4 2.0
    Kphy = 0.1,             # (nd) 0.7
    r_pom_het = 0.7 / day,  # 1/d 0.7
    Kpom = 2.0,     # (nd)
    Uz = 0.6,       # (nd)
    Hz = 0.5,       # (nd)
    r_het_nut = 0.15 / day,      # 1/d 0.05
    r_het_pom = 0.15 / day,      # 1/d 0.02
    r_pom_nut_oxy = 0.006 / day, # 1/d
    r_pom_dom = 0.05 / day,      # 1/d
    r_dom_nut_oxy = 0.10 / day,  # 1/d
    O2_suboxic = 30.0,    # mmol/m3
    r_pom_nut_nut = 0.010 / day, # 1/d
    r_dom_nut_nut = 0.003 / day, # 1/d
    OtoN = 8.625, # (nd)
    CtoN = 6.625, # (nd)
    NtoN = 5.3,   # (nd)
    NtoB = 0.016, # (nd)
    sinking_speeds = (P = 0.15 / day, HET = 4.0 / day, POM = 10.0 / day),
)

function setup_region(;
    # Grid
    grid_callable = grid_from_nc,
    grid_args = (
        arch = GPU(),
        halo = (7, 7, 7),
        filepath = joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_105to232_bathymetry_v2.nc"),
    ),
    # Buoyancy
    buoyancy = SeawaterBuoyancy(; equation_of_state = TEOS10EquationOfState(; reference_density)),
    # Closure
    # closure = regional_ocean_closure(),
    # closure = TKEDissipationVerticalDiffusivity(),
    closure = (
        TKEDissipationVerticalDiffusivity(),
        Oceananigans.TurbulenceClosures.HorizontalScalarBiharmonicDiffusivity(ν = 15, κ = 10),
    ),
    # closure = ConvectiveAdjustmentVerticalDiffusivity(
    #     convective_κz = 5e-4, background_κz = 1e-5,
    #     convective_νz = 5e-1, background_νz = 1e-2,
    #     ),

    # Tracer advection
    tracer_advection = (T = WENO(), S = WENO(), e = nothing, ϵ = nothing),
    # Momentum advection
    momentum_advection = WENOVectorInvariant(),
    # Tracers
    tracers = (:T, :S, :e, :ϵ),
    initial_conditions = (T = 5.0, S = 33.0),
    # Free surface
    free_surface_callable = free_surface_default,
    free_surface_args = (grid_ref,),
    # Coriolis
    coriolis = HydrostaticSphericalCoriolis(rotation_rate = Ω_Earth),
    # Forcing
    forcing_callable = forcing_from_file,
    forcing_args = (
        grid_ref = grid_ref,
        filepath = joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_105to232_forcing_v2.nc"),
        tracers = tracers,
    ),
    # Boundary conditions
    bc_callable = bc_ocean,
    bc_args = (grid_ref, bottom_drag_coefficient),
    # Atmosphere
    atmosphere_callable = JRA55PrescribedAtmosphere,
    atmosphere_args = (arch = grid_args.arch, latitude = (58.98, 59.94), longitude = (10.18, 11.03)),
    # Ocean emissivity from https://link.springer.com/article/10.1007/BF02233853
    # With suspended matter 0.96 https://www.sciencedirect.com/science/article/abs/pii/0034425787900095
    radiation = ClimaOcean.Radiation(grid_args.arch; ocean_emissivity = 0.96),
    atmosphere_ocean_flux_formulation = SimilarityTheoryFluxes(),
    # Biogeochemistry
    biogeochemistry_callable = nothing,
    biogeochemistry_args = (nothing,),
    # Output folder
    results_dir = joinpath(homedir(), "FjordsSim_results", "oslofjord"),
)

    return SetupModel(
        grid_callable,
        grid_args,
        grid_ref,
        buoyancy,
        closure,
        tracer_advection,
        momentum_advection,
        tracers,
        initial_conditions,
        free_surface_callable,
        free_surface_args,
        coriolis,
        forcing_callable,
        forcing_args,
        bc_callable,
        bc_args,
        atmosphere_callable,
        atmosphere_args,
        radiation,
        atmosphere_ocean_flux_formulation,
        biogeochemistry_callable,
        biogeochemistry_args,
        biogeochemistry_ref;
        results_dir,
    )
end

setup_region_3d() = setup_region()
setup_region_3d_45to82() = setup_region(
    grid_args = (
        arch = GPU(),
        halo = (7, 7, 7),
        filepath = joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_45to82_bathymetry.nc"),
    ),
    tracers = (:T, :S, :e, :ϵ),
    forcing_args = (
        grid_ref = grid_ref,
        filepath = joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_45to82_forcing.nc"),
        tracers = (:T, :S, :e, :ϵ),
    ),
)
setup_region_3d_OXYDEP_45to82() = setup_region(
    grid_args = (
        arch = GPU(),
        halo = (7, 7, 7),
        filepath = joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_45to82_bathymetry.nc"),
    ),
    tracers = (:T, :S, :e, :C, :NUT, :P, :HET, :POM, :DOM, :O₂),
    forcing_args = (
        grid_ref = grid_ref,
        filepath = joinpath(homedir(), "FjordsSim_data", "oslofjord", "OF_inner_45to82_forcing.nc"),
        tracers = tracers,
    ),
    initial_conditions = (T = 5.0, S = 33.0, C = 0.0, NUT = 10.0, P = 0.05, HET = 0.01, O₂ = 350.0, DOM = 1.0),
    biogeochemistry_callable = biogeochemistry_OXYDEP,
    biogeochemistry_args = (grid_ref, args_oxydep),
    bc_callable = bc_varna_bgh_oxydep,
    bc_args = (grid_ref, bottom_drag_coefficient, biogeochemistry_ref),
    tracer_advection = (
        T = WENO(),
        S = WENO(),
        C = WENO(),
        e = nothing,
        NUT = WENO(),
        P = WENO(),
        HET = WENO(),
        POM = WENO(),
        DOM = WENO(),
        O₂ = WENO(),
    ),
)
setup_region_3d_OXYDEP() = setup_region(
    tracers = (:T, :S, :e, :C, :NUT, :P, :HET, :POM, :DOM, :O₂),
    initial_conditions = (T = 5.0, S = 33.0, C = 0.0, NUT = 10.0, P = 0.05, HET = 0.01, O₂ = 350.0, DOM = 1.0),
    biogeochemistry_callable = biogeochemistry_OXYDEP,
    biogeochemistry_args = (grid_ref, args_oxydep),
    bc_callable = bc_varna_bgh_oxydep,
    bc_args = (grid_ref, bottom_drag_coefficient, biogeochemistry_ref),
    tracer_advection = (
        T = WENO(),
        S = WENO(),
        C = WENO(),
        e = nothing,
        NUT = WENO(),
        P = WENO(),
        HET = WENO(),
        POM = WENO(),
        DOM = WENO(),
        O₂ = WENO(),
    ),
)
setup_region_3d_LOBSTER() = setup_region(
    tracers = (:T, :S, :e, :NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON, :DIC, :Alk, :O₂, :sPOC, :bPOC, :DOC),
    initial_conditions = (T = 5.0, S = 33.0, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0),
    biogeochemistry_callable = biogeochemistry_LOBSTER,
    biogeochemistry_args = (grid_ref,),
    bc_callable = bc_lobster,
    bc_args = (grid_ref, bottom_drag_coefficient),
    tracer_advection = (
        T = WENO(),
        S = WENO(),
        e = nothing,
        NO₃ = WENO(),
        NH₄ = WENO(),
        P = WENO(),
        Z = WENO(),
        sPON = WENO(),
        bPON = WENO(),
        DON = WENO(),
        DIC = WENO(),
        Alk = WENO(),
        O₂ = WENO(),
        sPOC = WENO(),
        bPOC = WENO(),
        DOC = WENO(),
    ),
)
