using Oceananigans
using Oceananigans.Architectures
using Oceananigans.Units
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.OutputReaders: InMemory
using ClimaOcean
using OceanBioME
using ClimaOcean.OceanSimulations:
    default_ocean_closure, default_momentum_advection, default_tracer_advection
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Oceananigans.TurbulenceClosures: ConvectiveAdjustmentVerticalDiffusivity

include("../../src/FjordsSim.jl")

using .FjordsSim:
    SetupModel,
    grid_from_bathymetry_file!,
    grid_latitude_flat!,
    grid_column!,
    grid,
    forcing_varna,
    bc_varna_bgh_oxydep,
    bgh_oxydep_boundary_conditions,
    bc_varna,
    bc_ocean,
    PAR⁰,
    free_surface_default,
    atmosphere_JRA55,
    biogeochemistry_LOBSTER,
    biogeochemistry_OXYDEP

const bottom_drag_coefficient = 0.003
const reference_density = 1020

args_oxydep = (
    initial_photosynthetic_slope = 0.1953 / day, # 1/(W/m²)/s
    Iopt = 50.0,     # (W/m2)
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

# Values on the open boundary
external_values = (
    T = 10.0,
    S = 18.0,
    NUT = 2.0,
    O₂ = 300.0,
    P = 0.001,
    HET = 0.001,
)

function setup_varna(;
    # Grid
    grid_callable! = grid_from_bathymetry_file!,
    grid_parameters = (
        arch = GPU(),
        Nz = 10,
        halo = (7, 7, 7),
        datadir = joinpath(homedir(), "FjordsSim_data"),
        filename = "Varna_topo_channels.jld2",
        latitude = (43.177, 43.214),
        longitude = (27.640, 27.947),
    ),
    # Buoyancy
    buoyancy = SeawaterBuoyancy(;
        gravitational_acceleration = g_Earth,
        equation_of_state = TEOS10EquationOfState(; reference_density),
    ),
    # Closure
    closure = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 5e-5, background_κz=1e-5),

    # Tracer advection
    tracer_advection = (
        T = WENO(),
        S = WENO(),
        e = nothing,
    ),
    # Momentum advection
    momentum_advection = default_momentum_advection(),
    # Tracers
    tracers = (:T, :S, :e),
    initial_conditions = (T = 10, S = 15),
    # Free surface
    free_surface_callable = free_surface_default,
    free_surface_args = (grid,),
    # Coriolis
    coriolis = HydrostaticSphericalCoriolis(rotation_rate = Ω_Earth),
    # Forcing
    forcing_callable = forcing_varna,
    forcing_args = (bottom_drag_coefficient, grid_parameters.Nz, grid, external_values),
    # Boundary conditions
    bc_callable = bc_ocean,
    bc_args = (grid, bottom_drag_coefficient),
    ## Atmosphere
    atmosphere_callable = atmosphere_JRA55,
    atmosphere_args = (arch = grid_parameters.arch, backend = InMemory(), grid = grid),
    radiation = Radiation(grid_parameters.arch),
    ## Biogeochemistry
    biogeochemistry_callable = nothing,
    biogeochemistry_args = (nothing,),
)

    return SetupModel(
        grid_callable!,
        grid_parameters,
        grid,
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
        biogeochemistry_callable,
        biogeochemistry_args,
    )
end

setup_varna_3d() = setup_varna()
setup_varna_3d_Lobster() = setup_varna(
    biogeochemistry_callable = biogeochemistry_LOBSTER,
    biogeochemistry_args = (grid,),
    tracers = (:T, :S, :e, :NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM),
    initial_conditions = (T = 10, S = 15, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01),
    tracer_advection = (
        T = default_tracer_advection(),
        S = default_tracer_advection(),
        e = nothing,
        NO₃ = default_tracer_advection(),
        NH₄ = default_tracer_advection(),
        P = default_tracer_advection(),
        Z = default_tracer_advection(),
        sPOM = default_tracer_advection(),
        bPOM = default_tracer_advection(),
        DOM = default_tracer_advection(),
    ),
)
setup_varna_3d_OXYDEP() = setup_varna(
    tracers = (:T, :S, :e, :NUT, :P, :HET, :POM, :DOM, :O₂),
    initial_conditions = (T = 10, S = 15, NUT = 10.0, P = 0.01, HET = 0.05, O₂ = 350.0, DOM = 1.0),
    biogeochemistry_callable = biogeochemistry_OXYDEP,
    biogeochemistry_args = (grid, args_oxydep),
    bc_args = (grid, bottom_drag_coefficient, biogeochemistry_OXYDEP),
    bc_callable = bc_varna_bgh_oxydep,
    tracer_advection = (
        T = WENO(),
        S = WENO(),
        e = nothing,
        NUT = WENO(),
        P = WENO(),
        HET = WENO(),
        POM = WENO(),
        DOM = WENO(),
        O₂ = WENO(),
    ),
)
setup_varna_2d() = setup_varna(
    grid_callable! = grid_latitude_flat!,
    grid_parameters = (
        arch = GPU(),
        Nx = 30,
        Ny = 1,
        Nz = 20,
        halo = (1, 1, 1),
        latitude = (43.177, 43.214),
        longitude = (27.640, 27.947),
        depth = 20,
    ),
    closure = ScalarDiffusivity(ν = 1e-5, κ = 1e-5),
    tracer_advection = nothing,
    momentum_advection = nothing,
    tracers = (:T, :S),
    # Coriolis
    coriolis = nothing,
    # Forcing
    forcing_callable = NamedTuple,
    forcing_args = (),
    # Boundary conditions
    bc_callable = bc_ocean,
    bc_args = (grid, 0),
)
setup_varna_column() = setup_varna(
    grid_callable! = grid_column!,
    grid_parameters = (
        arch = GPU(),
        Nz = 20,
        halo = (3, 3, 3),
        latitude = 43.177,
        longitude = 27.640,
        depth = 20,
        h = 20,
    ),
    closure = ScalarDiffusivity(ν = 1e-5, κ = 1e-5),
    tracer_advection = nothing,
    momentum_advection = nothing,
    tracers = (:T, :S),
    # Coriolis
    coriolis = FPlane(latitude = 43.177),
    # Forcing
    forcing_callable = NamedTuple,
    forcing_args = (),
    # Boundary conditions
    bc_callable = bc_ocean,
    bc_args = (grid, 0),
)
