using Oceananigans.Architectures
using Oceananigans.Units
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Coriolis: Ω_Earth
using Oceananigans.OutputReaders: InMemory
using ClimaOcean
using ClimaOcean.OceanSimulations:
    default_ocean_closure, default_momentum_advection, default_tracer_advection
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using OceanBioME

include("../../src/FjordsSim.jl")

using .FjordsSim: grid_from_bathymetry_file!, forcing_varna, bc_varna

free_surface_default(grid) = SplitExplicitFreeSurface(grid[]; cfl = 0.7)
atmosphere_JRA55(arch, backend, grid) = JRA55_prescribed_atmosphere(arch; backend, grid = grid[])
biogeochemistry_LOBSTER(grid) = LOBSTER(; grid = grid[], carbonates = false, open_bottom = false)

# Grid
Nz = 10
grid_callable! = grid_from_bathymetry_file!
grid_parameters = (
    arch = GPU(),
    Nz = Nz,
    halo = (7, 7, 7),
    datadir = joinpath(homedir(), "BadgerArtifacts"),
    filename = "Varna_topo_channels.jld2",
    latitude = (43.177, 43.214),
    longitude = (27.640, 27.947),
)
grid = Ref{Any}(nothing)

mutable struct SetupVarna
    grid_callable!::Function
    grid_parameters::NamedTuple
    grid::Ref
    buoyancy::Any
    closure::Any
    tracer_advection::Any
    momentum_advection::Any
    tracers::Tuple
    free_surface_callable::Function
    free_surface_args::Tuple
    coriolis::Any
    forcing_callable::Union{Nothing,Function}
    forcing_args::Union{Tuple,NamedTuple}
    bc_callable::Function
    bc_args::Tuple
    atmosphere_callable::Function
    atmosphere_args::NamedTuple
    radiation::Any
    biogeochemistry_callable::Union{Nothing,Function}
    biogeochemistry_args::Union{Tuple,NamedTuple}

    function SetupVarna(;
        bottom_drag_coefficient = 0.003,
        reference_density = 1020,
        # Buoyancy
        buoyancy = SeawaterBuoyancy(;
            gravitational_acceleration = g_Earth,
            equation_of_state = TEOS10EquationOfState(; reference_density),
        ),
        # Closure
        closure = default_ocean_closure(),
        # Tracer advection
        tracer_advection = (
            T = default_tracer_advection(),
            S = default_tracer_advection(),
            e = nothing,
        ),
        # Momentum advection
        momentum_advection = default_momentum_advection(),
        # Tracers
        tracers = (:T, :S, :e),
        # Free surface
        free_surface_callable = free_surface_default,
        free_surface_args = (grid,),
        # Coriolis
        coriolis = HydrostaticSphericalCoriolis(rotation_rate = Ω_Earth),
        # Forcing
        forcing_callable = forcing_varna,
        forcing_args = (bottom_drag_coefficient, Nz),
        # Boundary conditions
        bc_callable = bc_varna,
        bc_args = (grid, bottom_drag_coefficient),
        ## Atmosphere
        atmosphere_callable = atmosphere_JRA55,
        atmosphere_args = (arch = grid_parameters.arch, backend = InMemory(), grid = grid),
        radiation = Radiation(grid_parameters.arch),
        ## Biogeochemistry
        biogeochemistry_callable = nothing,
        biogeochemistry_args = (nothing,),
    )

        return new(
            grid_callable!,
            grid_parameters,
            grid,
            buoyancy,
            closure,
            tracer_advection,
            momentum_advection,
            tracers,
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
end

setup_varna_3d() = SetupVarna()
setup_varna_3d_Lobster() =
    SetupVarna(biogeochemistry_callable = biogeochemistry_LOBSTER, biogeochemistry_args = (grid,))

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
    sinking_speeds = (PHY = 0.15 / day, HET = 4. / day, POM = 10.0 / day),
)
