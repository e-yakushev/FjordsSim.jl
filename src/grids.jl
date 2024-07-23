using FileIO
using Oceananigans
using Oceananigans.Architectures
using JLD2

struct SetupGridRegrid
    arch::AbstractSerialArchitecture
    size::Tuple
    latitude::Tuple
    longitude::Tuple
    z
    halo::Tuple
    datadir::String
    filename::String
    height_above_water
    minimum_depth::Number
end

function SetupGridRegrid(;
    Nx::Integer,
    Ny::Integer,
    latitude::Tuple,
    longitude::Tuple,
    datadir:: String = joinpath(homedir(), "data_fjords"),
    filename::String,
    depth,
    surface_layer_Δz,
    stretching_factor,
    surface_layer_height,
    arch::AbstractSerialArchitecture = CPU(),
    halo::Tuple = (4, 4, 4),
    height_above_water = nothing,
    minimum_depth::Number = 0,
    )

    z_faces = stretched_vertical_faces(;
        depth,
        surface_layer_Δz,
        stretching = PowerLawStretching(stretching_factor),
        surface_layer_height
        )

    Nz = length(z_faces) - 1

    return SetupGridRegrid(
        arch,
        (Nx, Ny, Nz),
        latitude,
        longitude,
        z_faces,
        halo,
        datadir,
        filename,
        height_above_water,
        minimum_depth,
    )
end

struct SetupGridPredefinedFromFile
    arch::AbstractSerialArchitecture
    Nx::Integer
    Ny::Integer
    dx::Integer
    dy::Integer
    z_levels::Vector
    z_middle::Vector
    halo::Tuple
    datadir::String
    filename::String
end

function ImmersedBoundaryGrid(setup_grid::SetupGridRegrid)
    grid = LatitudeLongitudeGrid(setup_grid.arch;
                                 size = setup_grid.size,
                                 latitude = setup_grid.latitude,
                                 longitude = setup_grid.longitude,
                                 z = setup_grid.z,
                                 halo = setup_grid.halo)
    h = regrid_bathymetry(
        grid;
        height_above_water=setup_grid.height_above_water,
        minimum_depth=setup_grid.minimum_depth,
        dir=setup_grid.datadir,
        filename = setup_grid.filename,
    )

    return ImmersedBoundaryGrid(grid, GridFittedBottom(h))
end

function ImmersedBoundaryGrid(setup_grid::SetupGridPredefinedFromFile)
    Nz = length(setup_grid.z_levels) - 1
    underlying_grid = RectilinearGrid(
        setup_grid.arch,
        topology = (Bounded, Bounded, Bounded),
        size = (setup_grid.Nx, setup_grid.Ny, Nz),
        x = (0, setup_grid.dx * setup_grid.Nx),
        y = (0, setup_grid.dy * setup_grid.Ny),
        z = setup_grid.z_levels,
        halo = setup_grid.halo,
    )
    filepath_topo = joinpath(setup_grid.datadir, setup_grid.filename)
    @load filepath_topo depth
    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)
end

struct PowerLawStretching{T}
    power :: T
end

function (stretching::PowerLawStretching)(Δz, z)
    γ = stretching.power
    return Δz^γ
end

struct LinearStretching{T}
    coefficient :: T
end

function (stretching::LinearStretching)(Δz, z)
    c = stretching.coefficient
    return (1 + c) * Δz
end

"""
    stretched_vertical_faces(; surface_layer_Δz = 5.0,
                               surface_layer_height = 100.0,
                               constant_bottom_spacing_depth = Inf,
                               maximum_Δz = Inf,
                               stretching = PowerLawStretching(1.02),
                               rounding_digits = 1,
                               depth = 5000)

Return an array of cell interfaces with `surface_layer_Δz` spacing in
a surface layer of height `surface_layer_height`, and stretched according to
the function `stretching(Δz_above, z_above)` down to `depth`.
The interfaces extends from `Lz = -z[1]` to `0 = z[end]`, where `Lz ≥ depth`.

The grid spacing `Δz` is limited to be less than `maximum_Δz`.
The grid is also uniformly-spaced below `constant_bottom_spacing_depth`.

`rounding_digits` controls the accuracy with which the grid face positions are saved.
"""
function stretched_vertical_faces(; surface_layer_Δz = 5.0,
                                    surface_layer_height = 100.0,
                                    constant_bottom_spacing_depth = Inf,
                                    maximum_Δz = Inf,
                                    stretching = PowerLawStretching(1.02),
                                    rounding_digits = 1,
                                    depth = 5000)

    Δz₀ = surface_layer_Δz
    h₀ = surface_layer_height

    # Generate surface layer grid
    z = [-Δz₀ * (k-1) for k = 1:ceil(h₀ / Δz₀)]

    # Generate stretched interior grid
    Lz₀ = depth

    while z[end] > - Lz₀
        Δz_above = z[end-1] - z[end]

        if z[end] > - constant_bottom_spacing_depth
            Δz = stretching(Δz_above, z[end])
            Δz = min(maximum_Δz, Δz)
        else
            Δz = Δz_above
        end

        push!(z, round(z[end] - Δz, digits=rounding_digits))
    end

    # Reverse grid to be right-side-up
    z = reverse(z)

    return z
end

@inline exponential_profile(z; Lz, h) = (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) 

function exponential_z_faces(; Nz, depth, h = Nz / 4.5)

    z_faces = exponential_profile.((1:Nz+1); Lz = Nz, h)

    # Normalize
    z_faces .-= z_faces[1]
    z_faces .*= - depth / z_faces[end]
    
    z_faces[1] = 0.0

    return reverse(z_faces)
end
