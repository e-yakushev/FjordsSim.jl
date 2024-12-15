using Oceananigans.Grids: LatitudeLongitudeGrid, ImmersedBoundaryGrid
using Oceananigans.ImmersedBoundaries: GridFittedBottom
using ClimaOcean.VerticalGrids: exponential_z_faces
using JLD2: @load

function grid_from_bathymetry_file(arch, halo, filepath, latitude, longitude)
    @load filepath depth z_faces
    Nx, Ny = size(depth)
    Nz = length(z_faces) - 1
    underlying_grid = LatitudeLongitudeGrid(
        arch;
        size = (Nx, Ny, Nz),
        halo = halo,
        z = z_faces,
        latitude,
        longitude,
    )
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)
    return grid
end

function grid_latitude_flat!(arch, Nx, Ny, Nz, halo, latitude, longitude, depth)
    z_faces = exponential_z_faces(; Nz = Nz, depth = depth, h = depth)
    grid = LatitudeLongitudeGrid(
        arch;
        size = (Nx, Ny, Nz),
        halo = halo,
        z = z_faces,
        latitude,
        longitude,
    )
    return grid
end

function grid_column!(arch, Nz, halo, latitude, longitude, depth, h)
    longitude = longitude .+ (-0.03, 0.03)
    latitude = latitude .+ (-0.03, 0.03)
    z_faces = exponential_z_faces(; Nz, depth, h)
    grid = LatitudeLongitudeGrid(
        arch;
        size = (3, 3, Nz),
        halo = halo,
        z = z_faces,
        latitude,
        longitude,
    )
    return grid
end
