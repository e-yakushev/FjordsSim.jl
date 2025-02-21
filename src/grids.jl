using Oceananigans.Grids: LatitudeLongitudeGrid, ImmersedBoundaryGrid
using Oceananigans.ImmersedBoundaries: GridFittedBottom
using ClimaOcean.VerticalGrids: exponential_z_faces
using JLD2: @load
using NCDatasets: Dataset

function compute_faces(centers)
    spacing = diff(centers)[1]  # Assuming uniform spacing
    faces = vcat([centers[1] - spacing / 2], (centers[1:end-1] .+ centers[2:end]) / 2, [centers[end] + spacing / 2])
    return faces
end

function grid_from_nc(arch, halo, filepath)
    ds = Dataset(filepath)
    z_faces = ds["z_faces"][:]
    # depths from an nc file are for grid centers
    depth = ds["h"][:, :]
    latitude = compute_faces(ds["lat"][:])
    longitude = compute_faces(ds["lon"][:])

    Nx, Ny = size(depth)
    Nz = length(z_faces)
    # Size should be for grid centers,
    # but z, latitude and langitude should be for faces
    underlying_grid =
        LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz - 1), halo = halo, z = z_faces, latitude, longitude)
    # depth should be for grid centers
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(coalesce.(depth, NaN)); active_cells_map = true)
    return grid
end

function grid_from_bathymetry_file(arch, halo, filepath, latitude, longitude)
    @load filepath depth z_faces
    Nx, Ny = size(depth)
    Nz = length(z_faces) - 1
    underlying_grid = LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz), halo = halo, z = z_faces, latitude, longitude)
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(depth); active_cells_map = true)
    return grid
end

function grid_latitude_flat!(arch, Nx, Ny, Nz, halo, latitude, longitude, depth)
    z_faces = exponential_z_faces(; Nz = Nz, depth = depth, h = depth)
    grid = LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz), halo = halo, z = z_faces, latitude, longitude)
    return grid
end

function grid_column!(arch, Nz, halo, latitude, longitude, depth, h)
    longitude = longitude .+ (-0.03, 0.03)
    latitude = latitude .+ (-0.03, 0.03)
    z_faces = exponential_z_faces(; Nz, depth, h)
    grid = LatitudeLongitudeGrid(arch; size = (3, 3, Nz), halo = halo, z = z_faces, latitude, longitude)
    return grid
end
