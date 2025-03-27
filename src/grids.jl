using Oceananigans.Grids: LatitudeLongitudeGrid, ImmersedBoundaryGrid
using Oceananigans.ImmersedBoundaries: GridFittedBottom
using JLD2: @load
using NCDatasets: Dataset

function compute_faces(centers)
    spacing = diff(centers)[1]  # Assuming uniform spacing
    faces = vcat([centers[1] - spacing / 2], (centers[1:end-1] .+ centers[2:end]) / 2, [centers[end] + spacing / 2])
    return faces
end

"""
This grid maker uses a netcdf file to create a grid.
A netcdf file should have 4 variables:
- z_faces - 1d
- h - 2d
- lat and lon - 1d
"""
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
    bathymetry = Field{Center, Center, Nothing}(underlying_grid)
    set!(bathymetry, coalesce.(depth, 0.0))
    fill_halo_regions!(bathymetry)
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry); active_cells_map = true)
    return grid
end

"""Takes depths from a jld2 file and preset lats and lon, return a grid."""
function grid_from_bathymetry_file(arch, halo, filepath, latitude, longitude)
    @load filepath depth z_faces
    Nx, Ny = size(depth)
    Nz = length(z_faces) - 1
    underlying_grid = LatitudeLongitudeGrid(arch; size = (Nx, Ny, Nz), halo = halo, z = z_faces, latitude, longitude)
    bathymetry = Field{Center, Center, Nothing}(underlying_grid)
    set!(bathymetry, depth)
    fill_halo_regions!(bathymetry)
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry); active_cells_map = true)
    return grid
end
