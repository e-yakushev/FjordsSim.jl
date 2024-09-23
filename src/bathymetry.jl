using ImageMorphology
using NCDatasets
using Downloads
using Oceananigans.Architectures: architecture
using Oceananigans.Grids: halo_size, λnodes, φnodes
using Oceananigans.Grids: x_domain, y_domain
using Oceananigans.Grids: topology
using Oceananigans.Utils: pretty_filesize, launch!
using Oceananigans.Fields: interpolate!
using Oceananigans.BoundaryConditions
using KernelAbstractions: @kernel, @index

"""
    regrid_bathymetry(target_grid;
                      url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf",
                      height_above_water = <none>,
                      minimum_depth = 0,
                      dir = joinpath(@__DIR__, "..", "data"),
                      filename = "ETOPO_2022_v1_60s_N90W180_surface.nc")

Regrid bathymetry associated with the NetCDF file at `path = joinpath(dir, filename)` to `target_grid`.
If `path` does not exist, then a download is attempted from `joinpath(url, filename)`.

Arguments:
==========

- target_grid: grid to interpolate onto

Keyword Arguments:
==================

- height_above_water: limits the maximum height of above-water topography (where h > 0). If
                      `nothing` the original topography is retained

- minimum_depth: minimum depth for the shallow regions. `h > minimum_depth` will be considered land

- dir: directory of the bathymetry-containing file

- filename: file containing bathymetric data. Must be netcdf with fields:
            (1) `lat` vector of latitude nodes
            (2) `lon` vector of longitude nodes
            (3) `z` matrix of depth values

- interpolation_passes: regridding/interpolation passes. The bathymetry is interpolated in
                        `interpolation_passes - 1` intermediate steps. With more steps the
                        final bathymetry will be smoother.
                        Example: interpolating from a 400x200 grid to a 100x100 grid in 4 passes will involve
                        - 400x200 -> 325x175
                        - 325x175 -> 250x150
                        - 250x150 -> 175x125
                        - 175x125 -> 100x100
                        If _coarsening_ the original grid, linear interpolation in passes is equivalent to
                        applying a smoothing filter, with more passes increasing the strength of the filter.
                        If _refining_ the original grid, additional passes will not help and no intermediate
                        steps will be performed.

- connected_regions_allowed: number of ``connected regions'' allowed in the bathymetry. Connected regions are fluid
                             regions that are fully encompassed by land (for example the ocean is one connected region).
                             Default is `Inf`. If a value < `Inf` is specified, connected regions will be preserved in order
                             of how many active cells they contain.
"""
function regrid_bathymetry_regional(
    target_grid;
    height_above_water = nothing,
    minimum_depth = 0,
    dir = joinpath(@__DIR__, "..", "data"),
    url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf",
    filename = "ETOPO_2022_v1_60s_N90W180_surface.nc",
    interpolation_passes = 1,
    connected_regions_allowed = 3,
) # Allow an `Inf` number of ``lakes''

    filepath = joinpath(dir, filename)

    if isfile(filepath)
        @info "Regridding bathymetry from existing file $filepath."
    else
        @info "Downloading bathymetry..."
        if !ispath(dir)
            @info "Making bathymetry directory $dir..."
            mkdir(dir)
        end

        fileurl = joinpath(url, filename)
        Downloads.download(fileurl, filepath; progress = download_progress, verbose = true)
    end

    dataset = Dataset(filepath)

    FT = eltype(target_grid)

    φ_data = dataset["lat"][:]
    λ_data = dataset["lon"][:]
    h_data = convert.(FT, dataset["z"][:, :])

    # Convert longitude to 0 - 360?
    # λ_data .+= 180
    # nhx      = size(h_data, 1)
    # h_data   = circshift(h_data, (nhx ÷ 2, 1))

    close(dataset)

    # Diagnose target grid information
    arch = architecture(target_grid)
    φ₁, φ₂ = y_domain(target_grid)
    λ₁, λ₂ = x_domain(target_grid)

    # Calculate limiting indices on the bathymetry grid
    i₁ = searchsortedfirst(λ_data, λ₁)
    i₂ = searchsortedfirst(λ_data, λ₂) - 1
    ii = i₁:i₂

    j₁ = searchsortedfirst(φ_data, φ₁)
    j₂ = searchsortedfirst(φ_data, φ₂) - 1
    jj = j₁:j₂

    # Remark if `target_grid` is not perfectly nested within the bathymetry grid
    Δλ = λ_data[2] - λ_data[1]
    Δφ = φ_data[2] - φ_data[1]

    λ₁_data = λ_data[i₁] - Δλ / 2
    λ₂_data = λ_data[i₂] + Δλ / 2
    φ₁_data = φ_data[j₁] - Δφ / 2
    φ₂_data = φ_data[j₂] + Δφ / 2

    λ₁ ≈ λ₁_data || @warn "The westernmost meridian of `target_grid` $λ₁ does not coincide with " *
          "the closest meridian of the bathymetry grid, $λ₁_data."
    λ₂ ≈ λ₂_data || @warn "The easternmost meridian of `target_grid` $λ₂ does not coincide with " *
          "the closest meridian of the bathymetry grid, $λ₂_data."
    φ₁ ≈ φ₁_data || @warn "The southernmost parallel of `target_grid` $φ₁ does not coincide with " *
          "the closest parallel of the bathymetry grid, $φ₁_data."
    φ₂ ≈ φ₂_data || @warn "The northernmost parallel of `target_grid` $φ₂ does not coincide with " *
          "the closest parallel of the bathymetry grid, $φ₂_data."

    # Restrict bathymetry _data to region of interest
    λ_data = λ_data[ii]
    φ_data = φ_data[jj]
    h_data = h_data[ii, jj]

    if !isnothing(height_above_water)
        # Overwrite the height of cells above water.
        # This has an impact on reconstruction. Greater height_above_water reduces total
        # wet area by biasing coastal regions to land during bathymetry regridding.
        land = h_data .> 0
        h_data[land] .= height_above_water
    end

    # Build the "native" grid of the bathymetry and make a bathymetry field.
    Nxn = length(λ_data)
    Nyn = length(φ_data)
    Nzn = 1

    native_grid = LatitudeLongitudeGrid(
        arch;
        size = (Nxn, Nyn, Nzn),
        latitude = (φ₁, φ₂),
        longitude = (λ₁, λ₂),
        z = (0, 1),
        halo = (10, 10, 1),
    )

    native_h = Field{Center,Center,Nothing}(native_grid)
    set!(native_h, h_data)

    target_h = interpolate_bathymetry_in_passes(
        native_h,
        target_grid;
        passes = interpolation_passes,
        connected_regions_allowed,
        minimum_depth,
    )

    return target_h
end
