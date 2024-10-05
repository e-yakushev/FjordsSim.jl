using Oceananigans
using Oceananigans.Architectures: architecture
using Oceananigans.DistributedComputations: child_architecture
using Oceananigans.Grids: halo_size, λnodes, φnodes
using Oceananigans.Grids: x_domain, y_domain
using Oceananigans.Grids: topology
using Oceananigans.Utils: pretty_filesize, launch!
using Oceananigans.Fields: interpolate!
using Oceananigans.BoundaryConditions
using ClimaOcean
using ClimaOcean.Bathymetry: interpolate_bathymetry_in_passes, remove_minor_basins!

using ImageMorphology
using KernelAbstractions: @kernel, @index

using NCDatasets
using Downloads
using Scratch

download_bathymetry_cache::String = ""
function __init__()
    global download_bathymetry_cache = @get_scratch!("Bathymetry")
end

function regrid_bathymetry_regional(
    target_grid;
    height_above_water = nothing,
    minimum_depth = 0,
    dir = download_bathymetry_cache,
    url = "https://www.ngdc.noaa.gov/thredds/fileServer/global/ETOPO2022/60s/60s_surface_elev_netcdf",
    filename = "ETOPO_2022_v1_60s_N90W180_surface.nc",
    interpolation_passes = 1,
    major_basins = Inf,
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

        try
            Downloads.download(fileurl, filepath; progress = download_progress, verbose = true)
        catch
            cmd = `wget --no-check-certificate -O $filepath $fileurl`
            run(cmd)
        end
    end

    dataset = Dataset(filepath)

    FT = eltype(target_grid)

    φ_data = dataset["lat"][:]
    λ_data = dataset["lon"][:]
    z_data = convert(Array{FT}, dataset["z"][:, :])

    # Convert longitude from (-180, 180) to (0, 360)
    # λ_data .+= 180
    # Nhx = size(z_data, 1)
    # z_data = circshift(z_data, (Nhx ÷ 2, 0))

    close(dataset)

    # Diagnose target grid information
    arch = child_architecture(architecture(target_grid))
    φ₁, φ₂ = y_domain(target_grid)
    λ₁, λ₂ = x_domain(target_grid)

    if λ₁ < 0 || λ₂ > 360
        throw(ArgumentError("Cannot regrid bathymetry between λ₁ = $(λ₁) and λ₂ = $(λ₂).
                             Bathymetry data is defined on longitudes spanning λ = (0, 360)."))
    end

    # Calculate limiting indices on the bathymetry grid
    i₁ = searchsortedfirst(λ_data, λ₁)
    i₂ = searchsortedfirst(λ_data, λ₂) - 1
    ii = i₁:i₂

    j₁ = searchsortedfirst(φ_data, φ₁)
    j₂ = searchsortedfirst(φ_data, φ₂) - 1
    jj = j₁:j₂

    # Restrict bathymetry _data to region of interest
    λ_data = λ_data[ii]
    φ_data = φ_data[jj]
    z_data = z_data[ii, jj]

    if !isnothing(height_above_water)
        # Overwrite the height of cells above water.
        # This has an impact on reconstruction. Greater height_above_water reduces total
        # wet area by biasing coastal regions to land during bathymetry regridding.
        land = z_data .> 0
        z_data[land] .= height_above_water
    end

    # Infer the "native grid" of the bathymetry data and make a bathymetry field.
    Δλ = λ_data[2] - λ_data[1]
    Δφ = φ_data[2] - φ_data[1]

    λ₁_data = λ_data[1] - Δλ / 2
    λ₂_data = λ_data[end] + Δλ / 2
    φ₁_data = φ_data[1] - Δφ / 2
    φ₂_data = φ_data[end] + Δφ / 2

    Nxn = length(λ_data)
    Nyn = length(φ_data)
    Nzn = 1

    native_grid = LatitudeLongitudeGrid(
        arch;
        size = (Nxn, Nyn, Nzn),
        latitude = (φ₁_data, φ₂_data),
        longitude = (λ₁_data, λ₂_data),
        z = (0, 1),
        halo = (10, 10, 1),
    )

    native_z = Field{Center,Center,Nothing}(native_grid)
    set!(native_z, z_data)

    target_z = interpolate_bathymetry_in_passes(
        native_z,
        target_grid;
        passes = interpolation_passes,
        minimum_depth,
    )

    if minimum_depth > 0
        zi = interior(target_z, :, :, 1)

        # Set the height of cells with z > -mininum_depth to z=0.
        # (In-place + GPU-friendly)
        zi .*= zi .<= -minimum_depth
    end

    if major_basins < Inf
        remove_minor_basins!(target_z, major_basins)
    end

    fill_halo_regions!(target_z)

    return target_z
end
