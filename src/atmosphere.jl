using Oceananigans
using ClimaOcean
using ClimaOcean.DataWrangling.JRA55:
    JRA55Metadata,
    JRA55_time_indices,
    JRA55_variable_names,
    TotallyInMemory,
    compute_bounding_indices,
    compute_bounding_nodes,
    infer_longitudinal_topology,
    download_dataset,
    short_name,
    location,
    metadata_path,
    native_times
using ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    ComponentInterfaces,
    default_gravitational_acceleration,
    TemperatureDependentAirViscosity,
    MomentumRoughnessLength,
    ScalarRoughnessLength,
    SimilarityScales

import ClimaOcean.DataWrangling.JRA55: JRA55FieldTimeSeries
import ClimaOcean.DataWrangling.JRA55: compute_bounding_indices

""" Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth """
@inline PAR⁰(x, y, t) =
    60 * (1 - cos((t + 15days) * 2π / 365days)) * (1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days)^2))) + 2

# ClimaOcean v0.5.4 fix for the custom longitude and latitude
# this is called from set! and uses grid to find the locations,
# which are 1 index more than necessary
function compute_bounding_indices(longitude::Nothing, latitude::Nothing, grid, LX, LY, λc, φc)
    λbounds = compute_bounding_nodes(longitude, grid, LX, λnodes)
    φbounds = compute_bounding_nodes(latitude, grid, LY, φnodes)

    i₁, i₂ = compute_bounding_indices(λbounds, λc)
    j₁, j₂ = compute_bounding_indices(φbounds, φc)
    TX = infer_longitudinal_topology(λbounds)

    # to prevent taking larger than grid areas
    i₁ = (i₂ - i₁ >= grid.Nx) ? (i₂ - grid.Nx + 1) : i₁
    j₁ = (j₂ - j₁ >= grid.Ny) ? (j₂ - grid.Ny + 1) : j₁

    return i₁, i₂, j₁, j₂, TX
end

# Julia does not dispatch on keyword argument, so this overwrites the original definition
# we need this to calculate the correct λc and φc values (ClimaOcean v0.5.4)
function JRA55FieldTimeSeries(
    metadata::JRA55Metadata,
    architecture = CPU(),
    FT = Float32;
    latitude = nothing,
    longitude = nothing,
    backend = InMemory(),
    time_indexing = Cyclical(),
    custom::Bool,
)

    # First thing: we download the dataset!
    download_dataset(metadata)

    # Unpack metadata details
    dataset = metadata.dataset
    name = metadata.name
    time_indices = JRA55_time_indices(dataset, metadata.dates, name)

    # Change the metadata to reflect the actual time indices
    dates = all_dates(dataset, name)[time_indices]
    metadata = Metadata(metadata.name, dates, metadata.dataset, metadata.dir)

    shortname = short_name(metadata)
    variable_name = metadata.name

    filepath = metadata_path(metadata) # Might be multiple paths!!!
    filepath = filepath isa AbstractArray ? first(filepath) : filepath

    # OnDisk backends do not support time interpolation!
    # Disallow OnDisk for JRA55 dataset loading
    if ((backend isa InMemory) && !isnothing(backend.length)) || backend isa OnDisk
        msg = string(
            "We cannot load the JRA55 dataset with a $(backend) backend. Use `InMemory()` or `JRA55NetCDFBackend(N)` instead.",
        )
        throw(ArgumentError(msg))
    end

    if !(variable_name ∈ JRA55_variable_names)
        variable_strs = Tuple("  - :$name \n" for name in JRA55_variable_names)
        variables_msg = prod(variable_strs)

        msg = string(
            "The variable :$variable_name is not provided by the JRA55-do dataset!",
            '\n',
            "The variables provided by the JRA55-do dataset are:",
            '\n',
            variables_msg,
        )

        throw(ArgumentError(msg))
    end

    # Record some important user decisions
    totally_in_memory = backend isa TotallyInMemory

    # Determine default time indices
    if totally_in_memory
        # In this case, the whole time series is in memory.
        # Either the time series is short, or we are doing a limited-area
        # simulation, like in a single column. So, we conservatively
        # set a default `time_indices = 1:2`.
        time_indices_in_memory = time_indices
        native_fts_architecture = architecture
    else
        # In this case, part or all of the time series will be stored in a file.
        # Note: if the user has provided a grid, we will have to preprocess the
        # .nc JRA55 data into a .jld2 file. In this case, `time_indices` refers
        # to the time_indices that we will preprocess;
        # by default we choose all of them. The architecture is only the
        # architecture used for preprocessing, which typically will be CPU()
        # even if we would like the final FieldTimeSeries on the GPU.
        time_indices_in_memory = 1:length(backend)
        native_fts_architecture = architecture
    end

    ds = Dataset(filepath)

    # Note that each file should have the variables
    #   - ds["time"]:     time coordinate
    #   - ds["lon"]:      longitude at the location of the variable
    #   - ds["lat"]:      latitude at the location of the variable
    #   - ds["lon_bnds"]: bounding longitudes between which variables are averaged
    #   - ds["lat_bnds"]: bounding latitudes between which variables are averaged
    #   - ds[shortname]: the variable data

    # Interfaces for the "native" JRA55 grid
    λn = ds["lon_bnds"][1, :]
    φn = ds["lat_bnds"][1, :]

    # The .nc coordinates lon_bnds and lat_bnds do not include
    # the last interface, so we push them here.
    push!(φn, 90)
    push!(λn, λn[1] + 360)

    # Nodes at the variable location
    # λc_from_file = ds["lon"][:]
    # φc_from_file = ds["lat"][:]  # these wont correspond to the oceananigans centers

    # ib₁, ib₂, jb₁, jb₂, TX = compute_bounding_indices(longitude, latitude, nothing, Face, Face, ds["lon"][:], ds["lat"][:])
    # create a grid to calculate the proper center nodes locations
    JRA55_entire_grid = LatitudeLongitudeGrid(
        CPU(),
        FT;
        halo = (3, 3),  # min.((ib₂ - ib₁ + 1, jb₂ - jb₁ + 1), (3, 3)),
        size = ((size(λn)...) - 1, (size(φn)...) - 1),
        longitude = λn,
        latitude = φn,
        topology = (Bounded, Bounded, Flat),
    )
    λc = λnodes(JRA55_entire_grid, Center())
    φc = φnodes(JRA55_entire_grid, Center())

    i₁, i₂, j₁, j₂, TX = compute_bounding_indices(longitude, latitude, nothing, Center, Center, λc, φc)
    data = ds[shortname][i₁:i₂, j₁:j₂, time_indices_in_memory]
    λr = λn[i₁:i₂+1]
    φr = φn[j₁:j₂+1]
    # λcr = λc[i₁:i₂]
    # φcr = φc[j₁:j₂]
    Nrx, Nry, Nt = size(data)
    close(ds)

    N = (Nrx, Nry)
    H = min.(N, (3, 3))

    JRA55_native_grid = LatitudeLongitudeGrid(
        native_fts_architecture,
        FT;
        halo = H,
        size = N,
        longitude = λr,
        latitude = φr,
        topology = (TX, Bounded, Flat),
    )

    boundary_conditions = FieldBoundaryConditions(JRA55_native_grid, (Center, Center, Nothing))
    times = native_times(metadata)

    if backend isa JRA55NetCDFBackend
        fts = FieldTimeSeries{Center,Center,Nothing}(
            JRA55_native_grid,
            times;
            backend,
            time_indexing,
            boundary_conditions,
            path = filepath,
            name = shortname,
        )

        # LX, LY, LZ = location(fts)
        # ii₁, ii₂, jj₁, jj₂, TX = compute_bounding_indices(nothing, nothing, fts.grid, LX, LY, λc, φc)
        # λc_bounded = λnodes(fts.grid, Center())
        # φc_bounded = φnodes(fts.grid, Center())
        # λc_bounded_face = λnodes(fts.grid, Face())
        # φc_bounded_face = φnodes(fts.grid, Face())

        # Fill the data in a GPU-friendly manner
        copyto!(interior(fts, :, :, 1, :), data)
        fill_halo_regions!(fts)

        return fts
    else
        native_fts = FieldTimeSeries{Center,Center,Nothing}(
            JRA55_native_grid,
            times;
            time_indexing,
            backend,
            boundary_conditions,
        )

        # Fill the data in a GPU-friendly manner
        copyto!(interior(native_fts, :, :, 1, :), data)
        fill_halo_regions!(native_fts)

        return native_fts
    end
end

function regional_roughness_lengths(FT = Oceananigans.defaults.FloatType)
    momentum = MomentumRoughnessLength(
        FT;
        gravitational_acceleration = default_gravitational_acceleration,
        maximum_roughness_length = 1.0, # An estimate?
        air_kinematic_viscosity = TemperatureDependentAirViscosity(FT),
        gravity_wave_parameter = 0.011,
        laminar_parameter = 0.11,
    )
    temperature = ScalarRoughnessLength(FT)
    water_vapor = ScalarRoughnessLength(FT)
    return SimilarityScales(momentum, temperature, water_vapor)
end
