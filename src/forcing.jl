using Statistics
using Printf
using Dates: DateTime, Year, Second

using CUDA
using Oceananigans
using Oceananigans: Forcing
using Oceananigans.Grids
using Oceananigans.Units
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans.OutputReaders: Cyclical

LX = LY = LZ = Center

function native_times_to_seconds(native_times, start_time = native_times[1])
    times = []
    for native_time in native_times
        time = native_time - start_time
        time = Second(time).value
        push!(times, time)
    end

    return times
end

function field_time_series_tracer_forcing_func(i, j, k, grid, clock, fields, parameters)
    tr = @inbounds parameters.fts[i, j, k, Time(clock.time)]
    condition = !(tr < -990.0)
    radiation_term = -parameters.λOpen * (fields.T[i, j, k] - tr)
    return @inbounds ifelse(condition, radiation_term, 0)
end

function save_fts(; jld2_filepath, fts_name, fts, grid, times, boundary_conditions)
    isfile(jld2_filepath) && rm(jld2_filepath)
    on_disk_fts = FieldTimeSeries{LX,LY,LZ}(
        grid,
        times;
        boundary_conditions,
        backend = OnDisk(),
        path = jld2_filepath,
        name = fts_name,
    )
    for i = 1:size(fts)[end]
        set!(on_disk_fts, fts[i], i, times[i])
    end
end

function forcing_fields_from_file(grid_callable, grid_args, datapath)
    # grid = grid_callable(grid_args...)
    grid_args = (grid_args..., arch = CPU())
    grid_cpu = grid_callable(grid_args...)
    time_indexing = Cyclical()

    # Load data
    file = jldopen(joinpath(datapath, "Varna_BRY.jld2"))
    forcing_times = file["data_dict"]["time"]
    u = file["data_dict"]["uo"]
    v = file["data_dict"]["vo"]
    temp = file["data_dict"]["thetao"]
    salt = file["data_dict"]["so"]
    depth = file["data_dict"]["depth"]
    close(file)

    indices_1990 = findall(dt -> Year(dt) == Year(DateTime(1990)), forcing_times)

    temp_1990 = temp[:, indices_1990]

    year_1990_dates = forcing_times[indices_1990]
    times = native_times_to_seconds(year_1990_dates)

    bathymetry = Array(interior(grid_cpu.immersed_boundary.bottom_height, :, :, :))
    east_border = bathymetry[end, :, 1]
    east_border_water_indices = findall(x -> x < 0, east_border)

    boundary_conditions = FieldBoundaryConditions(grid_cpu, (LX, LY, LZ))
    fts = FieldTimeSeries{LX,LY,LZ}(grid_cpu, times; time_indexing, boundary_conditions)
    data = fill(-999.0, size(fts)...)
    for i = 1:size(data)[end]
        temp_mean = mean(temp_1990[:, i])
        data[end, :, :, i] .= temp_mean
    end
    copyto!(interior(fts, :, :, :, :), data)
    fill_halo_regions!(fts)

    fts_name = "temp"
    jld2_filepath = joinpath(datapath, string("Varna_forcing_repeat_year_", fts_name, ".jld2"))
    save_fts(; jld2_filepath, fts_name, fts, grid = grid_cpu, times, boundary_conditions)

    forcing_fts = FieldTimeSeries(
        jld2_filepath,
        fts_name;
        # architecture = GPU(),
        backend = InMemory(),
        time_indexing,
    )
    # for i = 1:size(forcing_fts)[end]
    #     tr = forcing_fts[1, 1, 1, Time(5i)]
    #     println("Indicing check.")
    # end
    # fill_halo_regions!(forcing_fts)

    λOpen = 1 / (1hours)  # Relaxation timescale [s⁻¹] Open boundary
    TForcing = Forcing(
        field_time_series_tracer_forcing_func,
        discrete_form = true,
        parameters = (fts = forcing_fts, λOpen = λOpen),
    )

    return (T = TForcing,)
end

# BOTTOM DRAG
function forcing_bottom_drag(bottom_drag_coefficient)
    Fu = Forcing(u_immersed_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    Fv = Forcing(v_immersed_bottom_drag, discrete_form = true, parameters = bottom_drag_coefficient)
    return (u = Fu, v = Fv)
end

# radiation QuasiOpenBoundary
function radiation_QuasiOpenBoundary_TS(Nx, external_values)
    c = 1 / (30minutes)  # Relaxation timescale [s⁻¹].

    Text = external_values.T
    Sext = external_values.S

    Tradiation(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse(i == Nx, -c * (model_fields.T[i, j, k] - Text), 0)
    FTradiation = Forcing(Tradiation, field_dependencies = :T, discrete_form = true)

    Sradiation(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse(i == Nx, -c * (model_fields.S[i, j, k] - Sext), 0)
    FSradiation = Forcing(Sradiation, field_dependencies = :S, discrete_form = true)

    return (T = FTradiation, S = FSradiation)
end

# RIVERS ONLY
function forcing_rivers_S(Nz)
    λ = 1 / (30minutes)  # Relaxation timescale [s⁻¹].

    S_source = 0.1

    source_index = (1, 13, Nz)

    S_point_source(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.S[i, j, k] - S_source), 0)
    Sforcing = Forcing(S_point_source, field_dependencies = :S, discrete_form = true)

    return (S = Sforcing,)
end

function forcing_rivers_NO₃(Nz)
    λ = 1 / (30minutes)  # Relaxation timescale [s⁻¹].

    S_source = 0.1
    NO₃_source = 10

    source_index = (1, 13, Nz)

    S_point_source(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.S[i, j, k] - S_source), 0)
    Sforcing = Forcing(S_point_source, field_dependencies = :S, discrete_form = true)

    NO₃_point_source(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.NO₃[i, j, k] - NO₃_source), 0)
    NO₃forcing = Forcing(NO₃_point_source, field_dependencies = :NO₃, discrete_form = true)

    return (S = Sforcing, NO₃ = NO₃forcing)
end

#### COMBINED FOR EACH TRACER ####
#### Combine several forcing functions, f.e. rivers + open boundary conditions ####
# All should have discrete form

# Temperature
function combined_forcing_func_T(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.T[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.T[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

# Salinity
function combined_forcing_func_S(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.S[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.S[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

# NUT
function combined_forcing_func_NUT(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.NUT[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.NUT[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

# O2
function combined_forcing_func_O2(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.O₂[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.O₂[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

# PHY
function combined_forcing_func_P(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.P[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.P[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

# HET
function combined_forcing_func_HET(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.HET[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.HET[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

# DOM
function combined_forcing_func_DOM(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.DOM[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.DOM[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

# Contaminant
function combined_forcing_func_C(i, j, k, grid, clock, model_fields, parameters)
    # First condition (from Tradiation)
    condition1 = (i == parameters.Nx)

    # Second condition (from T_point_source)
    condition2 = ((i, j, k) == parameters.src_loc)

    # Calculate the two possible forcing terms
    radiation_term = -parameters.λOpen * (model_fields.C[i, j, k] - parameters.external_value)
    point_source_term = -parameters.λRiver * (model_fields.C[i, j, k] - parameters.Vsrc)

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds ifelse(
        condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0,
    )
end

#### FINAL FORCING ####

function forcing_sognefjord(bottom_drag_coefficient, Nz, grid, external_values)
    Nx = grid[].Nx            # eastern open boundary index
    λRiver = 1 / (30minutes)  # Relaxation timescale [s⁻¹] River
    λOpen = 1 / (12hours)       # Relaxation timescale [s⁻¹] Open boundary

    # values in the river, CAN BE MOVED TO SETUP
    src_loc = (1, 13, Nz) # river  # (i, j, k)
    # src_loc = (111, 42, Nz)   # factory
    Tsrc = 10.0
    Ssrc = 0.1
    NUTsrc = 10.0
    DOMsrc = 5.0
    O2src = 300.0
    Psrc = 0.001
    HETsrc = 0.001
    Csrc = 100.0

    # # CUDA.@allowscalar user_z = Array(znodes(grid[], Center()))
    # # U, V, T, S, z, tᶠ = forcing_fields_from_file(user_z=user_z, Nz_model=Nz)

    # TForcing = Forcing(combined_forcing_func_T,
    #         field_dependencies = :T, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=Tsrc, src_loc=src_loc, external_value=external_values.T),
    #         discrete_form = true)
    # SForcing = Forcing(combined_forcing_func_S,
    #         field_dependencies = :S, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=Ssrc, src_loc=(1, 13, Nz), external_value=external_values.S),
    #         discrete_form = true)
    # NUTForcing = Forcing(combined_forcing_func_NUT,
    #         field_dependencies = :NUT, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=NUTsrc, src_loc=(1, 13, Nz), external_value=external_values.NUT),
    #         discrete_form = true)
    # DOMForcing = Forcing(combined_forcing_func_DOM,
    #         field_dependencies = :DOM, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=DOMsrc, src_loc=src_loc, external_value=external_values.DOM),
    #         discrete_form = true)
    # O2Forcing = Forcing(combined_forcing_func_O2,
    #         field_dependencies = :O₂, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=O2src, src_loc=src_loc, external_value=external_values.O₂),
    #         discrete_form = true)
    # PForcing = Forcing(combined_forcing_func_P,
    #         field_dependencies = :P, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=Psrc, src_loc=src_loc, external_value=external_values.P),
    #         discrete_form = true)
    # HETForcing = Forcing(combined_forcing_func_HET,
    #         field_dependencies = :HET, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=HETsrc, src_loc=src_loc, external_value=external_values.HET),
    #         discrete_form = true)

    # ContForcing = Forcing(combined_forcing_func_C,
    #         field_dependencies = :C, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=Csrc, src_loc=src_loc, external_value=external_values.C),
    #         discrete_form = true)

    # n_reference, tᶠ, F, λOpen, Nx = params
    # TForcing = Forcing(interp_forcing_T, field_dependencies = :T,
    #                     parameters=(n_reference = n, tᶠ = tᶠ, F=T, λOpen = λOpen, Nx = Nx),
    #                     discrete_form = true)
    forcing_bottom = forcing_bottom_drag(bottom_drag_coefficient)

    # forcing_rivers = forcing_rivers_S(Nz)
    final_forcing = forcing_bottom

    return final_forcing
end

function forcing_varna(bottom_drag_coefficient, Nz, grid, external_values)
    Nx = grid[].Nx            # eastern open boundary index
    λRiver = 1 / (30minutes)  # Relaxation timescale [s⁻¹] River
    λOpen = 1 / (12hours)       # Relaxation timescale [s⁻¹] Open boundary

    # values in the river, CAN BE MOVED TO SETUP
    src_loc = (1, 13, Nz) # river  # (i, j, k)
    # src_loc = (111, 42, Nz)   # factory
    Tsrc = 10.0
    Ssrc = 0.1
    NUTsrc = 10.0
    DOMsrc = 5.0
    O2src = 300.0
    Psrc = 0.001
    HETsrc = 0.001
    Csrc = 100.0

    TForcing = Forcing(
        combined_forcing_func_T,
        field_dependencies = :T,
        parameters = (
            Nx = Nx,
            λRiver = 0,
            λOpen = λOpen,
            Vsrc = Tsrc,
            src_loc = src_loc,
            external_value = external_values.T,
        ),
        discrete_form = true,
    )
    SForcing = Forcing(
        combined_forcing_func_S,
        field_dependencies = :S,
        parameters = (
            Nx = Nx,
            λRiver = λRiver,
            λOpen = λOpen,
            Vsrc = Ssrc,
            src_loc = (1, 13, Nz),
            external_value = external_values.S,
        ),
        discrete_form = true,
    )
    NUTForcing = Forcing(
        combined_forcing_func_NUT,
        field_dependencies = :NUT,
        parameters = (
            Nx = Nx,
            λRiver = λRiver,
            λOpen = λOpen,
            Vsrc = NUTsrc,
            src_loc = (1, 13, Nz),
            external_value = external_values.NUT,
        ),
        discrete_form = true,
    )
    DOMForcing = Forcing(
        combined_forcing_func_DOM,
        field_dependencies = :DOM,
        parameters = (
            Nx = Nx,
            λRiver = λRiver,
            λOpen = λOpen,
            Vsrc = DOMsrc,
            src_loc = src_loc,
            external_value = external_values.DOM,
        ),
        discrete_form = true,
    )
    O2Forcing = Forcing(
        combined_forcing_func_O2,
        field_dependencies = :O₂,
        parameters = (
            Nx = Nx,
            λRiver = 0,
            λOpen = λOpen,
            Vsrc = O2src,
            src_loc = src_loc,
            external_value = external_values.O₂,
        ),
        discrete_form = true,
    )
    PForcing = Forcing(
        combined_forcing_func_P,
        field_dependencies = :P,
        parameters = (
            Nx = Nx,
            λRiver = 0,
            λOpen = λOpen,
            Vsrc = Psrc,
            src_loc = src_loc,
            external_value = external_values.P,
        ),
        discrete_form = true,
    )
    HETForcing = Forcing(
        combined_forcing_func_HET,
        field_dependencies = :HET,
        parameters = (
            Nx = Nx,
            λRiver = 0,
            λOpen = λOpen,
            Vsrc = HETsrc,
            src_loc = src_loc,
            external_value = external_values.HET,
        ),
        discrete_form = true,
    )

    # -λRiver here, because the poin source is on the upper boundary and should point down
    ContForcing = Forcing(
        combined_forcing_func_C,
        field_dependencies = :C,
        parameters = (
            Nx = Nx,
            λRiver = λRiver,
            λOpen = λOpen,
            Vsrc = Csrc,
            src_loc = src_loc,
            external_value = external_values.C,
        ),
        discrete_form = true,
    )

    forcing_bottom = forcing_bottom_drag(bottom_drag_coefficient)

    final_forcing = merge(
        forcing_bottom,
        (
            T = TForcing,
            S = SForcing,
            C = ContForcing,
            NUT = NUTForcing,
            DOM = DOMForcing,
            O₂ = O2Forcing,
            P = PForcing,
            HET = HETForcing,
        ),
    )

    return final_forcing
end
