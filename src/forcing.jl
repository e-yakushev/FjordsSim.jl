using Oceananigans: Forcing
using Oceananigans.Units
using Oceananigans.Architectures: architecture, on_architecture
using Oceananigans
using Oceananigans.Grids
using Printf
using CUDA


# THIS doesnt work
# # Callback that updates the forcing time index.
# # This callback runs on the CPU. # And GPU?
# forcing_time_index = n = Ref(1)
# function update_time_index(sim)
#     n = findfirst(t -> t > time(sim), forcing_times)

#     if isnothing(n)
#         msg = @sprintf("Simulation time %.2e was not found in forcing times.", time(sim))
#         error(msg)
#     end

#     forcing_time_index[] = n
#     return nothing
# end

#### READ FROM FILE AND INTERPOLATE ####
# from https://github.com/glwagner/TropicalTurbulantics.jl/blob/main/tropical_turbulence_setup.jl
@inline function interp_forcing_T(i, j, k, grid, clock, model_fields, params)
    n_reference, tᶠ, F, λOpen, Nx = params
    n = n_reference[]

    @inbounds begin
        F₁ = F[k, n-1]
        F₂ = F[k, n]
        t₁ = tᶠ[n-1]
        t₂ = tᶠ[n]
    end

    # Linear interpolation
    t = clock.time
    dFdt = (F₂ - F₁) / (t₂ - t₁)

    F = F₁ + dFdt * (t - t₁)

    # First condition (from radiation)
    condition1 = (i == Nx)
    
    # Calculate the forcing terms
    radiation_term = @inbounds ifelse(condition1, -λOpen * (model_fields.T[i, j, k] - F), 0) 

    # Apply logic: Combine both terms if both conditions are true
    return @inbounds radiation_term
end

# adapted from https://github.com/glwagner/TropicalTurbulantics.jl/blob/main/tropical_turbulence_setup.jl
# tropical_turbulence_setup
function forcing_fields_from_file(arch = GPU();
    user_z = -[19, 15, 10, 5, 0],  # Load the "original" model z-grid
    Nz_model = 5,
    datapath = joinpath(homedir(), "BadgerArtifacts", "Varna_BRY.jld2"))
        
    # Load data
    file = jldopen(datapath)


    zforc = file["data_dict"]["depth"]

    # Add one more cell interface at z = z_bottom
    # Not sure if this is correct. However, we must have
    # that length(z) = Nz + 1, where Nz is the number of cells.
    # Otherwise the grid is underdetermined.
    @info string("First cell interface is z[1] = ", zforc[1])
    @info string("Adding one more interface at z = z_bottom - 5m")

    zforc = vcat(zforc[1] - 5, zforc)

    @info string("First cell interface is now z[1] = ", zforc[1])
    
    # what is assert?
    # enforce a condition at runtime, primarily for debugging and validation purposes.
    # If the condition in the @assert statement evaluates to false, Julia raises an AssertionError and stops the program.
    # It is typically used to ensure that assumptions about the program's state hold true during execution.
    # @assert length(z) == Nz + 1  

    # Forcing time
    forcing_times = file["data_dict"]["time"]

    # Large-scale state
    # Can use this to implement restoring --- not implemented yet.
    U = file["data_dict"]["uo"]
    V = file["data_dict"]["vo"]
    T = file["data_dict"]["thetao"]
    S = file["data_dict"]["so"]

    println(size(T)) 
    close(file)

    # regrid if necessary
    if !isnothing(user_z)
        Nz_forc = length(zforc) - 1
        Nz_model = length(user_z) - 1

        original_vertical_grid = RectilinearGrid(size=Nz_forc; z=zforc, topology=(Flat, Flat, Bounded))
        new_vertical_grid      = RectilinearGrid(size=Nz_model; z=user_z, topology=(Flat, Flat, Bounded))

        orig_field = CenterField(original_vertical_grid)
        new_field = CenterField(new_vertical_grid)

        Nt = length(forcing_times)
        new_U = zeros(Nz_model, Nt)
        new_V = zeros(Nz_model, Nt)
        new_T = zeros(Nz_model, Nt)
        new_S = zeros(Nz_model, Nt)

        for (new, orig) in [
                            (new_U,  U),
                            (new_V,  V),
                            (new_T,  T),
                            (new_S,  S)]

        for n = 1:Nt
            orig_field .= reshape(orig[:, n], 1, 1, Nz_forc)
            regrid!(new_field, orig_field)
            new[:, n] .= interior(new_field, :)
        end
    end

    # Just pretend nothing happened
    U  = on_architecture(arch, new_U) 
    V  = on_architecture(arch, new_V) 
    T  = on_architecture(arch, new_T)   
    S  = on_architecture(arch, new_S) 
    z  = znodes(new_vertical_grid, Face())
    end

    # Convert CPU arrays to GPU arrays if necessary:
    tᶠ = on_architecture(arch, forcing_times)

    # u_forcing = Forcing(interp_forcing, discrete_form=true, parameters=(n, tᶠ))
    # v_forcing = Forcing(interp_forcing, discrete_form=true, parameters=(n, tᶠ))
    # T_forcing = Forcing(interp_forcing, discrete_form=true, parameters=(n, tᶠ))

    return U, V, T, S, z, tᶠ
end

#### BOTTOM DRAG ####
function forcing_bottom_drag(bottom_drag_coefficient)
        Fu = Forcing(
            u_immersed_bottom_drag,
            discrete_form = true,
            parameters = bottom_drag_coefficient,
        )
        Fv = Forcing(
            v_immersed_bottom_drag,
            discrete_form = true,
            parameters = bottom_drag_coefficient,
        )
        return (u = Fu, v = Fv, )
end

#### radiation QuasiOpenBoundary ONLY ####
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

#### RIVERS ONLY ####
function forcing_rivers_S(Nz)
    λ = 1 / (30minutes)  # Relaxation timescale [s⁻¹].

    S_source = 0.1

    source_index = (1, 13, Nz)

    S_point_source(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.S[i, j, k] - S_source), 0)
    Sforcing = Forcing(S_point_source, field_dependencies = :S, discrete_form = true)

    return (S = Sforcing, )
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

    return (S = Sforcing, NO₃ = NO₃forcing, )
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
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
    return @inbounds ifelse(condition1 || condition2,
        (condition1 ? radiation_term : 0) + (condition2 ? point_source_term : 0),
        0
    )
end

#### FINAL FORCING ####

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

    # CUDA.@allowscalar user_z = Array(znodes(grid[], Center()))
    # U, V, T, S, z, tᶠ = forcing_fields_from_file(user_z=user_z, Nz_model=Nz)

    TForcing = Forcing(combined_forcing_func_T,
            field_dependencies = :T, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=Tsrc, src_loc=src_loc, external_value=external_values.T),
            discrete_form = true)
    SForcing = Forcing(combined_forcing_func_S,
            field_dependencies = :S, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=Ssrc, src_loc=(1, 13, Nz), external_value=external_values.S),
            discrete_form = true)
    NUTForcing = Forcing(combined_forcing_func_NUT,
            field_dependencies = :NUT, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=NUTsrc, src_loc=(1, 13, Nz), external_value=external_values.NUT),
            discrete_form = true)
    DOMForcing = Forcing(combined_forcing_func_DOM,
            field_dependencies = :DOM, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=DOMsrc, src_loc=src_loc, external_value=external_values.DOM),
            discrete_form = true)
    O2Forcing = Forcing(combined_forcing_func_O2,
            field_dependencies = :O₂, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=O2src, src_loc=src_loc, external_value=external_values.O₂),
            discrete_form = true)
    PForcing = Forcing(combined_forcing_func_P,
            field_dependencies = :P, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=Psrc, src_loc=src_loc, external_value=external_values.P),
            discrete_form = true)
    HETForcing = Forcing(combined_forcing_func_HET,
            field_dependencies = :HET, parameters=(Nx = Nx, λRiver = 0, λOpen = λOpen, Vsrc=HETsrc, src_loc=src_loc, external_value=external_values.HET),
            discrete_form = true)
    
    # -λRiver here, because the poin source is on the upper boundary and should point down
    ContForcing = Forcing(combined_forcing_func_C,
            field_dependencies = :C, parameters=(Nx = Nx, λRiver = λRiver, λOpen = λOpen, Vsrc=Csrc, src_loc=src_loc, external_value=external_values.C),
            discrete_form = true)
            
    # n_reference, tᶠ, F, λOpen, Nx = params
    # TForcing = Forcing(interp_forcing_T, field_dependencies = :T,
    #                     parameters=(n_reference = n, tᶠ = tᶠ, F=T, λOpen = λOpen, Nx = Nx),
    #                     discrete_form = true)
    forcing_bottom = forcing_bottom_drag(bottom_drag_coefficient)

    # forcing_rivers = forcing_rivers_S(Nz)
    final_forcing = merge(forcing_bottom,
    (T = TForcing,
    S = SForcing,
    C = ContForcing,
    NUT = NUTForcing,
    DOM = DOMForcing,
    O₂ = O2Forcing,
    P = PForcing,
    HET = HETForcing,
    ))

    return final_forcing 
end


# # # TEST
# include("utils.jl")
# netcdf_to_jld2(joinpath(homedir(), "BadgerArtifacts", "Varna_BRY_1lat.nc"),
#                 joinpath(homedir(), "BadgerArtifacts", "Varna_BRY.jld2")
# )
# U, V, T, S, z, tᶠ = forcing_fields_from_file()