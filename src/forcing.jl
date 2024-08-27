using Oceananigans: Forcing

function rivers_forcing(Nz)
    ## River forcing
    λ = 1 / (30minutes)  # Relaxation timescale [s⁻¹].

    # Temperature and salinity of the meltwater outflow.
    T_source = 20
    S_source = 0.5   # (G. Shtereva & Dzhurova, 2006a).

    # Index of the point source at the middle of the western wall.
    source_index = (1, 13, Nz)

    # Point source

    # @inline T_point_source(i, j, k, grid, time, U, C, p) =
    #     @inbounds ifelse((i, j, k) == p.source_index, -p.λ * (C.T[i, j, k] - p.T_source), 0)

    # @inline S_point_source(i, j, k, grid, time, U, C, p) =
    #     @inbounds ifelse((i, j, k) == p.source_index, -p.λ * (C.S[i, j, k] - p.S_source), 0)

    T_point_source(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.T[i, j, k] - T_source), 0)

    S_point_source(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.S[i, j, k] - S_source), 0)

    # params = (source_index=source_index, T_source=T_source, S_source=S_source, λ=λ)
    params = (T_source = T_source, S_source = S_source, λ = λ)

    # Tforcing = Forcing(T_point_source, parameters=params)
    Tforcing = Forcing(T_point_source, field_dependencies = :T, discrete_form = true)
    Sforcing = Forcing(S_point_source, field_dependencies = :S, discrete_form = true)

    # return (T = Tforcing, S = Sforcing)  
    return (S = Sforcing,)  
end