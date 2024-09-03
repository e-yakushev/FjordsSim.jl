using Oceananigans: Forcing

function rivers_forcing(Nz)
    λ = 1 / (30minutes)  # Relaxation timescale [s⁻¹].

    S_source = 0.1

    source_index = (1, 13, Nz)

    S_point_source(i, j, k, grid, clock, model_fields) =
        @inbounds ifelse((i, j, k) == (1, 13, Nz), -λ * (model_fields.S[i, j, k] - S_source), 0)
    Sforcing = Forcing(S_point_source, field_dependencies = :S, discrete_form = true)

    return (S = Sforcing, )
end

function rivers_forcing_with_NO₃(Nz)
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
