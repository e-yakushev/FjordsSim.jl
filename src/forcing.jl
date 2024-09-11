using Oceananigans: Forcing

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

function forcing_varna(bottom_drag_coefficient, Nz)
    forcing_bottom = forcing_bottom_drag(bottom_drag_coefficient)
    forcing_rivers = forcing_rivers_S(Nz)
    return merge(forcing_bottom, forcing_rivers) 
end
