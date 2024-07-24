using Interpolations: LinearInterpolation

function initial_conditions_temp_salt_3d_predefined(setup_grid)
    z_middle = setup_grid.z_middle
    Nx, Ny, Nz = setup_grid.Nx, setup_grid.Ny, setup_grid.Nz

    z_ini = -reverse([0.5, 1, 2, 3, 4, 6, 8, 9, 10, 12, 16, 20])
    tprof = reverse([20, 20, 20, 20, 18, 15, 14, 13, 12, 12, 12, 12])
    itp = LinearInterpolation(z_ini, tprof)
    tprof_target = itp(z_middle)

    T₀ = Array{Float64}(undef, Nx, Ny, Nz)
    for i = 1:Nx
        for j = 1:Ny
            T₀[i, j, :] = tprof_target
        end
    end

    sprof = reverse([14, 14.1, 14.2, 14.5, 14.8, 15, 15.1, 15.2, 15.3, 15.4, 15.5, 15.5])
    itps = LinearInterpolation(z_ini, sprof)
    sprof_target = itps(z_middle)

    S₀ = Array{Float64}(undef, Nx, Ny, Nz)
    for i = 1:Nx
        for j = 1:Ny
            S₀[i, j, :] = sprof_target
        end
    end

    return T₀, S₀
end
