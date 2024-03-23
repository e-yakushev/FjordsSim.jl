using Oceananigans

include("../../src/FjordsSim.jl")

using .FjordsSim

fjords_setup = FjordsSetup(
    CPU(),
    (100, 100, 10),
    (58.8, 59.9),
    (10.1, 11.1),
    (-500, 0),
    (4, 4, 4),
    joinpath(homedir(), "fjords_data"),
    "ETOPO_2022_v1_15s_N60E000_surface.nc",
    nothing,
    0,
)
