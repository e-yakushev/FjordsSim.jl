# FjordsSim.jl

## Installation

1. Clone the git repository `git clone https://github.com/AquaBadgerForge/FjordsSim.jl.git`.
2. Move to the downloaded folder `cd FjordsSim.jl`.
3. Run Julia REPL and activate the FjordsSim environment `julia --project=.`.
4. Enter the Pkg REPL by pressing `]` from Julia REPL.
5. Type `instantiate` to 'resolve' a `Manifest.toml` from a `Project.toml` to install and precompile dependency packages.

## Run an example Oslofjord simulation

1. Download the [grid and forcing](https://www.dropbox.com/scl/fo/oidmjtcbn8k3o5h4wgje4/ABouiZ1gYqu8PKVWrHPpKGM?rlkey=d416l5k5k9ptli0q1qsrfm9qs&st=7mduj9dt&dl=0).
2. In `FjordsSim.jl/app/oslofjord/setup.jl` it is possible to specify the location of the input data files.
By default, the files should be in `$HOME/FjordsSim_data/oslofjord/`.
Also, it is possible to specify the result folder destination.
By default, the result will go to `$HOME/FjorsSim_results/oslofjord/`.
3. Run `julia --project app/oslofjord/simulation_105to232.jl` or `simulation_45to82.jl`.
This will generate netcdf (or for some outputs `*.jld2`) results files.

![example_result](./artifacts/phytoplankton_multi.png)
