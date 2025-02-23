# FjordsSim.jl
## Installation
1. Clone the git repository `git clone https://github.com/AquaBadgerForge/FjordsSim.jl.git`.
2. Move to the downloaded folder `cd FjordsSim.jl`.
3. Run Julia REPL and activate the FjordsSim environment `julia --project=.`.
4. Enter the Pkg REPL by pressing `]` from Julia REPL.
5. Type `instantiate` to 'resolve' a `Manifest.toml` from a `Project.toml` to install and precompile dependency packages.

## Folder Setup:
1. Create two directories in your home directory: Fjordssim_data and Fjordssim_results.
2. Inside each directory, create a subfolder for each case (e.g., Fjordssim_data/varna).
3. The file paths for each case are configured in setup.jl.
4. Place the bathymetry file and the boundary conditions file into the corresponding case folder (for example, Fjordssim_data/varna).

## Generating the Forcing File:
1. Run the Jupyter Notebook located at fjordssim-notebooks/varna/Varna_BRY.ipynb to generate the forcing file.
2. You can force any tracers as well as the velocity components (u and v). Make sure that the tracer names match those used in Fjordssim.
3. In the notebook, specify the forcing values and lambda parameters.
