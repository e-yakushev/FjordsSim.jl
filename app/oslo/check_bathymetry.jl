# Copyright 2024 The FjordsSim Authors.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

using CairoMakie

include("setup.jl")

using .FjordsSim: grid_bathymetry_from_lat_lon

## Model Setup
sim_setup = setup_oslo()
underlying_grid, bottom_height = grid_bathymetry_from_lat_lon(sim_setup.grid_parameters...)

# For plotting
bottom_height.data[bottom_height.data .>= 0] .= NaN

fig = Figure(size = (400, 400))
ax  = Axis(fig[1, 1])
hm = heatmap!(ax, bottom_height, colormap = :deep, colorrange = (-600, 0))
cb = Colorbar(fig[0, 1], hm, label = "Bottom height (m)", vertical = false)
hidedecorations!(ax)

save(joinpath(homedir(), "FjordsSim_results", "bathymetry.png"), fig)
