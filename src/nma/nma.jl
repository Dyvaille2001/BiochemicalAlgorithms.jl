module NMA

using LinearAlgebra
using StaticArrays
using Plots

import ..BiochemicalAlgorithms: atoms, update!, compute_energy!


# Include all your files
include("build_hessian.jl")
include("compute_modes.jl")
include("energy_along_mode.jl")
include("extract_coordinates_masses.jl")
include("plot_energy_displacement.jl")
include("vmd_mode_visuals.jl")
include("compute_frequencies.jl")
include("compute_b_factors.jl")
include("frequencies_thz.jl")


# Export public API
export build_hessian,
       compute_modes,
       energy_along_mode!,
       extract_coordinates_masses,
       plot_energy_displacement,
       write_vmd_mode,
       compute_frequencies,
       compute_b_factors,
       frequencies_thz

end