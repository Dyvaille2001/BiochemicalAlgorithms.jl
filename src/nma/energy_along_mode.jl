"""
    energy_along_mode!(ff, coords, vecs, mode_idx, inv_M; steps=40, max_disp=2.0)

Calculate the potential energy profile of a system as it moves along a specific normal mode.

This function "scans" the energy landscape by displacing the atomic coordinates
along the k-th eigenvector and evaluating the force field energy at each step.

Parameters
----------
ff : ForceField
    The molecular mechanics force field object (e.g., from Molly.jl).
coords : Matrix{Float64}
    The equilibrium (minimum energy) Cartesian coordinates.
vecs : Matrix{Float64}
    The mass-weighted eigenvectors (mode shapes).
mode_idx : Int
    The index of the mode to scan (λₖ).
inv_M : Diagonal{Float64}
    The inverse mass matrix (M⁻¹/² or M⁻¹) for Cartesian transformation.
steps : Int, optional
    Number of points to sample along the mode (default: 40).
max_disp : Float64, optional
    The maximum displacement amplitude in Å (default: 2.0).

Returns
-------
alphas : Vector{Float64}
    The displacement scaling factors (α) used during the scan.
energies : Vector{Float64}
    The potential energy values calculated at each displacement.

Notes
-----
- The resulting plot of Energy vs. α should be parabolic (U = ½kα²) if 
  the system is in a true harmonic minimum.
- This function modifies the underlying system coordinates (`ff.system`) 
  during the calculation.
"""
function energy_along_mode!(ff, coords, vecs, mode_idx, inv_M; steps=40, max_disp=2.0)
    energies = Float64[]
    
    # Create a range of displacement magnitudes from -max to +max
    # α = 0.0 represents the equilibrium structure
    alphas = range(-max_disp, max_disp, length=steps)

    # Transform mass-weighted eigenvector to Cartesian displacement: dₖ = M⁻¹/² * vₖ
    mode = inv_M * vecs[:, mode_idx]
    mode = mode / norm(mode)

    atoms_list = collect(atoms(ff.system))
    N = size(coords, 1)

    for α in alphas
        # Update the system coordinates: R_new = R_equilibrium + α * dₖ
        for i in 1:N
            dx, dy, dz = mode[3i-2], mode[3i-1], mode[3i]
            
            # Update the atom position in the ForceField system
            # SVector{3} is used for high performance in Molly.jl
            atoms_list[i].r = SVector{3, Float32}(
                coords[i, 1] + α * dx, 
                coords[i, 2] + α * dy, 
                coords[i, 3] + α * dz
            )
        end

        # Synchronize force field state and calculate the current Potential Energy
        update!(ff)
        push!(energies, compute_energy!(ff))
    end

    return collect(alphas), energies
end