"""
    plot_energy_modes(ff, coords, vecs, inv_M, modes)

Generate a comparison plot of potential energy profiles for multiple normal modes.

This function calls `energy_along_mode!` for each requested mode index,
normalizes the energy (sets the minimum to zero), and overlays them on 
a single plot to visualize the "stiffness" of different vibrations.

Parameters
----------
ff : ForceField
    The molecular mechanics system to evaluate.
coords : Matrix{Float64}
    The equilibrium Cartesian coordinates.
vecs : Matrix{Float64}
    The mass-weighted eigenvectors (mode shapes).
inv_M : Diagonal{Float64}
    The inverse mass matrix (M⁻¹/² or M⁻¹).
modes : Vector{Int}
    A list of mode indices to plot (e.g., [7, 8, 9, 20]).

Returns
-------
p : Plot
    The resulting Plots.jl object.

Notes
-----
- The energy is shifted using E - min(E) so that the equilibrium 
  structure always sits at 0.0 on the Y-axis.
- Steeper parabolas correspond to higher eigenvalues (λₖ).
- Flatter parabolas correspond to low-frequency, flexible modes.
"""
function plot_energy_displacement(ff, coords, vecs, inv_M, modes; filename="energy_vs_displacement.png")
    # Initialize a new plot with a clean theme
    p = plot(title="Energy Landscape along ANM Modes", 
             grid=true, 
             legend=:top)

    for m in modes
        # 1. Scan the energy landscape for mode 'm'
        α, E = energy_along_mode!(ff, coords, vecs, m, inv_M)
        
        # 2. Normalize Energy: Shift such that the minimum (equilibrium) is 0
        E_shifted = E .- minimum(E)
        
        # 3. Add to the existing plot
        # lw=2 makes the lines easier to see in the saved PNG
        plot!(p, α, E_shifted, label="Mode $m", lw=2)
    end

    # Labeling with units if known (typically Å and kcal/mol or kJ/mol)
    xlabel!(p,"Displacement Amplitude (α)")
    ylabel!(p,"Δ Potential Energy")

    # Save the result to disk
    savefig(p, filename)
    
    return p
end