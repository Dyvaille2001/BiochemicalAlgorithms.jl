"""
    compute_b_factors(vals, vecs, inv_sqrt_M, N)

Calculate atomic B-factors (Debye-Waller factors) from Normal Modes.

B-factors represent the mean-square displacement of atoms. This function
summons the contributions from all non-trivial vibrational modes (k ≥ 7).

Compute ANM-derived B-factors from vibrational modes.
(Structural Flexibility)B-factors (a measure of how much an atom "jiggles" in an X-ray structure)
The first six modes (rigid-body translation/rotation) are excluded.
Skip negative eigenvalues safely (if any) by treating them as zero contribution.

Parameters
----------
vals : Vector{Float64}
    The eigenvalues (λₖ = ωₖ²).
vecs : Matrix{Float64}
    The mass-weighted eigenvectors (mode shapes).
inv_sqrt_M : Diagonal{Float64}
    The inverse square root mass matrix (M⁻¹/²) used to transform
    modes back to Cartesian space.
N : Int
    Total number of atoms.

Returns
-------
fluctuations : Vector{Float64}
    The predicted B-factors for each atom (8π²/3 * <ΔR²>).

Notes
-----
- The loop starts at index 7 to skip the 6 rigid-body modes (λ₁₋₆ ≈ 0).
- Displacement for each mode is weighted by the inverse eigenvalue (1/λₖ).
"""
function compute_b_factors(vals, vecs, inv_sqrt_M, N)
    fluctuations = zeros(N)

    # Sum contributions from all vibrational modes (k = 7 to 3N)
    for k in 7:length(vals)
        # Weighting factor: 1/λₖ (equivalent to 1/ωₖ²)
        inv_freq_sq = 1.0 / vals[k]
        mode_k = vecs[:, k]

        for i in 1:N
            # Indices for the x, y, z components of atom i
            idx = (3i-2):(3i)
            
            # Transform mode from mass-weighted to Cartesian displacement: d = M⁻¹/² * v
            disp = inv_sqrt_M[idx, idx] * mode_k[idx]
            
            # Accumulate the mean-square fluctuation
            fluctuations[i] += sum(disp.^2) * inv_freq_sq
        end
    end

    # Apply the Debye-Waller scaling factor: (8π²/3)
    return fluctuations .* (8.0 * pi^2 / 3.0)
end