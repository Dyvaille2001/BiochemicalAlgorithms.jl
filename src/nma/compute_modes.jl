"""
    mass_weight_and_solve(H, masses)

Transform the Hessian to mass-weighted coordinates and solve the eigenvalue problem.

The mass-weighting transformation corrects the Hessian for atomic masses, 
turning the problem into a search for vibrational frequencies (eigenvalues) 
and normalized mode shapes (eigenvectors).

Parameters
----------
H : AbstractMatrix{Float64}
    The 3N × 3N Hessian matrix (either Sparse or Dense).
masses : Vector{Float64}
    A vector of length N containing the mass of each atom.

Returns
-------
vals : Vector{Float64}
    The eigenvalues:λ (square of the angular frequencies, λ=ω², represent the energetic cost or stiffness of a motion.).
vecs : Matrix{Float64}
    The mass-weighted eigenvectors (columns of matrix describing the mode shapes, & are mass-weighted normal modes

Notes
-----
- The first 6 eigenvalues should be approximatelyzero (Rigid Body Modes, degrees of freedom (3 translations, 3 rotations)).
- The eigenvectors returned are orthonormal in mass-weighted space.
- The function uses `Symmetric()` to ensure numerical stability during 
  decomposition.
"""
function compute_modes(H, masses)
    # 1. Expand N masses to 3N degrees of freedom (x, y, z for each atom)
    m_weighted = repeat(masses, inner=3)
    
    # 2. Create the inverse square root diagonal matrix: M^(-1/2)
    # This acts as the transformation matrix
    inv_sqrt_m = 1.0 ./ sqrt.(m_weighted)
    inv_sqrt_M = Diagonal(inv_sqrt_m)
    
    # 3. Apply the transformation: H_weighted = M^(-1/2) * H * M^(-1/2)
    # This effectively scales rows and columns by 1/sqrt(m_i * m_j)
    H_weighted = inv_sqrt_M * H * inv_sqrt_M
    
    # 4. Eigen-decomposition
    # We wrap in Symmetric() to help the solver recognize the matrix type
    res = eigen(Symmetric(H_weighted))
    
    return res.values, res.vectors, inv_sqrt_M
end
