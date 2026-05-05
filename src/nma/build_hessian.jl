"""
    build_hessian(coords; cutoff=12.0, gamma=1.0)

Construct the Hessian matrix using a harmonic potential and anisotropic interactions.

The Hessian describes harmonic interactions between atoms within a
distance cutoff. Each interacting pair contributes a 3×3 block
corresponding to second derivatives of the harmonic potential.

Parameters
----------
coords : Matrix{Float64}
    N×3 matrix of Cartesian coordinates (Å).
cutoff : Float64, optional
    Interaction cutoff distance in Å (default: 12.0).
gamma : Float64, optional
    Uniform spring constant scaling factor (default: 1.0).

Returns
-------
Matrix{Float64}
    3N × 3N Hessian matrix in Cartesian coordinates.

Notes
-----
- Only atom pairs within `cutoff` contribute.
- The Hessian is symmetric by construction.
- The interaction has the form: H_ij ∝ -γ * (r_ij r_ijᵀ) / r_ij⁴. 
    corresponding to a directional (anisotropic) spring.
- Units depend on the chosen `gamma`.
"""

function build_hessian(coords; cutoff=12.0, gamma=1.0)
    n_atoms = size(coords, 1)
    H = zeros(3n_atoms, 3n_atoms) # Hessian is 3N x 3N (x,y,z for each atom)
    
    for i in 1:n_atoms-1
        for j in i+1:n_atoms
            rij = coords[i, :] - coords[j, :] # Vector between atom i and j
            dist_sq = sum(rij.^2)             # Squared distance
            
            if dist_sq < cutoff^2             # Only interact if within cutoff
                dist = sqrt(dist_sq)
                kij = gamma / dist_sq         # Stiffness scales inversely with distance
                
                # The 3x3 sub-matrix of second derivatives (Anisotropic interaction)
                # Formula: (kij / dist^2) * (rij * rij')
                outer = (rij * rij') * (kij / dist_sq)
                
                idx_i = (3i-2):(3i)           # Indices for atom i (x,y,z)
                idx_j = (3j-2):(3j)           # Indices for atom j (x,y,z)
                
                H[idx_i, idx_j] -= outer      # Interaction between i and j
                H[idx_j, idx_i] -= outer      # Symmetric counterpart
                H[idx_i, idx_i] += outer      # Diagonal (sum of all forces on i)
                H[idx_j, idx_j] += outer      # Diagonal (sum of all forces on j)
            end
        end
    end
    return H
end
