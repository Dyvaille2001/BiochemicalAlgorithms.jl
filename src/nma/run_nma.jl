"""
    NormalModeAnalysis

Container for a complete vibrational normal mode calculation.

The fields follow the same calculation flow as MMTK's `VibrationalModes`:
Cartesian coordinates and masses are used to build a Hessian, the Hessian is
mass-weighted and diagonalized, and the eigenvectors are transformed back to
Cartesian displacements.
"""
struct NormalModeAnalysis{T<:Real}
    coords::Matrix{T}
    masses::Vector{T}
    hessian::Matrix{T}
    eigenvalues::Vector{T}
    mass_weighted_modes::Matrix{T}
    inv_sqrt_mass_matrix::Diagonal{T, Vector{T}}
    cartesian_modes::Matrix{T}
    frequencies::Vector{ComplexF64}
end

"""
    normal_mode_analysis(sys; cutoff=12.0, gamma=1.0, unit=:thz)

Run an anisotropic-network normal mode analysis for a BiochemicalAlgorithms
system-like object.

This is the Julia analogue of the MMTK example:

1. extract Cartesian coordinates and atomic masses,
2. construct the Cartesian Hessian,
3. solve the mass-weighted eigenvalue problem,
4. convert eigenvalues to frequencies,
5. transform modes back to Cartesian coordinates.

The first six modes are the rigid-body translation/rotation modes for a free
molecule. Therefore the first chemically meaningful vibrational mode is mode 7
in Julia's 1-based indexing, corresponding to `modes[6]` in the MMTK Python
example.
"""
function normal_mode_analysis(sys; cutoff=12.0, gamma=1.0, unit=:thz)
    coords, masses = extract_coordinates_masses(sys)
    H = build_hessian(coords; cutoff=cutoff, gamma=gamma)
    vals, vecs, inv_sqrt_M = compute_modes(H, masses)
    freqs = frequencies_thz(vals, masses; unit=unit)
    cartesian_modes = Matrix(inv_sqrt_M * vecs)

    T = promote_type(eltype(coords), eltype(masses), eltype(vals))
    return NormalModeAnalysis{T}(
        Matrix{T}(coords),
        Vector{T}(masses),
        Matrix{T}(H),
        Vector{T}(vals),
        Matrix{T}(vecs),
        Diagonal(Vector{T}(diag(inv_sqrt_M))),
        Matrix{T}(cartesian_modes),
        freqs,
    )
end

"""
    mode_displacements(nma, mode_index; normalize=true, amplitude=1.0)

Return the Cartesian displacement vectors for one normal mode as an `N x 3`
matrix. Use `mode_index=7` for the first non-rigid vibrational mode of a free
molecule.
"""
function mode_displacements(nma::NormalModeAnalysis, mode_index::Integer;
                            normalize::Bool=true, amplitude::Real=1.0)
    1 <= mode_index <= length(nma.eigenvalues) ||
        throw(BoundsError(nma.eigenvalues, mode_index))

    mode = copy(nma.cartesian_modes[:, mode_index])
    if normalize
        mode_norm = norm(mode)
        mode_norm > 0 && (mode ./= mode_norm)
    end
    mode .*= amplitude
    return permutedims(reshape(mode, 3, :))
end

"""
    first_vibrational_mode(nma) -> Int

Return the conventional first non-rigid mode index. This is 7 for a nonlinear
free molecule in Julia's 1-based indexing.
"""
first_vibrational_mode(::NormalModeAnalysis) = 7
