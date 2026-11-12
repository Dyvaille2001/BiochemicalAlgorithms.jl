############################################################
# nma_core.jl  -  core ENM / NMA functions (final version)
#
# Conventions
#   H        : 3N x 3N, kJ/(mol·Å²)
#   masses   : amu
#   kfun(r)  : spring constant k_ij in kJ/(mol·Å²)
#              distance-weighted (this thesis): r -> γ / r²   (γ in kJ/mol)
#              standard ANM                   : r -> γ        (γ in kJ/(mol·Å²))
#   modes    : Cartesian displacement of mode k is invsqrtM * vecs[:,k]
############################################################
using LinearAlgebra, Statistics, Printf

# ---------- units -------------------------------------------------
# sqrt((kJ/mol)/(amu·Å²)) = 9.99999992e12 rad/s; 1e13 is used as the
# reference value (relative error 8e-9, far below model uncertainty).
const INTERNAL_TO_RADPS = 1.0e13
const C_CM_PER_S        = 2.99792458e10

const UNIT_FACTORS = Dict(
    :radps => INTERNAL_TO_RADPS,
    :hz    => INTERNAL_TO_RADPS / (2π),
    :thz   => INTERNAL_TO_RADPS / (2π) / 1.0e12,
    :cm1   => INTERNAL_TO_RADPS / (2π * C_CM_PER_S),
    :arb   => 1.0,
)

const ATOMIC_MASSES = Dict("H" => 1.008, "C" => 12.011, "N" => 14.007,
                           "O" => 15.999, "S" => 32.06)

const RESIDUE_MASSES = Dict(
    "ALA"=>71.08, "ARG"=>156.19, "ASN"=>114.10, "ASP"=>115.09, "CYS"=>103.14,
    "GLN"=>128.13, "GLU"=>129.12, "GLY"=>57.05, "HIS"=>137.14, "ILE"=>113.16,
    "LEU"=>113.16, "LYS"=>128.17, "MET"=>131.19, "PHE"=>147.18, "PRO"=>97.12,
    "SER"=>87.08, "THR"=>101.10, "TRP"=>186.21, "TYR"=>163.18, "VAL"=>99.13)

# ---------- spring laws -------------------------------------------
"Distance-weighted ENM used in this thesis: k_ij = γ / r²  (γ in kJ/mol)."
kfun_distance_weighted(gamma=1.0) = r -> gamma / r^2

"Standard (constant-spring) ANM: k_ij = γ  (γ in kJ/(mol·Å²))."
kfun_constant(gamma=1.0) = r -> gamma

# ---------- Hessian -----------------------------------------------
"""
    build_hessian(coords, kfun; cutoff=12.0, pairweight=nothing)

ENM Hessian (3N x 3N). Pair block = kfun(r) · r rᵀ / r² (× pairweight(i,j)).
Pairs with r ≥ cutoff are excluded. `cutoff = Inf` builds a complete graph
(O(N²) pairs) and is intended only for small test systems.
`pairweight(i, j)` (optional) multiplies the stiffness of pair (i, j); it is
used for stiffness-perturbation scans.
"""
function build_hessian(coords::AbstractMatrix{<:Real}, kfun::Function;
                       cutoff=12.0, pairweight=nothing)
    N  = size(coords, 1)
    size(coords, 2) == 3 || throw(ArgumentError("coords must be N×3"))
    H  = zeros(3N, 3N)
    c2 = cutoff^2
    d  = zeros(3)
    @inbounds for i in 1:N-1, j in i+1:N
        d[1] = coords[i,1] - coords[j,1]
        d[2] = coords[i,2] - coords[j,2]
        d[3] = coords[i,3] - coords[j,3]
        d2 = d[1]^2 + d[2]^2 + d[3]^2
        (d2 > 0 && d2 < c2) || continue
        s = kfun(sqrt(d2)) / d2
        pairweight === nothing || (s *= pairweight(i, j))
        for a in 1:3, b in 1:3
            v = s * d[a] * d[b]
            H[3i-3+a, 3j-3+b] -= v
            H[3j-3+b, 3i-3+a] -= v
            H[3i-3+a, 3i-3+b] += v
            H[3j-3+a, 3j-3+b] += v
        end
    end
    return H
end

"Backwards-compatible shorthand for the distance-weighted model."
build_hessian(coords; cutoff=12.0, gamma=1.0) =
    build_hessian(coords, kfun_distance_weighted(gamma); cutoff=cutoff)

# ---------- zero modes --------------------------------------------
"""
    count_zero_modes(vals; tol=1e-8)

Number of eigenvalues with |λ| < tol · max|λ|. Returns length(vals) if every
eigenvalue is exactly zero (degenerate/empty system).
"""
function count_zero_modes(vals; tol=1e-8)
    isempty(vals) && return 0
    m = maximum(abs, vals)
    m == 0 && return length(vals)
    return count(x -> abs(x) < tol * m, vals)
end

"""
    first_vibrational_mode(vals; tol=1e-8, expected=6)

Index of the first non-zero eigenvalue (1-based). Warns if the count differs
from `expected` — treat that warning as a hard failure in production runs.
"""
function first_vibrational_mode(vals; tol=1e-8, expected=6)
    nz = count_zero_modes(vals; tol=tol)
    nz == expected || @warn "Expected $expected zero modes, found $nz " *
                            "(disconnected network? cutoff too small? strained input?)"
    return nz + 1
end

# ---------- modes --------------------------------------------------
"""
    compute_modes(H, masses) -> vals, vecs, invsqrtM

Solves the mass-weighted problem Hmw v = λ v with Hmw = M^(-1/2) H M^(-1/2).
`vecs` columns are orthonormal mass-weighted eigenvectors.
Cartesian displacement of mode k: invsqrtM * vecs[:,k].
Pass `ones(N)` for an unweighted ENM.
"""
function compute_modes(H::AbstractMatrix, masses::AbstractVector)
    N = length(masses)
    size(H, 1) == 3N || throw(DimensionMismatch("H must be $(3N)×$(3N)"))
    all(masses .> 0) || throw(ArgumentError("masses must be positive"))
    invsqrtM = Diagonal(1.0 ./ sqrt.(repeat(Float64.(masses), inner=3)))
    res = eigen(Symmetric(invsqrtM * H * invsqrtM))
    return res.values, res.vectors, invsqrtM
end

"""
    frequencies(vals; unit=:cm1, tol=1e-8)

Real vector. Zero modes -> 0. Negative eigenvalue -> negative frequency
(= imaginary). Physical only if H is in kJ/(mol·Å²) and masses in amu.
"""
function frequencies(vals; unit=:cm1, tol=1e-8)
    isempty(vals) && return Float64[]
    haskey(UNIT_FACTORS, unit) || throw(ArgumentError("unsupported unit: $unit"))
    f    = UNIT_FACTORS[unit]
    λmax = maximum(abs, vals)
    λmax == 0 && return zeros(Float64, length(vals))
    return map(vals) do λ
        abs(λ) <= tol * λmax ? 0.0 : (λ > 0 ? sqrt(λ) * f : -sqrt(-λ) * f)
    end
end

# ---------- fluctuations and correlations -------------------------
"""
    covariance(H; nzero=nothing, nmodes=nothing, kT=1.0, tol=1e-8)

Pseudo-inverse of the UNWEIGHTED Hessian: ⟨dr drᵀ⟩ = kT · H⁺ (mass-independent).
Use `covariance_mw` when the result is combined with mass-weighted modes.
"""
function covariance(H; nzero=nothing, nmodes=nothing, kT=1.0, tol=1e-8)
    F  = eigen(Symmetric(H))
    nz = isnothing(nzero) ? first_vibrational_mode(F.values; tol=tol) - 1 : nzero
    li = isnothing(nmodes) ? length(F.values) : min(nz + nmodes, length(F.values))
    λmax = maximum(abs, F.values)
    idx = [j for j in nz+1:li if F.values[j] > tol * λmax]
    isempty(idx) && return zeros(size(H))
    λ = F.values[idx]
    V = F.vectors[:, idx]
    return kT .* ((V ./ λ') * V')
end

"""
    covariance_mw(H, masses; kT=1.0, nzero=6, tol=1e-8)

Cartesian covariance from the MASS-WEIGHTED Hessian.

    Hmw = M^(-1/2) H M^(-1/2),  Hmw v = λ v
    u   = M^(-1/2) v             (Cartesian displacement)
    C   = kT Σ_k u_k u_kᵀ / λ_k

Rigid-body modes are removed in mass-weighted space (a centre-of-mass type
gauge). Classical equilibrium fluctuations do not depend on masses; the
mass-weighted route only changes which rigid-body motions are projected out.
With equal masses it equals `covariance`.
"""
function covariance_mw(H::AbstractMatrix, masses::AbstractVector;
                       kT=1.0, nzero::Integer=6, tol=1e-8)
    N = length(masses)
    size(H) == (3N, 3N) || throw(DimensionMismatch("H must be $(3N)×$(3N)"))
    s   = 1.0 ./ sqrt.(repeat(Float64.(masses), inner=3))
    Hmw = Diagonal(s) * H * Diagonal(s)
    F   = eigen(Symmetric(Hmw))
    nzero < length(F.values) || error("nzero=$nzero ≥ 3N=$(length(F.values))")
    λmax = maximum(abs, F.values)
    idx  = [j for j in nzero+1:length(F.values) if F.values[j] > tol * λmax]
    isempty(idx) && return zeros(size(H))
    λ = F.values[idx]
    V = F.vectors[:, idx]
    U = s .* V
    return kT .* ((U ./ λ') * U')
end

"Per-atom mean-square fluctuation (Å²) from a 3N x 3N covariance."
function msf_from_cov(C::AbstractMatrix)
    N = size(C, 1) ÷ 3
    return [C[3i-2,3i-2] + C[3i-1,3i-1] + C[3i,3i] for i in 1:N]
end

"Dynamical cross-correlation matrix from a 3N x 3N covariance."
function dccm_from_cov(C::AbstractMatrix)
    N   = size(C, 1) ÷ 3
    msf = msf_from_cov(C)
    D   = zeros(N, N)
    for i in 1:N, j in 1:N
        cij = C[3i-2,3j-2] + C[3i-1,3j-1] + C[3i,3j]
        den = sqrt(msf[i] * msf[j])
        D[i, j] = den > 0 ? cij / den : 0.0
    end
    return D
end

# ---------- gamma calibration (B-factors) -------------------------
"""
    calibrate_gamma(H1, bexp; kT=2.494, mask=nothing, weighted=false, masses=nothing, nzero=6)

H1   : Hessian built with γ = 1 (spring law γ·f(r))
bexp : experimental Cα B-factors (Å²), same atom order. Non-finite and
       non-positive entries are masked automatically.
Returns (γ, Pearson r, slope s). γ has the units that make γ·f(r) a spring
constant in kJ/(mol·Å²): kJ/mol for the distance-weighted law k = γ/r², and
kJ/(mol·Å²) for the constant-spring law.

B_pred = s · g,  s = Σ g b / Σ g²,  γ = (8π²/3) kT / s
Use ONE γ for WT and mutant.
"""
function calibrate_gamma(H1, bexp; kT=2.494, mask=nothing,
                         weighted::Bool=false, masses=nothing, nzero=6)
    N = size(H1, 1) ÷ 3
    length(bexp) == N || throw(DimensionMismatch("expected $N B-factors, got $(length(bexp))"))

    m = mask === nothing ? trues(N) : copy(mask)
    m .&= isfinite.(bexp) .& (bexp .> 0)
    count(m) >= 3 || throw(ArgumentError("fewer than 3 valid B-factors"))

    C = weighted ? covariance_mw(H1, masses; kT=1.0, nzero=nzero) :
                   covariance(H1; kT=1.0)
    g = msf_from_cov(C)

    s = sum(g[m] .* bexp[m]) / sum(abs2, g[m])
    γ = (8π^2 / 3) * kT / s
    return γ, cor(g[m], bexp[m]), s
end

# ---------- subspace comparison -----------------------------------
"""
    rmsip(V, W, k; orth_warn=true)

Root-mean-square inner product over the first k columns of V and W.
Both inputs must be orthonormal for the 1/k normalisation to be meaningful;
a warning is emitted if they are not.
"""
function rmsip(V::AbstractMatrix, W::AbstractMatrix, k::Integer; orth_warn=true)
    size(V, 1) == size(W, 1) || throw(DimensionMismatch("V and W have different rows"))
    k = min(k, size(V, 2), size(W, 2))
    k >= 1 || throw(ArgumentError("need at least one mode"))
    A = V[:, 1:k]; B = W[:, 1:k]
    if orth_warn
        eV = maximum(abs, A'A - I(k))
        eW = maximum(abs, B'B - I(k))
        (eV > 1e-6 || eW > 1e-6) &&
            @warn "rmsip: input modes are not orthonormal (dev $eV / $eW)"
    end
    return sqrt(sum(abs2, A'B) / k)
end

"""
    overlap_matrix(V, W, k)

|v_i · w_j| for the first k modes of each set (k × k, entries in [0, 1]
provided V and W are orthonormal).
"""
function overlap_matrix(V::AbstractMatrix, W::AbstractMatrix, k::Integer)
    size(V, 1) == size(W, 1) || throw(DimensionMismatch("V and W have different rows"))
    k = min(k, size(V, 2), size(W, 2))
    return abs.(V[:, 1:k]' * W[:, 1:k])
end

# ---------- VMD / PDB output --------------------------------------
"""
    write_mode_pdb(path, coords, disp3N, atomnames, resnames, chains, resnums;
                   n_frames=30, rms_amp=1.5, elements=nothing)

Multi-model PDB animation of one mode. `disp3N` is a Cartesian displacement
vector of length 3N; it is rescaled so the RMS per-atom displacement equals
`rms_amp` (Å) at the frame where sin = 1. Pass 3-letter residue names.
"""
function write_mode_pdb(path, coords, disp3N, atomnames, resnames, chains, resnums;
                        n_frames=30, rms_amp=1.5, elements=nothing)
    N = size(coords, 1)
    length(disp3N) == 3N || throw(DimensionMismatch("disp3N must have length 3N"))
    d  = permutedims(reshape(disp3N, 3, :))              # N x 3
    nd = norm(d)
    nd > 0 || throw(ArgumentError("zero displacement vector"))
    d  = d .* (sqrt(N) / nd)                             # RMS per atom = 1 Å

    els = elements === nothing ?
          [lpad(string(first(strip(n))), 2) for n in atomnames] : elements

    open(path, "w") do io
        for f in 1:n_frames
            s = rms_amp * sin(2π * f / n_frames)
            @printf(io, "MODEL     %4d\n", f)
            for i in 1:N
                nm = length(atomnames[i]) < 4 ? " " * rpad(atomnames[i], 3) : atomnames[i]
                @printf(io, "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n",
                        i, nm, resnames[i], chains[i], resnums[i],
                        coords[i,1] + s*d[i,1], coords[i,2] + s*d[i,2],
                        coords[i,3] + s*d[i,3], 1.0, 0.0, els[i])
            end
            println(io, "ENDMDL")
        end
        println(io, "END")
    end
    return path
end

# ---------- rank correlation --------------------------------------
"Tie-aware (average) ranks."
function ranks(x::AbstractVector{<:Real})
    n   = length(x)
    idx = sortperm(x)
    r   = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && x[idx[j+1]] == x[idx[i]]
            j += 1
        end
        avg = (i + j) / 2
        for k in i:j
            r[idx[k]] = avg
        end
        i = j + 1
    end
    return r
end

"Spearman rank correlation with tie handling; NaN if undefined."
function spearman(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(DimensionMismatch("x and y differ in length"))
    length(x) >= 3 || return NaN
    rx, ry = ranks(x), ranks(y)
    (std(rx) > 0 && std(ry) > 0) || return NaN
    return cor(rx, ry)
end
