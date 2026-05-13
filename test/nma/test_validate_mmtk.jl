using Test
using Statistics
using LinearAlgebra

using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA: extract_coordinates_masses,
                                 build_hessian,
                                 compute_modes,
                                 frequencies_thz

# ============================================================
# FILE READERS
# ============================================================

"""
Read Julia frequency file of the form:

# Mode  Eigenvalue  Frequency_THz
1   -5.95e-17   1.22e-8im
2   ...
"""
function read_julia_freq_file(filepath::AbstractString)
    modes = Int[]
    eigvals = Float64[]
    freqs = ComplexF64[]

    for line in eachline(filepath)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue

        parts = split(s)
        @assert length(parts) >= 3 "Invalid Julia frequency file line: $line"

        mode = parse(Int, parts[1])
        eigv = parse(Float64, parts[2])
        fstr = parts[3]

        freq = if endswith(fstr, "im")
            val = parse(Float64, fstr[1:end-2])
            complex(0.0, val)
        else
            complex(parse(Float64, fstr), 0.0)
        end

        push!(modes, mode)
        push!(eigvals, eigv)
        push!(freqs, freq)
    end

    return (; modes, eigvals, freqs)
end

"""
Read MMTK frequency file.

Accepted formats:
1) plain values:
   8.09e-06j
   1.12e-05

2) full console lines:
   Mode 0 with frequency 8.09e-06j
"""
function read_mmtk_frequencies(filepath::AbstractString)
    freqs = ComplexF64[]

    for line in eachline(filepath)
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue

        if occursin("frequency", s)
            s = strip(split(s, "frequency")[end])
        end

        if endswith(s, "j")
            val = parse(Float64, s[1:end-1])
            push!(freqs, complex(0.0, val))
        else
            val = parse(Float64, s)
            push!(freqs, complex(val, 0.0))
        end
    end

    return freqs
end

# ============================================================
# COMPUTE CURRENT JULIA ALAALA NMA RESULTS
# ============================================================

function compute_alaala_modes()
    pdb_path = ball_data_path("../test/data/AlaAla.pdb")
    mol = load_pdb(pdb_path, Float64)

    fdb = FragmentDB{Float64}()
    normalize_names!(mol, fdb)
    reconstruct_fragments!(mol, fdb)
    build_bonds!(mol, fdb)

    ff = AmberFF(mol)
    update!(ff)
    optimize_structure!(ff)

    coords, masses = extract_coordinates_masses(mol)
    H = build_hessian(coords)
    vals, vecs, inv_sqrt_M = compute_modes(H, masses)
    freqs = frequencies_thz(vals, masses; unit=:thz)

    return (; coords, masses, vals, vecs, inv_sqrt_M, freqs)
end

# ============================================================
# HELPERS
# ============================================================

is_imag_mode(z::ComplexF64; tol=1e-12) = abs(imag(z)) > tol && abs(real(z)) ≤ tol
is_real_mode(z::ComplexF64; tol=1e-12) = abs(real(z)) > tol && abs(imag(z)) ≤ tol
is_near_zero_mode(z::ComplexF64; tol=1e-3) = abs(z) < tol

"""
Return positive real frequencies from mode `start_mode` onward.
Used to ignore rigid-body / trivial modes at the beginning.
"""
function nontrivial_real_modes(freqs::AbstractVector{ComplexF64}; start_mode::Int=7)
    out = Float64[]
    for f in freqs[start_mode:end]
        if is_real_mode(f) && real(f) > 0
            push!(out, real(f))
        end
    end
    return out
end

"""
Robust estimate of constant scale factor between two spectra.
If ref ≈ c * test, then c is estimated by median(ref[i]/test[i]).
"""
function median_scale_factor(ref::AbstractVector{<:Real}, test::AbstractVector{<:Real})
    @assert length(ref) == length(test)
    ratios = [ref[i] / test[i] for i in eachindex(ref) if abs(test[i]) > 1e-14]
    return median(ratios)
end

# ============================================================
# TESTS
# ============================================================

@testset "MMTK validation for AlaAla NMA" begin
    pkgroot = pkgdir(BiochemicalAlgorithms)
    outdir = joinpath(pkgroot, "tutorials", "nma", "outputs", "alaala_nma")

    julia_file = joinpath(outdir, "all_eigenvalues_frequencies_thz.txt")
    mmtk_file  = joinpath(outdir, "ala_mmtk_freqthz.txt")

    @test isfile(julia_file)
    @test isfile(mmtk_file)

    julia_saved = read_julia_freq_file(julia_file)
    mmtk_freqs  = read_mmtk_frequencies(mmtk_file)
    computed    = compute_alaala_modes()

    vals  = computed.vals
    freqs = computed.freqs
    N     = size(computed.coords, 1)

    # ========================================================
    # 1) STRUCTURAL CORRECTNESS TESTS
    #
    # These tests check whether the Julia NMA implementation
    # is internally consistent and physically plausible,
    # independent of exact agreement with MMTK.
    #
    # If these fail, the problem is inside the Julia NMA pipeline
    # itself: mode count, saving/loading, rigid-body treatment,
    # or numerical stability.
    # ========================================================

    @testset "Structural correctness: mode count" begin
        # Tests whether the implementation produces the correct
        # total number of normal modes for N atoms: 3N.
        #
        # Implication if this fails:
        # - wrong atom count
        # - wrong Hessian dimension
        # - bug in mode computation pipeline
        @test length(vals) == 3N
        @test length(freqs) == 3N
        @test length(julia_saved.freqs) == 3N
        @test length(mmtk_freqs) == 3N
    end

    @testset "Structural correctness: saved Julia output matches fresh computation" begin
        # Tests whether the saved file is a faithful export of
        # the currently computed frequencies/eigenvalues.
        #
        # Implication if this fails:
        # - export function is wrong
        # - stale output file
        # - mismatch between workflow and current implementation
        @test length(julia_saved.eigvals) == length(vals)
        @test length(julia_saved.freqs) == length(freqs)

        for i in eachindex(vals)
            @test isapprox(julia_saved.eigvals[i], vals[i]; atol=1e-12, rtol=1e-10)
            @test isapprox(real(julia_saved.freqs[i]), real(freqs[i]); atol=1e-12, rtol=1e-10)
            @test isapprox(imag(julia_saved.freqs[i]), imag(freqs[i]); atol=1e-12, rtol=1e-10)
        end
    end

    @testset "Structural correctness: rigid-body / near-zero mode sanity" begin
        # For a minimized non-linear molecule, the first ~6 modes
        # should be near zero (translations + rotations).
        #
        # Implication if this fails:
        # - minimization may be insufficient
        # - Hessian may be wrong
        # - rigid-body structure may be mishandled
        n_zero_julia = count(f -> is_near_zero_mode(f; tol=1e-2), freqs[1:6])
        n_zero_mmtk  = count(f -> is_near_zero_mode(f; tol=1e-2), mmtk_freqs[1:6])

        @test n_zero_julia >= 4
        @test n_zero_mmtk >= 2
    end

    @testset "Structural correctness: imaginary mode count sanity" begin
        # Imaginary frequencies indicate negative curvature directions.
        # Small numbers among the first modes can happen numerically,
        # but many imaginary modes would indicate a non-minimum or a bug.
        #
        # Implication if this fails:
        # - minimization mismatch
        # - Hessian sign convention issue
        # - force-field / coordinate inconsistency
        n_imag_julia = count(f -> is_imag_mode(f), freqs)
        n_imag_mmtk  = count(f -> is_imag_mode(f), mmtk_freqs)

        @info "Imaginary mode counts" julia=n_imag_julia mmtk=n_imag_mmtk

        @test n_imag_julia == 2
        @test n_imag_mmtk == 2
    end

    @testset "Structural correctness: nontrivial frequencies are sorted ascending" begin
        # Tests whether the nontrivial positive frequencies are sorted,
        # as expected after eigenvalue decomposition.
        #
        # Implication if this fails:
        # - eigenpairs not sorted
        # - comparison to MMTK becomes unreliable
        julia_real = nontrivial_real_modes(freqs; start_mode=7)
        @test issorted(julia_real)
    end

    # ========================================================
    # 2) MMTK COMPARISON: BROAD AGREEMENT TEST
    #
    # This asks:
    # "Does the Julia spectrum broadly resemble the MMTK spectrum?"
    #
    # This is NOT yet asking for perfect agreement.
    # It checks whether the general spectral shape is similar.
    #
    # If this passes:
    # - your implementation is broadly consistent with MMTK
    #
    # If this fails:
    # - there is likely a deeper mismatch in Hessian definition,
    #   mass-weighting, minimization endpoint, or model convention.
    # ========================================================

    @testset "MMTK comparison: broad spectral agreement" begin
        # Compare sorted nontrivial positive frequencies.
        # Sorting avoids over-penalizing small mode reordering.
        julia_real = sort(nontrivial_real_modes(freqs; start_mode=7))
        mmtk_real  = sort(nontrivial_real_modes(mmtk_freqs; start_mode=7))

        ncmp = min(length(julia_real), length(mmtk_real))
        @test ncmp > 20

        julia_real = julia_real[1:ncmp]
        mmtk_real  = mmtk_real[1:ncmp]

        corrval = cor(mmtk_real, julia_real)

        @info "Broad spectral agreement" corr=corrval compared_modes=ncmp

        # Broad agreement threshold:
        # high correlation means the overall spectral shape is similar,
        # even if scaling differs.
        @test corrval > 0.95
    end

    # ========================================================
    # 3) MMTK COMPARISON: DIAGNOSTIC TEST
    #
    # This asks:
    # "Are Julia and MMTK frequencies related mostly by a constant factor?"
    #
    # If YES:
    # - likely mainly a unit conversion / prefactor issue
    #
    # If NO:
    # - likely a semantic/model mismatch:
    #   Hessian construction, mass weighting, minimization,
    #   atom ordering, or force-field convention differences.
    # ========================================================

    @testset "MMTK comparison: constant scaling vs semantic mismatch diagnostic" begin
        julia_real = sort(nontrivial_real_modes(freqs; start_mode=7))
        mmtk_real  = sort(nontrivial_real_modes(mmtk_freqs; start_mode=7))

        ncmp = min(length(julia_real), length(mmtk_real))
        julia_real = julia_real[1:ncmp]
        mmtk_real  = mmtk_real[1:ncmp]

        # Fit best constant scale factor
        scale = median_scale_factor(mmtk_real, julia_real)
        scaled_julia = scale .* julia_real

        # Relative errors after best constant scaling
        relerrs = abs.(mmtk_real .- scaled_julia) ./ max.(abs.(mmtk_real), 1e-12)

        # Ratio spread:
        # If ratios are nearly constant, then difference is mostly a scale factor.
        ratios = mmtk_real ./ julia_real
        ratio_cv = std(ratios) / mean(ratios)   # coefficient of variation

        @info "Scaling diagnostic" scale=scale ratio_cv=ratio_cv median_relerr=median(relerrs) max_relerr=maximum(relerrs)

        # We keep these thresholds moderate, because this is a diagnostic test.
        #
        # Interpretation:
        # - ratio_cv small  -> mostly constant scaling
        # - ratio_cv large  -> deeper semantic mismatch
        # - median_relerr small after scaling -> good consistency
        @test isfinite(scale)
        @test median(relerrs) < 0.30
        @test ratio_cv < 0.50
    end

    # ========================================================
    # 4) DIAGNOSTIC PRINT BLOCK
    #
    # This block does not decide validity by itself.
    # It prints evidence to help identify where mismatch arises.
    # ========================================================

    @testset "MMTK diagnostic printout" begin
        julia_real = sort(nontrivial_real_modes(freqs; start_mode=7))
        mmtk_real  = sort(nontrivial_real_modes(mmtk_freqs; start_mode=7))

        ncmp = min(length(julia_real), length(mmtk_real))
        julia_real = julia_real[1:ncmp]
        mmtk_real  = mmtk_real[1:ncmp]

        scale = median_scale_factor(mmtk_real, julia_real)
        scaled_julia = scale .* julia_real
        ratios = mmtk_real ./ julia_real
        relerrs = abs.(mmtk_real .- scaled_julia) ./ max.(abs.(mmtk_real), 1e-12)

        @info "Diagnostic summary" corr=cor(mmtk_real, julia_real) scale=scale ratio_min=minimum(ratios) ratio_median=median(ratios) ratio_max=maximum(ratios)

        for i in 1:min(10, ncmp)
            @info "Mode $(i+6)" julia_thz=julia_real[i] scaled_julia_thz=scaled_julia[i] mmtk_thz=mmtk_real[i] relerr=relerrs[i]
        end

        @test true
    end

    @info "Interpretation hint" begin
    if ratio_cv < 0.1
        "Mostly constant scaling → unit/prefactor issue"
    elseif ratio_cv < 0.3
        "Moderate variation → likely minor model differences"
    else
        "Strong variation → Hessian/model mismatch likely"
    end
end
end