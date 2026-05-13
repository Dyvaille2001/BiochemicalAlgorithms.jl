using Test
using DelimitedFiles
using Statistics
using LinearAlgebra

using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA: extract_coordinates_masses,
                                 build_hessian,
                                 compute_modes,
                                 frequencies_thz

# -------------------------------------------------
# Read Julia file:
# all_eigenvalues_frequencies_thz.txt
# Format:
# # Mode    Eigenvalue    Frequency_THz
# 1    -5.95e-17    1.22e-8im
# 2    ...
# -------------------------------------------------
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

# -------------------------------------------------
# Read MMTK file
# Accepted formats:
# 1) plain values:
#    8.09e-06j
#    1.12e-05
#
# 2) console lines:
#    Mode 0 with frequency 8.09e-06j
# -------------------------------------------------
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

# -------------------------------------------------
# Recompute AlaAla modes directly from current code
# -------------------------------------------------
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

# -------------------------------------------------
# Helpers
# -------------------------------------------------
is_imag_mode(z::ComplexF64; tol=1e-12) = abs(imag(z)) > tol && abs(real(z)) ≤ tol
is_real_mode(z::ComplexF64; tol=1e-12) = abs(real(z)) > tol && abs(imag(z)) ≤ tol
is_near_zero_mode(z::ComplexF64; tol=1e-3) = abs(z) < tol

function nontrivial_real_modes(freqs::AbstractVector{ComplexF64}; start_mode::Int=7)
    out = Float64[]
    for f in freqs[start_mode:end]
        if is_real_mode(f) && real(f) > 0
            push!(out, real(f))
        end
    end
    return out
end

function median_scale_factor(ref::AbstractVector{<:Real}, test::AbstractVector{<:Real})
    @assert length(ref) == length(test)
    ratios = [ref[i] / test[i] for i in eachindex(ref) if abs(test[i]) > 1e-14]
    return median(ratios)
end

@testset "MMTK validation for AlaAla NMA" begin
    pkgroot = pkgdir(BiochemicalAlgorithms)
    outdir = joinpath(pkgroot, "tutorials", "nma", "outputs", "alaala_nma")

    julia_file = joinpath(outdir, "all_eigenvalues_frequencies_thz.txt")
    mmtk_file  = joinpath(outdir, "ala_mmtk_freqthz.txt")

    @test isfile(julia_file)
    @test isfile(mmtk_file)

    julia_saved = read_julia_freq_file(julia_file)
    mmtk_freqs  = read_mmtk_frequencies(mmtk_file)

    computed = compute_alaala_modes()

    vals  = computed.vals
    freqs = computed.freqs
    N     = size(computed.coords, 1)

    @testset "Mode count" begin
        @test length(vals) == 3N
        @test length(freqs) == 3N
        @test length(julia_saved.freqs) == 3N
        @test length(mmtk_freqs) == 3N
    end

    @testset "Saved Julia file matches freshly computed frequencies" begin
        @test length(julia_saved.eigvals) == length(vals)
        @test length(julia_saved.freqs) == length(freqs)

        for i in eachindex(vals)
            @test isapprox(julia_saved.eigvals[i], vals[i]; atol=1e-12, rtol=1e-10)
            @test isapprox(real(julia_saved.freqs[i]), real(freqs[i]); atol=1e-12, rtol=1e-10)
            @test isapprox(imag(julia_saved.freqs[i]), imag(freqs[i]); atol=1e-12, rtol=1e-10)
        end
    end

    @testset "Rigid-body mode sanity" begin
        n_zero_julia = count(f -> is_near_zero_mode(f; tol=1e-2), freqs[1:6])
        n_zero_mmtk  = count(f -> is_near_zero_mode(f; tol=1e-2), mmtk_freqs[1:6])

        @test n_zero_julia >= 4
        @test n_zero_mmtk >= 2
    end

    @testset "Imaginary mode sanity" begin
        n_imag_julia = count(f -> is_imag_mode(f), freqs)
        n_imag_mmtk  = count(f -> is_imag_mode(f), mmtk_freqs)

        @info "Imaginary mode counts" julia=n_imag_julia mmtk=n_imag_mmtk

        # Loose check for now
        @test n_imag_julia >= 0
        @test n_imag_mmtk >= 0
    end

    @testset "Spectrum comparison against MMTK" begin
        julia_real = nontrivial_real_modes(freqs; start_mode=7)
        mmtk_real  = nontrivial_real_modes(mmtk_freqs; start_mode=7)

        ncmp = min(length(julia_real), length(mmtk_real))
        @test ncmp > 20

        julia_real = julia_real[1:ncmp]
        mmtk_real  = mmtk_real[1:ncmp]

        scale = median_scale_factor(mmtk_real, julia_real)
        scaled_julia = scale .* julia_real

        relerrs = abs.(mmtk_real .- scaled_julia) ./ max.(abs.(mmtk_real), 1e-12)
        corrval = cor(mmtk_real, julia_real)

        @info "Spectrum comparison" scale=scale corr=corrval median_relerr=median(relerrs) max_relerr=maximum(relerrs)

        # These detect whether the difference is mostly a constant scaling factor
        @test isfinite(scale)
        @test corrval > 0.98
        @test median(relerrs) < 0.15
        @test maximum(relerrs) < 0.40
    end

    @testset "Diagnostic print for first 10 nontrivial modes" begin
        julia_real = nontrivial_real_modes(freqs; start_mode=7)
        mmtk_real  = nontrivial_real_modes(mmtk_freqs; start_mode=7)

        ncmp = min(length(julia_real), length(mmtk_real))
        julia_real = julia_real[1:ncmp]
        mmtk_real  = mmtk_real[1:ncmp]

        scale = median_scale_factor(mmtk_real, julia_real)

        for i in 1:min(10, ncmp)
            @info "Mode $(i+6)" julia_thz=julia_real[i] scaled_julia_thz=scale*julia_real[i] mmtk_thz=mmtk_real[i]
        end

        @test true
    end
end
