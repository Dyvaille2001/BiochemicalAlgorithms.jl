#Comparing computed AlaAla modes to MMTK reference modes using a mode-overlap heatmap.
#This script reads a cleaned MMTK mode export, computes modes for AlaAla using
#BiochemicalAlgorithms.jl, and generates a heatmap of mode overlaps between the two sets.


using Statistics
using LinearAlgebra
using Plots

using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA: extract_coordinates_masses,
                                 build_hessian,
                                 compute_modes,
                                 frequencies_thz

# ============================================================
# Read clean MMTK mode export
# Format per line:
# mode_index frequency v1 v2 v3 ... v(3N)
# ============================================================

function read_mmtk_modes_clean(filepath::AbstractString)
    mode_indices = Int[]
    freqs = ComplexF64[]
    modes = Vector{Vector{Float64}}()

    for line in eachline(filepath)
        s = strip(line)
        isempty(s) && continue

        parts = split(s)
        length(parts) >= 3 || error("Invalid line in MMTK modes file: $line")

        mode_idx = parse(Int, parts[1])
        freq_str = parts[2]

        freq = if endswith(freq_str, "j")
            complex(0.0, parse(Float64, freq_str[1:end-1]))
        else
            complex(parse(Float64, freq_str), 0.0)
        end

        vec = parse.(Float64, parts[3:end])

        push!(mode_indices, mode_idx)
        push!(freqs, freq)
        push!(modes, vec)
    end

    return (; mode_indices, freqs, modes)
end

# ============================================================
# Compute Julia AlaAla modes
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
# Helpers
# ============================================================

is_imag_mode(z::ComplexF64; tol=1e-12) = abs(imag(z)) > tol && abs(real(z)) ≤ tol
is_real_mode(z::ComplexF64; tol=1e-12) = abs(real(z)) > tol && abs(imag(z)) ≤ tol

function unweight_modes(vecs::AbstractMatrix, inv_sqrt_M::AbstractMatrix)
    modes = Vector{Vector{Float64}}(undef, size(vecs, 2))
    for i in 1:size(vecs, 2)
        modes[i] = Vector{Float64}(inv_sqrt_M * vecs[:, i])
    end
    return modes
end

function mode_overlap(u::Vector{Float64}, v::Vector{Float64})
    nu = norm(u)
    nv = norm(v)
    (nu == 0 || nv == 0) && return 0.0
    return abs(dot(u, v)) / (nu * nv)
end

function nontrivial_real_modes_with_vectors(freqs::Vector{ComplexF64},
                                            modes::Vector{Vector{Float64}};
                                            start_mode::Int=7)
    kept_freqs = Float64[]
    kept_modes = Vector{Vector{Float64}}()
    kept_indices = Int[]

    for i in start_mode:length(freqs)
        f = freqs[i]
        if is_real_mode(f) && real(f) > 0
            push!(kept_freqs, real(f))
            push!(kept_modes, modes[i])
            push!(kept_indices, i)
        end
    end

    return kept_freqs, kept_modes, kept_indices
end

function overlap_matrix(julia_modes::Vector{Vector{Float64}},
                        mmtk_modes::Vector{Vector{Float64}})
    M = Matrix{Float64}(undef, length(julia_modes), length(mmtk_modes))
    for i in eachindex(julia_modes)
        for j in eachindex(mmtk_modes)
            M[i, j] = mode_overlap(julia_modes[i], mmtk_modes[j])
        end
    end
    return M
end

# ============================================================
# Main heatmap generation
# ============================================================

function generate_mode_overlap_heatmap(; nmodes::Int=20)
    pkgroot = pkgdir(BiochemicalAlgorithms)
    outdir = joinpath(pkgroot, "tutorials", "nma", "outputs", "alaala_nma")

    mmtk_file = joinpath(outdir, "ala_mmtk_modes_clean.txt")
    isfile(mmtk_file) || error("MMTK mode file not found: $mmtk_file")

    mmtk_data = read_mmtk_modes_clean(mmtk_file)
    computed = compute_alaala_modes()

    julia_modes_all = unweight_modes(computed.vecs, computed.inv_sqrt_M)

    julia_freqs, julia_modes, julia_indices =
        nontrivial_real_modes_with_vectors(computed.freqs, julia_modes_all; start_mode=7)

    mmtk_freqs, mmtk_modes, mmtk_indices =
        nontrivial_real_modes_with_vectors(mmtk_data.freqs, mmtk_data.modes; start_mode=7)

    n = min(nmodes, length(julia_modes), length(mmtk_modes))

    julia_freqs = julia_freqs[1:n]
    julia_modes = julia_modes[1:n]
    julia_indices = julia_indices[1:n]

    mmtk_freqs = mmtk_freqs[1:n]
    mmtk_modes = mmtk_modes[1:n]
    mmtk_indices = mmtk_indices[1:n]

    O = overlap_matrix(julia_modes, mmtk_modes)

    xlabels = ["$(i)\n$(round(f, digits=2))" for (i, f) in zip(mmtk_indices, mmtk_freqs)]
    ylabels = ["$(i)\n$(round(f, digits=2))" for (i, f) in zip(julia_indices, julia_freqs)]

    p = heatmap(
        1:n, 1:n, O;
        xlabel = "MMTK mode index / THz",
        ylabel = "Julia mode index / THz",
        title = "Mode-overlap heatmap (first $n nontrivial modes)",
        xticks = (1:n, xlabels),
        yticks = (1:n, ylabels),
        colorbar_title = "overlap",
        clim = (0.0, 1.0),
        xrotation = 45,
        size = (1200, 900)
    )

    pngfile = joinpath(outdir, "mode_overlap_heatmap_first_$(n)_nontrivial_modes.png")
    savefig(p, pngfile)

    println("Saved heatmap to: $pngfile")

    # Print strongest match in each Julia row
    println("\nBest MMTK match for each Julia mode:")
    for i in 1:n
        bestval, bestj = findmax(O[i, :])
        println("Julia mode $(julia_indices[i]) (THz=$(round(julia_freqs[i], digits=3))) " *
                "-> MMTK mode $(mmtk_indices[bestj]) (THz=$(round(mmtk_freqs[bestj], digits=3))) " *
                "overlap=$(round(bestval, digits=3))")
    end

    return p, O
end

# Run
p, O = generate_mode_overlap_heatmap(nmodes=20)
display(p)