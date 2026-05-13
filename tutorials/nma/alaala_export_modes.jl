using Pkg
Pkg.add("DelimitedFiles")


using LinearAlgebra
using Plots
using DelimitedFiles
using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA: extract_coordinates_masses,
                                 build_hessian,
                                 compute_modes,
                                 compute_frequencies, 
                                 vmd_mode_visuals,
                                 energy_along_mode!
                                 frequencies_thz
                                 plot_energy_displacement


println("Running AlaAla NMA export script...")
println("========================================")
println("   AlaAla NMA Analysis & Visualization   ")

# 1. Load and preprocess system
pdb_path = ball_data_path("../test/data/AlaAla.pdb")
mol = load_pdb(pdb_path, Float64)

fdb = FragmentDB{Float64}()
normalize_names!(mol, fdb)
reconstruct_fragments!(mol, fdb)
build_bonds!(mol, fdb)

println("System loaded and preprocessed.")

# 2. Force field and minimization
ff = AmberFF(mol)
update!(ff)

println("Initial energy: ", compute_energy!(ff), " kJ/mol")

optimize_structure!(ff)

E_min = compute_energy!(ff)
println("Energy after minimization: ", E_min, " kJ/mol")

# 3. Extract coordinates & masses
coords, masses = extract_coordinates_masses(mol)
N_atoms = size(coords, 1)
total_modes = 3 * N_atoms

println("Atoms: $N_atoms | Modes: $total_modes")
println("Total vibrational modes: ", total_modes)

# 4. Build Hessian and compute modes
H = build_hessian(coords)
vals, vecs, inv_sqrt_M = compute_modes(H, masses)

# 5. Compute frequencies (optional, for interpretation)
freqs = frequencies_thz(vals, masses; unit=:thz)

println("First 12 Eigenvalues and Frequencies")
println("====================================")

for i in 1:12
    λ = vals[i]
    f = freqs[i]

    if imag(f) != 0
        println("Mode $i : λ = $λ    freq = $(imag(f))im THz")
    else
        println("Mode $i : λ = $λ    freq = $(real(f)) THz")
    end
end

# Output directory
outdir = joinpath(@__DIR__, "outputs", "alaala_nma")
mkpath(outdir)

eigfreq_file = joinpath(outdir, "all_eigenvalues_frequencies_thz.txt")

open(eigfreq_file, "w") do io
    println(io, "# Mode\tEigenvalue\tFrequency_THz")
    for i in eachindex(vals)
        λ = vals[i]
        f = freqs[i]

        if imag(f) != 0
            println(io, "$i\t$λ\t$(imag(f))im")
        else
            println(io, "$i\t$λ\t$(real(f))")
        end
    end
end

println("Eigenvalues + THz frequencies saved → $eigfreq_file")



# Optional plain numeric export too
#writedlm(joinpath(outdir, "all_eigenvalues_raw.txt"), vals)

# --- Export VMD modes ---

println("\nExporting modes for VMD...")
for m in 7:10
    filename = joinpath(outdir, "mode_$(m).xyz")
    vmd_mode_visuals(filename, coords, vecs, m, inv_sqrt_M)
    println("Exported mode $m → $filename")
end

# --- Plot energy vs displacement ---
println("\nGenerating energy landscape plot...")

plotfile = joinpath(outdir, "energy_vs_displacement_modes_7_12.png")

p = plot_energy_displacement(
    ff,
    coords,
    vecs,
    inv_sqrt_M,
    7:12;
    filename = plotfile
)

display(p)

println("Plot saved → $plotfile")
# 8. Done
println("   NMA analysis completed successfully   ")
println("Done.")