using Pkg
Pkg.add("DelimitedFiles")


using LinearAlgebra
using Plots
using DelimitedFiles
using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA: extract_coordinates_masses,
                                 build_hessian,
                                 compute_modes,
                                 vmd_mode_visuals,
                                 energy_along_mode!
                                 plot_energy_displacement


println("Running bala1 NMA export script...")
println("========================================")
println("   bala1 NMA Analysis & Visualization   ")

# 1. Load and preprocess system
pdb_path = ball_data_path("../test/data/bALA1.pdb")
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

println("First 12 eigenvalues:")
println(vals[1:12])


# Output directory
outdir = joinpath(@__DIR__, "outputs", "bala1_nma")
mkpath(outdir)

# 6. Save all eigenvalues

eigenvalues_file = joinpath(outdir, "bala1_eigenvalues.txt")
open(eigenvalues_file, "w") do io
    println(io, "# Mode\tEigenvalue")
    for (i, λ) in enumerate(vals)
        println(io, "$i\t$λ")
    end
end

println("All eigenvalues saved → $eigenvalues_file")

# Optional plain numeric export too
#writedlm(joinpath(outdir, "bala1_eigenvalues_raw.txt"), vals)

# --- Export VMD modes ---

println("\nExporting modes for VMD...")
for m in 7:10
    filename = joinpath(outdir, "mode_$(m).xyz")
    vmd_mode_visuals(filename, coords, vecs, m, inv_sqrt_M)
    println("Exported mode $m → $filename")
end

# --- Plot energy vs displacement ---
println("\nGenerating energy landscape plot...")

plotfile = joinpath(outdir, "energy_vs_displacement_modes_bala1_7_12.png")

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