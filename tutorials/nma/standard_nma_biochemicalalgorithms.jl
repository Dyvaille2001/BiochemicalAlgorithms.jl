using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA

# Standard normal mode calculation in BiochemicalAlgorithms.jl.
# This follows the MMTK workflow:
# load structure -> preprocess -> build Amber force field -> minimize -> NMA.

pdb_path = ball_data_path("../test/data/bALA1.pdb")

mol = load_pdb(pdb_path, Float64)
fragments = FragmentDB{Float64}()
normalize_names!(mol, fragments)
reconstruct_fragments!(mol, fragments)
build_bonds!(mol, fragments)

ff = AmberFF(mol)
println("Initial energy: ", compute_energy!(ff), " kJ/mol")

optimize_structure!(ff)
println("Minimized energy: ", compute_energy!(ff), " kJ/mol")

nma = normal_mode_analysis(mol; cutoff=12.0, gamma=1.0, unit=:thz)

println()
println("Mode    eigenvalue              frequency (THz)")
for i in eachindex(nma.eigenvalues)
    freq = nma.frequencies[i]
    label = iszero(imag(freq)) ? string(real(freq)) : string(imag(freq), "im")
    println(rpad(i, 7), rpad(nma.eigenvalues[i], 24), label)
end

mode_index = first_vibrational_mode(nma)
println()
println("First non-rigid vibrational mode is mode ", mode_index)
println("Frequency: ", nma.frequencies[mode_index], " THz")

outdir = joinpath(@__DIR__, "outputs", "standard_nma")
mkpath(outdir)

mode_file = joinpath(outdir, "mode_$(mode_index).xyz")
write_vmd_mode(
    mode_file,
    nma.coords,
    nma.mass_weighted_modes,
    mode_index,
    nma.inv_sqrt_mass_matrix;
    n_frames=30,
)

println("Mode animation written to: ", mode_file)
