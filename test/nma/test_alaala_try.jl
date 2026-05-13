using BiochemicalAlgorithms.NMA
using LinearAlgebra
using Statistics
using BiochemicalAlgorithms.NMA: extract_coordinates_masses

println("Running AlaAla NMA pipeline...")

pdb_path = ball_data_path("../test/data/AlaAla.pdb")
mol = load_pdb(pdb_path, Float64)

fdb = FragmentDB{Float64}()
normalize_names!(mol, fdb)
reconstruct_fragments!(mol, fdb)
build_bonds!(mol, fdb)

ff = AmberFF(mol)
update!(ff)

@show compute_energy!(ff)

optimize_structure!(ff)

@show compute_energy!(ff)
coords, masses = BiochemicalAlgorithms.NMA.extract_coordinates_masses(mol)

H = BiochemicalAlgorithms.NMA.build_hessian(coords)
vals, vecs, invM = BiochemicalAlgorithms.NMA.compute_modes(H, masses)

println("Eigenvalues:")
println(vals[1:20])

println("Done.")