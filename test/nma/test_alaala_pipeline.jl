using Test
using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA: extract_coordinates_masses, build_hessian, compute_modes

@testset "AlaAla full NMA pipeline" begin
    pdb_path = ball_data_path("../test/data/AlaAla.pdb")
    mol = load_pdb(pdb_path, Float64)

    fdb = FragmentDB{Float64}()
    normalize_names!(mol, fdb)
    reconstruct_fragments!(mol, fdb)
    build_bonds!(mol, fdb)

    ff = AmberFF(mol)
    update!(ff)
    optimize_structure!(ff)
    @show compute_energy!(ff)
    
    coords, masses = extract_coordinates_masses(mol)
    H = build_hessian(coords)
    vals, vecs, invM = compute_modes(H, masses)

    N = size(coords, 1)

    @test compute_energy!(ff) ≈ -374.3136 atol=1e-2
    @test size(H) == (3N, 3N)   #TEST the function build_hessian returns the correct size Hessian
    @test length(vals) == 3N    #TEST the number of eigenvalues matches the expected 3N modes for N atoms  #69 Modes for AlaAla (23 atoms)
    @test size(vecs) == (3N, 3N)  #TEST FUNCTION compute_modes returns the correct size eigenvector matrix

    @test all(abs.(vals[1:6]) .< 1e-5)  #TEST the first 6 eigenvalues are approximately zero (rigid body modes)
    @test all(vals[7:end] .> 0)    #TEST all non-rigid body modes have positive eigenvalues (stable vibrations)

    @test issorted(vals)   #TEST the eigenvalues are sorted in ascending order
end