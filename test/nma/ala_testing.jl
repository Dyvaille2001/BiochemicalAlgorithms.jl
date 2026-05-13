@testitem "AlaAla full NMA" begin
    using BiochemicalAlgorithms
    using LinearAlgebra
    using Statistics

    # --- Load system (USE EXISTING PACKAGE FUNCTIONS) ---
    pdb_path = ball_data_path("../test/data/AlaAla.pdb")
    mol = load_pdb(pdb_path, Float64)

    fdb = FragmentDB()
    normalize_names!(mol, fdb)
    reconstruct_fragments!(mol, fdb)
    build_bonds!(mol, fdb)

    # --- Force field ---
    ff = AmberFF(mol)
    update!(ff)

    # --- Minimize ---
    optimize_structure!(ff)

    # --- Your NMA functions ---
    coords, masses = extract_coordinates_masses(mol)

    H = build_hessian(coords)

    vals, vecs, invM = compute_modes(H, masses)

    # --- Tests ---
    N = size(coords,1)

    @test size(H) == (3N, 3N)
    @test length(vals) == 3N

    # 6 rigid body modes
    @test all(abs.(vals[1:6]) .< 1e-5)

    # Physical modes positive
    @test all(vals[7:end] .> 0)

    # B-factors
    #b = compute_b_factors(vals, vecs, invM, N)

    @test length(b) == N
    @test all(b .>= 0)
    @test std(b) > 0
end