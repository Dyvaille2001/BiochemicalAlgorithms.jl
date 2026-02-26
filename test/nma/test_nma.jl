@testitem "TrivialAtomBijection" begin
    for T in [Float32, Float64]
        mol = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
        fdb = FragmentDB()
        normalize_names!(mol, fdb)
        reconstruct_fragments!(mol, fdb)
        build_bonds!(mol, fdb)
        ff = AmberFF(mol)

        @test compute_energy!(ff) 
        optimize_structure!(ff)
        @test compute_energy!(ff) ≈ -374.3136f0

        
    end
end