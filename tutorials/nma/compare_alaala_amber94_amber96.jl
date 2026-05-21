using Printf
using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA: extract_coordinates_masses,
                                 build_hessian,
                                 compute_modes,
                                 frequencies_thz

function run_alaala_nma(label::AbstractString, ff_path::Union{Nothing, String})
    pdb_path = ball_data_path("../test/data/AlaAla.pdb")
    mol = load_pdb(pdb_path, Float64)

    fdb = FragmentDB{Float64}()
    normalize_names!(mol, fdb)
    reconstruct_fragments!(mol, fdb)
    build_bonds!(mol, fdb)

    ff = isnothing(ff_path) ? AmberFF(mol) : AmberFF(mol, ff_path)
    update!(ff)
    initial_energy = compute_energy!(ff)
    optimize_structure!(ff)
    update!(ff)
    minimized_energy = compute_energy!(ff)

    coords, masses = extract_coordinates_masses(mol)
    H = build_hessian(coords)
    vals, vecs, inv_sqrt_M = compute_modes(H, masses)
    freqs = frequencies_thz(vals, masses; unit=:thz)

    return (;
        label,
        initial_energy,
        minimized_energy,
        coords,
        masses,
        vals,
        vecs,
        inv_sqrt_M,
        freqs,
    )
end

function real_frequency(freq::ComplexF64)
    return iszero(imag(freq)) ? real(freq) : NaN
end

function write_frequency_file(path::AbstractString, result)
    open(path, "w") do io
        println(io, "# Mode\tEigenvalue\tFrequency_THz")
        for i in eachindex(result.vals)
            f = result.freqs[i]
            ftxt = iszero(imag(f)) ? string(real(f)) : string(imag(f), "im")
            println(io, i, '\t', result.vals[i], '\t', ftxt)
        end
    end
    return path
end

function main()
    amber94 = ball_data_path("forcefields/AMBER/amber94.ini")

    println("Running AlaAla AMBER96/default NMA...")
    amber96_result = run_alaala_nma("AMBER96/default", nothing)

    println("Running AlaAla AMBER94 NMA...")
    amber94_result = run_alaala_nma("AMBER94", amber94)

    n_atoms96 = size(amber96_result.coords, 1)
    n_atoms94 = size(amber94_result.coords, 1)
    n_modes96 = length(amber96_result.vals)
    n_modes94 = length(amber94_result.vals)

    outdir = joinpath(@__DIR__, "outputs", "alaala_amber_compare")
    mkpath(outdir)
    write_frequency_file(joinpath(outdir, "alaala_amber96_frequencies.txt"), amber96_result)
    write_frequency_file(joinpath(outdir, "alaala_amber94_frequencies.txt"), amber94_result)

    println()
    println("AlaAla AMBER94 vs AMBER96/default comparison")
    println("============================================")
    println("AMBER94 parameter file: ", amber94)
    println()
    @printf("%-18s %10s %10s %22s %22s\n", "Force field", "Atoms", "Modes", "Initial E (kJ/mol)", "Minimized E (kJ/mol)")
    @printf("%-18s %10d %10d %22.10f %22.10f\n",
            amber96_result.label, n_atoms96, n_modes96,
            amber96_result.initial_energy, amber96_result.minimized_energy)
    @printf("%-18s %10d %10d %22.10f %22.10f\n",
            amber94_result.label, n_atoms94, n_modes94,
            amber94_result.initial_energy, amber94_result.minimized_energy)

    println()
    println("First non-trivial modes (7-12)")
    println("==============================")
    @printf("%6s %18s %18s %18s %14s\n", "Mode", "AMBER96 THz", "AMBER94 THz", "Delta THz", "% change")
    for i in 7:12
        f96 = real_frequency(amber96_result.freqs[i])
        f94 = real_frequency(amber94_result.freqs[i])
        delta = f94 - f96
        pct = 100 * delta / f96
        @printf("%6d %18.12f %18.12f %18.12e %13.6f\n", i, f96, f94, delta, pct)
    end

    println()
    println("All frequencies written to: ", outdir)
end

main()

