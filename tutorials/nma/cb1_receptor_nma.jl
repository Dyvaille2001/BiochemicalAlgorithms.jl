using BiochemicalAlgorithms
using BiochemicalAlgorithms.NMA

function usage()
    println("""
    Usage:
        julia --project=. tutorials/nma/cb1_receptor_nma.jl <CB1.pdb> [selection] [--no-minimize] [--maxiters=N]

    selection:
        ca     C-alpha atoms only (default, recommended for receptor-scale ANM)
        heavy  all non-hydrogen protein atoms
        all    all atoms loaded from the PDB

    Example:
        julia --project=. tutorials/nma/cb1_receptor_nma.jl data/cb1/CB1.pdb ca

    By default, the selected CB1 receptor part is preprocessed and minimized
    with AmberFF before NMA. Use --no-minimize only for quick debugging.
    """)
end

function parse_maxiters(args; default=500)
    for arg in args
        if startswith(arg, "--maxiters=")
            return parse(Int, split(arg, "=", limit=2)[2])
        end
    end
    return default
end

function cb1_ranges_for_pdb(pdb_path::AbstractString)
    pdb_id = uppercase(splitext(basename(pdb_path))[1])
    startswith(pdb_id, "5TGZ") && return [(99, 306), (332, 414)]
    startswith(pdb_id, "5U09") && return [(90, 301), (333, 421)]
    startswith(pdb_id, "6KPG") && return [(71, 425)]
    return [(typemin(Int), typemax(Int))]
end

function cb1_receptor_atom(atom, cb1_ranges)
    frag = parent_fragment(atom)
    isnothing(frag) && return false
    is_amino_acid(frag) || return false
    return any(first(range) <= frag.number <= last(range) for range in cb1_ranges)
end

function select_atoms_for_nma(mol, selection::AbstractString, cb1_ranges)
    all_atoms = atoms(mol)

    if selection == "ca"
        return filter(atom -> cb1_receptor_atom(atom, cb1_ranges) && strip(atom.name) == "CA", all_atoms)
    elseif selection == "heavy"
        return filter(atom -> cb1_receptor_atom(atom, cb1_ranges) && atom.element != Elements.H, all_atoms)
    elseif selection == "all"
        return filter(atom -> cb1_receptor_atom(atom, cb1_ranges), all_atoms)
    else
        throw(ArgumentError("Unknown selection '$selection'. Use ca, heavy, or all."))
    end
end

function load_pdb_coordinates_for_nma(pdb_path::AbstractString)
    cleaned_path = tempname() * ".pdb"
    open(cleaned_path, "w") do io
        for line in eachline(pdb_path)
            record = length(line) >= 6 ? line[1:6] : line
            if record ∉ ("HELIX ", "SHEET ")
                println(io, line)
            end
        end
    end

    try
        return load_pdb(cleaned_path, Float64; keep_metadata=false, create_coils=false)
    finally
        isfile(cleaned_path) && rm(cleaned_path; force=true)
    end
end

function keep_cb1_receptor_only!(mol, cb1_ranges)
    for frag in collect(fragments(mol))
        keep = is_amino_acid(frag) &&
               any(first(range) <= frag.number <= last(range) for range in cb1_ranges)
        keep || delete!(frag)
    end
    for chain in collect(chains(mol))
        length(fragments(chain)) == 0 && delete!(chain)
    end
    return mol
end

function preprocess_and_minimize!(mol, cb1_ranges; do_minimize::Bool=true, maxiters::Int=500)
    fdb = FragmentDB{Float64}()
    println("Keeping CB1 receptor fragments only...")
    keep_cb1_receptor_only!(mol, cb1_ranges)
    println("Normalizing residue and atom names...")
    normalize_names!(mol, fdb)
    println("Reconstructing fragments and building bonds...")
    reconstruct_fragments!(mol, fdb)
    build_bonds!(mol, fdb)

    if do_minimize
        println("Building AmberFF and minimizing CB1 receptor coordinates...")
        ff = AmberFF(mol)
        update!(ff)
        initial_energy = compute_energy!(ff)
        println("Initial potential energy: ", initial_energy, " kJ/mol")
        optimize_structure!(ff; maxiters=maxiters)
        update!(ff)
        minimized_energy = compute_energy!(ff)
        println("Minimized potential energy: ", minimized_energy, " kJ/mol")
        @show compute_energy!(ff)
    else
        println("Skipping AmberFF minimization (--no-minimize).")
    end

    return mol
end

if isempty(ARGS)
    usage()
    error("Missing CB1 receptor PDB path.")
end

pdb_path = abspath(ARGS[1])
extra_args = ARGS[2:end]
do_minimize = !("--no-minimize" in extra_args)
maxiters = parse_maxiters(extra_args)
selection_args = filter(arg -> arg != "--no-minimize" && !startswith(arg, "--maxiters="), extra_args)
selection = isempty(selection_args) ? "ca" : lowercase(first(selection_args))

isfile(pdb_path) || error("CB1 receptor PDB file not found: $pdb_path")

println("Running CB1 receptor NMA")
println("========================")
println("Input:     ", pdb_path)
println("Selection: ", selection)
println("Minimize:  ", do_minimize)
println("Max iters: ", maxiters)
output_prefix = splitext(basename(pdb_path))[1]
cb1_ranges = cb1_ranges_for_pdb(pdb_path)
println("CB1 ranges: ", cb1_ranges)

outdir = joinpath(@__DIR__, "outputs", "cb1_receptor_nma")
mkpath(outdir)

mol = load_pdb_coordinates_for_nma(pdb_path)
preprocess_and_minimize!(mol, cb1_ranges; do_minimize=do_minimize, maxiters=maxiters)

coordinate_label = do_minimize ? "minimized" : "preprocessed"
minimized_pdb = joinpath(outdir, "$(output_prefix)_$(coordinate_label)_cb1.pdb")
write_pdb(minimized_pdb, mol)
println("CB1 receptor PDB used for NMA written to: ", minimized_pdb)

nma_atoms = select_atoms_for_nma(mol, selection, cb1_ranges)
n_selected = length(nma_atoms)
n_selected > 0 || error("No atoms selected for NMA. Try selection='heavy' or check the PDB contents.")

println("Selected atoms: ", n_selected)
println("Modes:          ", 3 * n_selected)

nma = normal_mode_analysis(nma_atoms; cutoff=12.0, gamma=1.0, unit=:thz)

freq_file = joinpath(outdir, "$(output_prefix)_$(selection)_eigenvalues_frequencies.txt")
open(freq_file, "w") do io
    println(io, "# mode\teigenvalue\tfrequency_thz")
    for i in eachindex(nma.eigenvalues)
        freq = nma.frequencies[i]
        freq_text = iszero(imag(freq)) ? string(real(freq)) : string(imag(freq), "im")
        println(io, i, '\t', nma.eigenvalues[i], '\t', freq_text)
    end
end

println("Eigenvalues and frequencies written to: ", freq_file)

first_mode = first_vibrational_mode(nma)

if length(nma.eigenvalues) >= first_mode
    last_mode = min(first_mode + 3, length(nma.eigenvalues))

    println("Exporting mode animations: ", first_mode, ":", last_mode)
    for mode_index in first_mode:last_mode
        mode_file = joinpath(outdir, "$(output_prefix)_$(selection)_mode_$(mode_index).xyz")
        write_vmd_mode(
            mode_file,
            nma.coords,
            nma.mass_weighted_modes,
            mode_index,
            nma.inv_sqrt_mass_matrix;
            n_frames=30,
        )
        println("  mode ", mode_index, " -> ", mode_file)
    end

    println()
    println("First non-rigid vibrational mode: ", first_mode)
    println("Frequency: ", nma.frequencies[first_mode], " THz")
else
    println("Skipping mode animation: selected structure has fewer than 7 modes.")
end

println("Done.")
