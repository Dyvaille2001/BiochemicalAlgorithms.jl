const INPUT_PDB = abspath(joinpath(@__DIR__, "..", "..", "test", "data", "5TGZ.pdb"))
const OUTPUT_PDB = abspath(joinpath(@__DIR__, "..", "..", "test", "data", "5TGZ_F238L.pdb"))

const TARGET_CHAIN = 'A'
const TARGET_RESI = 237
const SOURCE_RESNAME = "PHE"
const MUTANT_RESNAME = "LEU"

# In 5TGZ/P21554 numbering, the structural residue corresponding to the
# custom-numbered F238 site is PHE A 237. This script makes a lightweight
# PHE->LEU model for ANM/NMA by retaining the common backbone/CB atoms, mapping
# PHE CG/CD1/CD2 onto LEU CG/CD1/CD2, and deleting the extra aromatic atoms.
const PHE_TO_LEU_ATOMS = Dict(
    "N" => "N",
    "CA" => "CA",
    "C" => "C",
    "O" => "O",
    "CB" => "CB",
    "CG" => "CG",
    "CD1" => "CD1",
    "CD2" => "CD2",
)

function pdb_field(line::AbstractString, first::Int, last::Int)
    length(line) < first && return ""
    return line[first:min(last, lastindex(line))]
end

function replace_pdb_field(line::String, first::Int, last::Int, value::AbstractString)
    padded = rpad(line, max(last, length(line)))
    return padded[1:first-1] * value * padded[last+1:end]
end

function mutate_atom_line(line::String)
    record = pdb_field(line, 1, 6)
    record in ("ATOM  ", "HETATM", "ANISOU") || return line

    resname = strip(pdb_field(line, 18, 20))
    chain = only(pdb_field(line, 22, 22))
    resi = parse(Int, strip(pdb_field(line, 23, 26)))

    if resname != SOURCE_RESNAME || chain != TARGET_CHAIN || resi != TARGET_RESI
        return line
    end

    atom_name = strip(pdb_field(line, 13, 16))
    haskey(PHE_TO_LEU_ATOMS, atom_name) || return nothing

    new_atom_name = PHE_TO_LEU_ATOMS[atom_name]
    mutated = replace_pdb_field(line, 13, 16, lpad(new_atom_name, 4))
    mutated = replace_pdb_field(mutated, 18, 20, MUTANT_RESNAME)
    return mutated
end

function main()
    isfile(INPUT_PDB) || error("Input PDB not found: $INPUT_PDB")

    mutated_atoms = 0
    skipped_atoms = 0

    open(OUTPUT_PDB, "w") do io
        for line in eachline(INPUT_PDB)
            mutated = mutate_atom_line(line)
            if isnothing(mutated)
                skipped_atoms += 1
                continue
            end
            mutated != line && (mutated_atoms += 1)
            println(io, mutated)
        end
    end

    println("Wrote mutant PDB: ", OUTPUT_PDB)
    println("Mutated retained atoms: ", mutated_atoms)
    println("Deleted aromatic-only PHE atoms: ", skipped_atoms)
    println("Mutation modeled as: PHE A 237 -> LEU A 237")
end

main()
