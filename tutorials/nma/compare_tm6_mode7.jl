using Printf

const OUTDIR = joinpath(@__DIR__, "outputs", "cb1_receptor_nma")

struct AtomRecord
    name::String
    resname::String
    chain::Char
    resi::Int
    element::String
end

function pdb_field(line::AbstractString, first::Int, last::Int)
    length(line) < first && return ""
    return line[first:min(last, lastindex(line))]
end

function read_heavy_atoms(pdb_file::AbstractString)
    atoms = AtomRecord[]
    for line in eachline(pdb_file)
        startswith(line, "ATOM") || continue
        name = strip(pdb_field(line, 13, 16))
        resname = strip(pdb_field(line, 18, 20))
        chain = only(pdb_field(line, 22, 22))
        resi = parse(Int, strip(pdb_field(line, 23, 26)))
        element = strip(pdb_field(line, 77, 78))
        element == "H" && continue
        push!(atoms, AtomRecord(name, resname, chain, resi, element))
    end
    return atoms
end

function read_xyz_frames(xyz_file::AbstractString)
    lines = readlines(xyz_file)
    n_atoms = parse(Int, strip(lines[1]))
    frame_size = n_atoms + 2
    n_frames = div(length(lines), frame_size)
    frames = Vector{Matrix{Float64}}(undef, n_frames)

    for frame in 1:n_frames
        offset = (frame - 1) * frame_size
        coords = zeros(Float64, n_atoms, 3)
        for i in 1:n_atoms
            parts = split(strip(lines[offset + 2 + i]))
            coords[i, 1] = parse(Float64, parts[2])
            coords[i, 2] = parse(Float64, parts[3])
            coords[i, 3] = parse(Float64, parts[4])
        end
        frames[frame] = coords
    end

    return frames
end

function mode_displacements_from_xyz(xyz_file::AbstractString; frame_a::Int=8, frame_ref::Int=15)
    frames = read_xyz_frames(xyz_file)
    scale_a = 2.0 * sin(2π * frame_a / length(frames))
    scale_ref = 2.0 * sin(2π * frame_ref / length(frames))
    scale_delta = scale_a - scale_ref
    return (frames[frame_a] .- frames[frame_ref]) ./ scale_delta
end

function summarize_region(label, pdb_name, xyz_name; chain::Char='A', first_res::Int=332, last_res::Int=368)
    pdb_file = joinpath(OUTDIR, pdb_name)
    xyz_file = joinpath(OUTDIR, xyz_name)

    atoms = read_heavy_atoms(pdb_file)
    disp = mode_displacements_from_xyz(xyz_file)

    length(atoms) == size(disp, 1) ||
        error("Atom count mismatch for $label: $(length(atoms)) atoms in PDB, $(size(disp, 1)) atoms in XYZ")

    all_magnitudes = [sqrt(sum(abs2, disp[i, :])) for i in eachindex(atoms)]
    region_indices = findall(a -> a.chain == chain && first_res <= a.resi <= last_res, atoms)
    isempty(region_indices) && error("No atoms found for $label region $chain:$first_res-$last_res")

    region_magnitudes = all_magnitudes[region_indices]

    return (
        label = label,
        total_atoms = length(atoms),
        region_atoms = length(region_indices),
        receptor_mean = sum(all_magnitudes) / length(all_magnitudes),
        tm6_mean = sum(region_magnitudes) / length(region_magnitudes),
        tm6_max = maximum(region_magnitudes),
    )
end

wt = summarize_region(
    "WT 5TGZ",
    "5TGZ_minimized_cb1.pdb",
    "5TGZ_heavy_mode_7.xyz",
)

mut = summarize_region(
    "F238L repacked",
    "5TGZ_F238L_repacked_minimized_cb1.pdb",
    "5TGZ_F238L_repacked_heavy_mode_7.xyz",
)

println("Mode 7 displacement comparison for TM6 (chain A residues 332-368)")
println("Displacements are from the normalized exported mode vector.")
println()
@printf("%-16s %12s %12s %18s %18s %18s\n",
        "System", "Atoms", "TM6 atoms", "Mean receptor", "Mean TM6", "Max TM6")
for row in (wt, mut)
    @printf("%-16s %12d %12d %18.10f %18.10f %18.10f\n",
            row.label, row.total_atoms, row.region_atoms,
            row.receptor_mean, row.tm6_mean, row.tm6_max)
end

delta = mut.tm6_mean - wt.tm6_mean
percent = 100 * delta / wt.tm6_mean
println()
@printf("TM6 mean displacement change: %.10f (%+.4f%%)\n", delta, percent)

