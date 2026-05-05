"""
    extract_coords_masses(sys)

Extract Cartesian coordinates and atomic masses from a molecular system.

This function iterates through the atoms in the provided system object, 
mapping element symbols to their respective atomic weights and 
collecting their spatial positions.

Parameters
----------
sys : System
    A molecular system object (e.g., from BiochemicalAlgorithms.jl, Molly.jl or BioStructures.jl) 
    that supports the `atoms(sys)` and `atom.r` interface.

Returns
-------
coords : Matrix{Float64}
    N×3 matrix of atomic positions in Ångströms.
masses : Vector{Float64}
    Vector of length N containing atomic masses in atomic mass units (amu).

Notes
-----
- Unknown elements default to the mass of Carbon (12.011 amu).
- The order of atoms in the output matrices matches the order in the system.
"""
const ATOMIC_MASSES = Dict(
    "H" => 1.008, "C" => 12.011, "N" => 14.007,
    "O" => 15.999, "S" => 32.06
    )


function extract_coordinates_masses(sys)
    atoms_list = collect(atoms(sys))
    N = length(atoms_list)

    coords = zeros(N, 3)
    masses = zeros(N)

    
    for (i, atom) in enumerate(atoms_list)
        coords[i, :] = atom.r
        masses[i] = get(ATOMIC_MASSES, string(atom.element), 12.011)
    end

    return coords, masses
end