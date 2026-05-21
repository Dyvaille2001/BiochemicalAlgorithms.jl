module MMTKLikeForceFieldNMA

using LinearAlgebra
using Printf

using BiochemicalAlgorithms

export MMTKLikeModes,
       load_prepared_system,
       forcefield_gradient!,
       finite_difference_force_constant_matrix!,
       mmtk_like_nma,
       write_frequency_table,
       write_mode_table

const DEFAULT_MASSES = Dict(
    "H" => 1.008,
    "C" => 12.011,
    "N" => 14.007,
    "O" => 15.999,
    "S" => 32.06,
)

"""
    MMTKLikeModes

Result container for a force-field normal-mode calculation that follows the
same conceptual steps as MMTK `VibrationalModes`:

1. minimize the structure with a molecular mechanics force field,
2. construct the Cartesian force-constant matrix from force-field gradients,
3. mass-weight the force-constant matrix,
4. diagonalize it,
5. transform eigenvectors back to Cartesian mode displacements.

Unlike the ANM helper in `src/nma`, this code does not use a contact-network
Hessian. It numerically differentiates the AMBER force-field gradient.
"""
struct MMTKLikeModes
    coords::Matrix{Float64}
    masses::Vector{Float64}
    force_constants::Matrix{Float64}
    mass_weighted_force_constants::Matrix{Float64}
    eigenvalues::Vector{Float64}
    mass_weighted_modes::Matrix{Float64}
    cartesian_modes::Matrix{Float64}
    frequencies_thz::Vector{ComplexF64}
    initial_energy::Float64
    minimized_energy::Float64
end

function atomic_mass(atom)
    return get(DEFAULT_MASSES, string(atom.element), 12.011)
end

function coordinates_and_masses(sys)
    atom_list = collect(atoms(sys))
    n_atoms = length(atom_list)
    coords = zeros(Float64, n_atoms, 3)
    masses = zeros(Float64, n_atoms)

    for (i, atom) in enumerate(atom_list)
        coords[i, :] .= atom.r
        masses[i] = atomic_mass(atom)
    end

    return coords, masses
end

coordinate_vector(sys) = collect(Float64, Iterators.flatten(atoms(sys).r))

function set_coordinate_vector!(sys, x::AbstractVector{<:Real})
    length(x) == 3 * length(atoms(sys)) ||
        throw(DimensionMismatch("coordinate vector has length $(length(x)); expected $(3 * length(atoms(sys)))"))
    atoms(sys).r .= eachcol(reshape(Float64.(x), 3, :))
    return sys
end

"""
    load_prepared_system(pdb_path; minimize=true, maxiters=100)

Load a PDB, apply the same preprocessing style used in the package tests, build
an AMBER force field, and optionally minimize the coordinates.
"""
function load_prepared_system(pdb_path::AbstractString; minimize::Bool=true, maxiters::Int=100)
    mol = load_pdb(pdb_path, Float64)

    fdb = FragmentDB{Float64}()
    normalize_names!(mol, fdb)
    reconstruct_fragments!(mol, fdb)
    build_bonds!(mol, fdb)

    ff = AmberFF(mol)
    update!(ff)
    initial_energy = Float64(compute_energy!(ff))

    if minimize
        optimize_structure!(ff; maxiters=maxiters)
        update!(ff)
    end

    minimized_energy = Float64(compute_energy!(ff))
    return mol, ff, initial_energy, minimized_energy
end

"""
    forcefield_gradient!(ff, x)

Return the Cartesian potential-energy gradient at coordinate vector `x`.
BiochemicalAlgorithms stores forces, so the gradient is `-force`.
"""
function forcefield_gradient!(ff, x::AbstractVector{<:Real})
    set_coordinate_vector!(ff.system, x)
    update!(ff)
    compute_forces!(ff)
    return -collect(Float64, Iterators.flatten(atoms(ff.system).F))
end

"""
    finite_difference_force_constant_matrix!(ff, x0; delta=1e-4)

Construct the Cartesian force-constant matrix by central finite differences of
the force-field gradient:

`H[:, i] = (g(x0 + δ eᵢ) - g(x0 - δ eᵢ)) / (2δ)`

This is the full-coordinate analogue of the numerical differentiation branch
in MMTK's normal-mode code.
"""
function finite_difference_force_constant_matrix!(ff, x0::AbstractVector{<:Real};
                                                  delta::Float64=1e-4,
                                                  verbose::Bool=true)
    n = length(x0)
    H = zeros(Float64, n, n)
    xplus = copy(Float64.(x0))
    xminus = copy(Float64.(x0))

    for i in 1:n
        xplus[i] += delta
        xminus[i] -= delta

        gplus = forcefield_gradient!(ff, xplus)
        gminus = forcefield_gradient!(ff, xminus)
        H[:, i] .= (gplus .- gminus) ./ (2delta)

        xplus[i] = x0[i]
        xminus[i] = x0[i]

        if verbose && (i == 1 || i == n || i % 25 == 0)
            @printf("finite-difference column %d / %d\n", i, n)
        end
    end

    set_coordinate_vector!(ff.system, x0)
    update!(ff)

    return Symmetric(0.5 .* (H .+ H'))
end

function mass_weight_force_constants(H::AbstractMatrix, masses::AbstractVector)
    masses3 = repeat(Float64.(masses), inner=3)
    inv_sqrt_m = 1.0 ./ sqrt.(masses3)
    Hmw = Diagonal(inv_sqrt_m) * Matrix(H) * Diagonal(inv_sqrt_m)
    return Symmetric(0.5 .* (Hmw .+ Hmw')), inv_sqrt_m
end

function forcefield_frequencies_thz(eigenvalues::AbstractVector)
    NA = 6.02214076e23
    amu = 1.66053906660e-27
    angstrom = 1.0e-10
    internal_to_hz = sqrt((1000.0 / NA) / (amu * angstrom^2)) / (2π)
    factor = internal_to_hz / 1.0e12

    return ComplexF64[
        λ < 0 ? complex(0.0, sqrt(abs(λ)) * factor) : complex(sqrt(λ) * factor, 0.0)
        for λ in eigenvalues
    ]
end

"""
    mmtk_like_nma(pdb_path; delta=1e-4, minimize=true, maxiters=100)

Run a standalone MMTK-like normal-mode calculation using AMBER gradients and a
finite-difference force-constant matrix.

This is intended for small validation systems such as AlaAla and bALA1. It is
O((3N)^2) in storage and requires `2 * 3N` force evaluations, so it is not the
right tool for whole CB1 heavy-atom calculations.
"""
function mmtk_like_nma(pdb_path::AbstractString; delta::Float64=1e-4,
                       minimize::Bool=true, maxiters::Int=100,
                       verbose::Bool=true)
    mol, ff, initial_energy, minimized_energy =
        load_prepared_system(pdb_path; minimize=minimize, maxiters=maxiters)

    x0 = coordinate_vector(mol)
    coords, masses = coordinates_and_masses(mol)

    H = finite_difference_force_constant_matrix!(ff, x0; delta=delta, verbose=verbose)
    Hmw, inv_sqrt_m = mass_weight_force_constants(H, masses)

    eig = eigen(Hmw)
    cartesian_modes = Diagonal(inv_sqrt_m) * eig.vectors
    frequencies = forcefield_frequencies_thz(eig.values)

    return MMTKLikeModes(
        coords,
        masses,
        Matrix(H),
        Matrix(Hmw),
        eig.values,
        eig.vectors,
        cartesian_modes,
        frequencies,
        initial_energy,
        minimized_energy,
    )
end

function write_frequency_table(path::AbstractString, result::MMTKLikeModes)
    open(path, "w") do io
        println(io, "# mode\teigenvalue\tfrequency_thz")
        for i in eachindex(result.eigenvalues)
            f = result.frequencies_thz[i]
            ftxt = iszero(imag(f)) ? string(real(f)) : string(imag(f), "im")
            println(io, i, '\t', result.eigenvalues[i], '\t', ftxt)
        end
    end
    return path
end

function write_mode_table(path::AbstractString, result::MMTKLikeModes)
    open(path, "w") do io
        for i in eachindex(result.frequencies_thz)
            f = result.frequencies_thz[i]
            ftxt = iszero(imag(f)) ? string(real(f)) : string(imag(f), "im")
            print(io, i - 1, ' ', ftxt)
            for x in result.cartesian_modes[:, i]
                print(io, ' ', x)
            end
            println(io)
        end
    end
    return path
end

function parse_cli(args)
    pdb_path = isempty(args) ? joinpath(pkgdir(BiochemicalAlgorithms), "test", "data", "AlaAla.pdb") : args[1]
    outdir = joinpath(pkgdir(BiochemicalAlgorithms), "tutorials", "nma", "outputs", "mmtk_like")
    delta = 1e-4
    maxiters = 100
    minimize = true

    for arg in args[2:end]
        if startswith(arg, "--outdir=")
            outdir = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--delta=")
            delta = parse(Float64, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--maxiters=")
            maxiters = parse(Int, split(arg, "=", limit=2)[2])
        elseif arg == "--no-minimize"
            minimize = false
        else
            error("Unknown argument: $arg")
        end
    end

    return (; pdb_path, outdir, delta, maxiters, minimize)
end

function main(args=ARGS)
    opts = parse_cli(args)
    mkpath(opts.outdir)

    println("Running MMTK-like force-field NMA")
    println("Input:     ", abspath(opts.pdb_path))
    println("Outdir:    ", abspath(opts.outdir))
    println("delta:     ", opts.delta)
    println("minimize:  ", opts.minimize)
    println("maxiters:  ", opts.maxiters)

    result = mmtk_like_nma(
        opts.pdb_path;
        delta=opts.delta,
        minimize=opts.minimize,
        maxiters=opts.maxiters,
        verbose=true,
    )

    stem = splitext(basename(opts.pdb_path))[1]
    freq_file = joinpath(opts.outdir, "$(stem)_mmtk_like_frequencies.txt")
    mode_file = joinpath(opts.outdir, "$(stem)_mmtk_like_modes_clean.txt")

    write_frequency_table(freq_file, result)
    write_mode_table(mode_file, result)

    println("Initial energy:   ", result.initial_energy, " kJ/mol")
    println("Minimized energy: ", result.minimized_energy, " kJ/mol")
    println("Modes:            ", length(result.eigenvalues))
    println("Frequencies:      ", freq_file)
    println("Modes table:      ", mode_file)

    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end

