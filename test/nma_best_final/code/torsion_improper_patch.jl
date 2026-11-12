# torsion_improper_patch.jl  -- patch for src/forcefields/common/torsion_component.jl
#
# PROBLEM (seen in output): for a 3-bonded centre the code builds ONE atom
# arrangement, with a4 = partner of the first bond. Amber impropers are
# "X-X-C-O" / "X-X-N-H": the 4th atom must be O (or H). If O/H is not first in
# the bond list the improper is silently missed (your run: a4 was CH3 for ACE C).
#
# FIX: try each of the 3 neighbours as a4 (and both orders of a1,a2), keep the
# first arrangement that has parameters, report only if none matches.
#
# HOW TO APPLY (2 edits in torsion_component.jl):
#  1. Replace the whole function _try_assign_torsion! with version (1) below.
#  2. In _setup_improper_torsions!, replace the block
#         bs = bonds(atom)
#         if length(bs) == 3
#             if get_full_name(atom) ∈ impropers.name
#                 ... (the three nested i_1/i_2/i_3 loops) ...
#             end
#         end
#     with block (2) below.
# The call in _setup_proper_torsions! stays unchanged (report defaults to true).
# NOT RUN by me: after patching, validate_ff.jl must show "improper: 4".

# ---------------- (1) replacement for _try_assign_torsion! ---------------
function _try_assign_torsion!(
    tc::TorsionComponent{T},
    torsions::Vector{CosineTorsion{T}},
    torsion_combinations::Dict{K, V},
    a1::Atom{T},
    a2::Atom{T},
    a3::Atom{T},
    a4::Atom{T},
    V_factor::T,
    ϕ₀_factor::T,
    is_proper::Bool;
    report::Bool=true
)::Bool where {T, K, V}

    ff = tc.ff

    type_a1::String = a1.atom_type
    type_a2::String = a2.atom_type
    type_a3::String = a3.atom_type
    type_a4::String = a4.atom_type

    pt = @coalesce(
        get(torsion_combinations, (type_a1, type_a2, type_a3, type_a4,), missing),
        get(torsion_combinations, (type_a4, type_a3, type_a2, type_a1,), missing),
        get(torsion_combinations, ("*",     type_a2, type_a3, type_a4,), missing),
        get(torsion_combinations, ("*",     type_a3, type_a2, type_a1,), missing),
        get(torsion_combinations, (type_a1, type_a2, type_a3, "*",    ), missing),
        get(torsion_combinations, (type_a4, type_a3, type_a2, "*",    ), missing),
        get(torsion_combinations, ("*",     type_a2, type_a3, "*"    ,), missing),
        get(torsion_combinations, ("*",     type_a3, type_a2, "*"    ,), missing),
        get(torsion_combinations, ("*",     "*",     type_a3, type_a4,), missing)
    )

    if ismissing(pt)
        if report
            push!(tc.unassigned_torsions, (a1, a2, a3, a4, is_proper))
            push!(ff.unassigned_atoms, a1)
            push!(ff.unassigned_atoms, a2)
            push!(ff.unassigned_atoms, a3)
            push!(ff.unassigned_atoms, a4)
            if length(ff.unassigned_atoms) > ff.options[:max_number_of_unassigned_atoms]
                throw(TooManyErrors())
            end
        end
        return false
    end

    n = only(filter(d -> d.n == "N", pt)).div
    pts = filter(d -> d.n != "N", pt)

    if n != nrow(pts)
        @warn "Invalid torsion configuration for atom types $(type_a1)-$(type_a2)-$(type_a3)-$(type_a4)!"
        return false
    end

    push!(torsions,
        CosineTorsion(
            V_factor  .* getproperty.(pts, :V),
            ϕ₀_factor .* getproperty.(pts, :phi0),
            getproperty.(pts, :f),
            getproperty.(pts, :div),
            a1, a2, a3, a4
        )
    )
    return true
end

# ---------------- (2) replacement block inside _setup_improper_torsions! ---
#    for atom in atoms(ff.system)
#        bs = bonds(atom)
#        # the ini may list centres as "ALA:N" while the atom's full name is "ALA-M:N",
#        # so accept the name with OR without the variant suffix:
#        if length(bs) == 3 &&
#           (get_full_name(atom) ∈ impropers.name ||
#            get_full_name(atom, FullNameType.NO_VARIANT_EXTENSIONS) ∈ impropers.name)
#            nbrs = [get_partner(b, atom) for b in bs]
#            assigned = false
#            for i4 in 1:3
#                a4 = nbrs[i4]
#                o  = [nbrs[k] for k in 1:3 if k != i4]
#                for (a1, a2) in ((o[1], o[2]), (o[2], o[1]))
#                    if _try_assign_torsion!(tc, improper_torsions, torsion_dict,
#                                            a1, a2, atom, a4, V_factor, ϕ₀_factor, false;
#                                            report=false)
#                        assigned = true
#                        break
#                    end
#                end
#                assigned && break
#            end
#            # nothing matched: report once
#            assigned || _try_assign_torsion!(tc, improper_torsions, torsion_dict,
#                                             nbrs[3], nbrs[2], atom, nbrs[1],
#                                             V_factor, ϕ₀_factor, false)
#        end
#    end


#FROM  old working modified code; REPLACE with the above block (2) in torsion_component.jl
for atom in atoms(ff.system)
        bs = bonds(atom)
        # the ini may list centres as "ALA:N" while the atom's full name is "ALA-M:N",
        # so accept the name with OR without the variant suffix:
        if length(bs) == 3 &&
           (get_full_name(atom) ∈ impropers.name ||
            get_full_name(atom, FullNameType.NO_VARIANT_EXTENSIONS) ∈ impropers.name)
            nbrs = [get_partner(b, atom) for b in bs]
            assigned = false
            for i4 in 1:3
                a4 = nbrs[i4]
                o  = [nbrs[k] for k in 1:3 if k != i4]
                for (a1, a2) in ((o[1], o[2]), (o[2], o[1]))
                    if _try_assign_torsion!(tc, improper_torsions, torsion_dict,
                                            a1, a2, atom, a4, V_factor, ϕ₀_factor, false;
                                            report=false)
                        assigned = true
                        break
                    end
                end
                assigned && break
            end
            # nothing matched: report once
            assigned || _try_assign_torsion!(tc, improper_torsions, torsion_dict,
                                             nbrs[3], nbrs[2], atom, nbrs[1],
                                             V_factor, ϕ₀_factor, false)
        end
    end
