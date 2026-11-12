############################################################
# validate_forcefield.jl
# All-atom force-field NMA (AmberFF Hessian) versus MMTK (Hinsen), blocked
# alanine and Ala-Ala, with BOTH improper conventions:
#   :amber  Julia/BALL default (Amber parameter-file convention)
#   :mmtk   MMTK's convention: centre first, others sorted by (type, index),
#           force constant halved   [MMForceField._order, AmberImproperParameters]
#
# Run from the project root; needs the BiochemicalAlgorithms project.
#   julia --project=. test\nma_best_final\code\validate_forcefield.jl [bala1 alaala ...]
#
# Inputs are in ../mmtk (all from the same mmtk_compare.py run):
#   <tag>_min.pdb  mmtk_<tag>_freqs.txt  mmtk_<tag>_gradients.txt  mmtk_<tag>_energies.txt
#
# Also needs the improper patch applied to torsion_component.jl
# (torsion_improper_patch.jl), otherwise impropers are mis-assigned.
# Exit code 1 if any :mmtk row fails its criteria.
############################################################
using BiochemicalAlgorithms, StaticArrays, LinearAlgebra, Printf, Statistics
include(joinpath(@__DIR__, "nma_core.jl"))
include(joinpath(@__DIR__, "nma_report.jl"))

const DATA_ROOT = normpath(joinpath(@__DIR__, ".."))
const MMTK_DIR = joinpath(DATA_ROOT, "mmtk")
const RESULTS_DIR = joinpath(DATA_ROOT, "validation_forcefield")
mkpath(RESULTS_DIR)

const TAGS      = isempty(ARGS) ? ["bala1", "alaala"] : String.(ARGS)
const AMBER94   = ball_data_path("forcefields/AMBER/amber94.ini")
const AMBER96   = ball_data_path("forcefields/AMBER/amber96.ini")
const MD        = joinpath(RESULTS_DIR, "validation_forcefield.md")
isfile(MD) && rm(MD)

# criteria for the :mmtk rows (kJ/mol, kJ/mol/Å, relative)
const TOL_ETOT  = 1e-4
const TOL_GRAD  = 1e-5
const TOL_EMIN  = 1e-4
const TOL_FMEAN = 1e-4
const TOL_FMAX  = 1e-3

# ============================================================
# parameter-file patch: NME charges/types under the key "NME:"
# (the ini keys them as "NME-C:" but the NME fragment has no variant)
# ============================================================
function ini_with_nme(src)
    nme = [replace(l, "NME-C:" => "NME:") for l in eachline(AMBER96)
           if occursin(r"^\s*1\.0\s+NME-C:", l)]
    length(nme) == 6 || error("expected 6 NME-C lines in amber96.ini, found $(length(nme))")
    out = String[]; done = false
    for l in eachline(src)
        push!(out, l)
        if !done && occursin(r"^\s*1\.0\s+ACE:O\s", l)
            append!(out, nme); done = true
        end
    end
    done || error("could not find the 'ACE:O' charge line in $src")
    dst = joinpath(mktempdir(), basename(src))
    write(dst, join(out, "\n") * "\n")
    return dst
end

# build_bonds! can leave a hydrogen without a bond (template names differ from the PDB)
function fix_orphan_hydrogens!(mol)
    digits_ = ['0','1','2','3','4','5','6','7','8','9']
    is_h(a) = startswith(lstrip(strip(a.name), digits_), "H")
    atms0 = collect(atoms(mol))
    bonded = Set{Int}()
    for b in bonds(mol); push!(bonded, b.a1); push!(bonded, b.a2); end
    for a in atms0
        (a.idx in bonded || !is_h(a)) && continue
        best, bd = nothing, Inf
        for c in atms0
            (c.idx == a.idx || is_h(c)) && continue
            d = norm(c.r - a.r)
            d < bd && ((best, bd) = (c, d))
        end
        if best !== nothing && bd < 1.3
            println("    adding missing bond ", get_full_name(best), " - ", get_full_name(a),
                    "  d = ", round(bd, digits=3))
            Bond(mol, best.idx, a.idx, BondOrder.Single)
        end
    end
end

# ============================================================
# system construction
# ============================================================
function build_system(pdb_path)
    mol = load_pdb(pdb_path, Float64)
    fdb = FragmentDB{Float64}()
    normalize_names!(mol, fdb)

    # MMTK writes N-terminal hydrogens as HT1..3 (normalized 1HT..3HT); the ini uses 1H..3H.
    # Rename BEFORE reconstruction so MMTK's coordinates are kept and nothing is added.
    for a in collect(atoms(mol))
        nm = String(strip(a.name))
        if nm in ("1HT", "2HT", "3HT")
            a.name = replace(nm, "HT" => "H")
        end
    end
    reconstruct_fragments!(mol, fdb)

    # capped molecules (ACE/NME): reconstruct_fragments! adds free-terminus atoms; remove them
    has_caps = any(strip(parent_fragment(a).name) in ("ACE", "NME") for a in atoms(mol))
    if has_caps
        for a in collect(atoms(mol))
            strip(a.name) in ("1H", "2H", "3H", "OXT") && delete!(a)
        end
    end
    build_bonds!(mol, fdb)
    fix_orphan_hydrogens!(mol)

    src = isfile(AMBER94) ? AMBER94 : AMBER96
    ff  = AmberFF(mol, ini_with_nme(src))
    tc  = only(filter(c -> c isa BiochemicalAlgorithms.TorsionComponent, ff.components))
    return mol, ff, tc, has_caps
end

# MMTK improper: central atom first, the other three sorted by (atom type, atom index);
# force constant halved.
function mmtk_order(t; scale=0.5)
    cand = [t.a1, t.a2, t.a4]
    o = cand[sortperm([(string(a.atom_type), a.idx) for a in cand])]
    return CosineTorsion(scale .* t.V, t.ϕ₀, t.f, t.div, t.a3, o[1], o[2], o[3])
end

function element_symbol(a)                       # a.element is Unknown for HH31-type names
    e = string(a.element)
    haskey(ATOMIC_MASSES, e) && return e
    nm = lstrip(strip(a.name), ['0','1','2','3','4','5','6','7','8','9'])
    return string(first(nm))
end

# ============================================================
# reference data from MMTK
# ============================================================
function read_mmtk_energies(path)
    d = Dict{String,Float64}()
    for l in eachline(path)
        p = split(strip(l)); isempty(p) && continue
        if p[1] == "minimum_energy"
            d["minimum_energy"] = parse(Float64, p[2])
        elseif p[1] == "reread" && length(p) >= 3
            d[p[2]] = parse(Float64, p[3])
        end
    end
    return d
end

function read_mmtk_gradients(path)
    mm = NamedTuple[]
    for l in eachline(path)
        p = split(strip(l)); length(p) == 7 || continue
        push!(mm, (name=p[1], r=parse.(Float64, p[2:4]), g=parse.(Float64, p[5:7])))
    end
    return mm
end

read_vector(path) = [parse(Float64, strip(l)) for l in eachline(path) if !isempty(strip(l))]

# ============================================================
# measurements
# ============================================================
function energy_terms(ff)
    E = compute_energy!(ff)
    d = ff.energy
    return (bond    = d["Bond Stretches"],
            angle   = d["Angle Bends"],
            vdw     = d["Van der Waals"] + get(d, "Hydrogen Bonds", 0.0),
            elec    = d["Electrostatic"],
            torsion = d["Proper Torsion"] + d["Improper Torsion"],
            total   = E)
end

function julia_gradient(ff)
    atms = collect(atoms(ff.system))
    update!(ff); compute_forces!(ff)
    return reduce(vcat, [-collect(a.F) for a in atms])
end

"RMS over atoms of |g_MMTK - g_Julia| (kJ/mol/Å); atoms matched by position."
function gradient_rms(ff, mm)
    atms = collect(atoms(ff.system))
    gJ = julia_gradient(ff)
    gM = zeros(3 * length(atms))
    for (i, a) in enumerate(atms)
        d = [norm(collect(a.r) - m.r) for m in mm]
        k = argmin(d)
        d[k] < 0.02 || error("no MMTK atom within 0.02 Å of $(get_full_name(a))")
        gM[3i-2:3i] .= mm[k].g
    end
    if norm(gM + gJ) < norm(gM - gJ)             # sign convention
        gM .*= -1
    end
    return sqrt(sum(abs2, gM .- gJ) / length(atms))
end

"Max |analytic - finite difference| / force scale per component, at a jittered geometry."
function force_ratios(ff; h=1e-5, ntest=20, jitter=0.08)
    atms = collect(atoms(ff.system))
    r_orig = [a.r for a in atms]
    for a in atms
        a.r = a.r + SVector{3,Float64}(jitter .* (rand(3) .- 0.5))
    end
    update!(ff)
    out = Dict{String,Float64}()
    for comp in ff.components
        maxerr, maxF = 0.0, 0.0
        for _ in 1:ntest
            i = rand(eachindex(atms)); d = rand(1:3)
            r0 = atms[i].r
            shifted(δ) = (v = collect(r0); v[d] += δ; SVector{3,Float64}(v))
            atms[i].r = shifted(+h); update!(ff); Ep = compute_energy!(comp)
            atms[i].r = shifted(-h); update!(ff); Em = compute_energy!(comp)
            atms[i].r = r0;          update!(ff)
            for a in atms; a.F = zero(a.F); end
            compute_forces!(comp)
            fd = -(Ep - Em) / (2h)
            maxerr = max(maxerr, abs(atms[i].F[d] - fd))
            maxF   = max(maxF, abs(fd))
        end
        out[comp.name] = maxerr / max(maxF, eps())
    end
    for (a, r) in zip(atms, r_orig); a.r = r; end
    update!(ff)
    return out
end

"Hessian by central differences of the analytic gradient (kJ/mol/Å²)."
function amber_hessian(ff; h=1e-4)
    atms = collect(atoms(ff.system))
    N = length(atms)
    x0 = reduce(vcat, [collect(a.r) for a in atms])
    function grad(x)
        for (i, a) in enumerate(atms)
            a.r = SVector{3,Float64}(x[3i-2], x[3i-1], x[3i])
        end
        update!(ff); compute_forces!(ff)
        g = zeros(3N)
        for (i, a) in enumerate(atms); g[3i-2:3i] .= -a.F; end
        return g
    end
    H = zeros(3N, 3N)
    for k in 1:3N
        xp = copy(x0); xp[k] += h
        xm = copy(x0); xm[k] -= h
        H[:, k] = (grad(xp) .- grad(xm)) ./ (2h)
    end
    grad(x0); update!(ff)
    asym = norm(H - H') / norm(H)
    return (H .+ H') ./ 2, asym
end

"Newton steps (zero modes skipped, bounded step) to drive |g| below gtol."
function newton_polish!(ff; iters=10, h=1e-4, maxstep=0.05, gtol=1e-5)
    atms = collect(atoms(ff.system))
    for _ in 1:iters
        g = julia_gradient(ff)
        norm(g) < gtol && return norm(g)
        H, _ = amber_hessian(ff; h=h)
        F = eigen(Symmetric(H)); λ, V = F.values, F.vectors
        tol = 1e-7 * maximum(abs, λ)
        keep = findall(>(tol), λ)
        dx = -V[:, keep] * ((V[:, keep]' * g) ./ λ[keep])
        s = min(1.0, maxstep / maximum(abs, dx))
        for (i, a) in enumerate(atms)
            a.r = a.r + SVector{3,Float64}(s * dx[3i-2], s * dx[3i-1], s * dx[3i])
        end
        update!(ff)
    end
    return norm(julia_gradient(ff))
end

# ============================================================
# one case = (system tag, improper convention)
# ============================================================
function run_case(tag::String, conv::Symbol)
    println("\n", "="^70, "\n  ", tag, "  |  improper convention: ", conv, "\n", "="^70)
    files = Dict(:pdb   => joinpath(MMTK_DIR, "$(tag)_min.pdb"),
                 :freq  => joinpath(MMTK_DIR, "mmtk_$(tag)_freqs.txt"),
                 :grad  => joinpath(MMTK_DIR, "mmtk_$(tag)_gradients.txt"),
                 :ener  => joinpath(MMTK_DIR, "mmtk_$(tag)_energies.txt"))
    for (k, f) in files
        isfile(f) || error("missing $(basename(f)) (run mmtk_compare.py $tag $tag and copy the files)")
    end

    mol, ff, tc, has_caps = build_system(files[:pdb])
    conv == :mmtk && (tc.improper_torsions = [mmtk_order(t) for t in tc.improper_torsions])
    nat = natoms(mol)
    println("  atoms: $nat   bonds: $(nbonds(mol))   capped: $has_caps   ",
            "impropers: $(length(tc.improper_torsions))   warnings: $(count_warnings(ff))")
    count_warnings(ff) == 0 || error("force field has setup warnings; fix before validating")

    # --- energies at MMTK's (re-read) geometry ---------------------------
    E   = energy_terms(ff)
    M   = read_mmtk_energies(files[:ener])
    ref = (bond=M["harmonic_bond"], angle=M["harmonic_bond_angle"], vdw=M["Lennard-Jones"],
           elec=M["electrostatic"], torsion=M["cosine_dihedral_angle"], total=M["TOTAL"])
    dE  = (bond=E.bond-ref.bond, angle=E.angle-ref.angle, vdw=E.vdw-ref.vdw,
           elec=E.elec-ref.elec, torsion=E.torsion-ref.torsion, total=E.total-ref.total)

    # --- gradient, force/energy consistency ------------------------------
    grms   = gradient_rms(ff, read_mmtk_gradients(files[:grad]))
    fr     = force_ratios(ff)
    fratio = maximum(values(fr))

    # --- minimise, polish ------------------------------------------------
    optimize_structure!(ff; maxiters=5000)
    gfin = newton_polish!(ff)
    Emin = compute_energy!(ff)

    # --- Hessian, modes ---------------------------------------------------
    H, asym = amber_hessian(ff; h=1e-4)
    atms = collect(atoms(mol))
    masses = [get(ATOMIC_MASSES, element_symbol(a), NaN) for a in atms]
    any(isnan, masses) && error("unknown element mass")
    vals, _, _ = compute_modes(H, masses)
    nz = count_zero_modes(vals; tol=1e-6)
    f  = frequencies(vals; unit=:cm1)
    fj = sort(f[nz+1:end])
    fm = sort(read_vector(files[:freq]))
    n  = min(length(fj), length(fm))
    length(fj) == length(fm) || @warn "mode count differs: Julia $(length(fj)), MMTK $(length(fm))"
    rel = abs.(fj[1:n] .- fm[1:n]) ./ fm[1:n]
    lowm  = fm[1:n] .< 500
    highm = fm[1:n] .> 2500
    nimag = count(<(0), f[nz+1:end])

    println(@sprintf("  ΔE_total = %+.2e   grad RMS = %.2e   Emin: Julia %.6f  MMTK %.6f",
                     dE.total, grms, Emin, M["minimum_energy"]))
    println(@sprintf("  frequencies: mean rel. diff %.2e, max %.2e, lowest mode %.4f vs %.4f cm⁻¹",
                     mean(rel), maximum(rel), fj[1], fm[1]))

    return (tag=tag, conv=conv, natoms=nat, E=E, ref=ref, dE=dE, grms=grms, fratio=fratio,
            fr=fr, gfin=gfin, asym=asym, Emin=Emin, Emin_mmtk=M["minimum_energy"],
            nmodes=n, nzero=nz, nimag=nimag, mean_rel=mean(rel), max_rel=maximum(rel),
            low_rel=any(lowm) ? mean(rel[lowm]) : NaN,
            high_rel=any(highm) ? mean(rel[highm]) : NaN,
            ratio=median(fj[1:n] ./ fm[1:n]), f1_j=fj[1], f1_m=fm[1])
end

# ============================================================
# run all cases, report
# ============================================================
results = Any[]
failed_cases = String[]
for tag in TAGS, conv in (:amber, :mmtk)
    try
        push!(results, run_case(tag, conv))
    catch err
        println("  CASE FAILED: ", tag, " / ", conv, " -> ", sprint(showerror, err))
        push!(failed_cases, "$tag/$conv")
    end
end

isempty(results) && (println("no case completed"); exit(1))

println("\n\n", "#"^70, "\n  VALIDATION TABLES\n", "#"^70)

report_table("Table 1 — energy terms at the SAME geometry, Julia − MMTK (kJ/mol)",
    ["System", "Improper convention", "Atoms", "Δ bond", "Δ angle", "Δ vdW", "Δ electrostatic",
     "Δ torsion", "Δ total", "MMTK total"],
    [[r.tag, string(r.conv), r.natoms, r.dE.bond, r.dE.angle, r.dE.vdw, r.dE.elec,
      r.dE.torsion, r.dE.total, r.ref.total] for r in results];
    md=MD, tsv=joinpath(RESULTS_DIR, "validation_forcefield_energies.tsv"))

report_table("Table 2 — gradients, force/energy consistency and minimum",
    ["System", "Improper convention", "gradient RMS vs MMTK (kJ/mol/Å)",
     "max force/energy ratio", "gradient norm after polish", "Hessian asymmetry",
     "E_min Julia", "E_min MMTK", "Δ E_min"],
    [[r.tag, string(r.conv), r.grms, r.fratio, r.gfin, r.asym, r.Emin, r.Emin_mmtk,
      r.Emin - r.Emin_mmtk] for r in results];
    md=MD, tsv=joinpath(RESULTS_DIR, "validation_forcefield_gradients.tsv"))

report_table("Table 3 — normal-mode frequencies vs MMTK (cm⁻¹)",
    ["System", "Improper convention", "Modes", "zero modes", "imaginary", "mean rel. diff",
     "max rel. diff", "low modes (<500)", "high modes (>2500)", "median ratio",
     "lowest mode Julia", "lowest mode MMTK"],
    [[r.tag, string(r.conv), r.nmodes, r.nzero, r.nimag, r.mean_rel, r.max_rel, r.low_rel,
      r.high_rel, r.ratio, r.f1_j, r.f1_m] for r in results];
    md=MD, tsv=joinpath(RESULTS_DIR, "validation_forcefield_modes.tsv"))

# ---- pass/fail for the MMTK-convention rows ------------------------------
rows4 = Vector{Any}[]
nfail = Ref(0)
for r in results
    r.conv == :mmtk || continue
    checks = [("abs(Δ E_total) < $(TOL_ETOT)",     abs(r.dE.total) < TOL_ETOT,     abs(r.dE.total)),
              ("gradient RMS < $(TOL_GRAD)",      r.grms < TOL_GRAD,              r.grms),
              ("abs(Δ E_min) < $(TOL_EMIN)",       abs(r.Emin - r.Emin_mmtk) < TOL_EMIN, abs(r.Emin - r.Emin_mmtk)),
              ("mean rel. freq diff < $(TOL_FMEAN)", r.mean_rel < TOL_FMEAN,       r.mean_rel),
              ("max rel. freq diff < $(TOL_FMAX)", r.max_rel < TOL_FMAX,           r.max_rel),
              ("6 zero modes, 0 imaginary",       r.nzero == 6 && r.nimag == 0,   r.nzero)]
    for (name, ok, val) in checks
        push!(rows4, [r.tag, name, val, ok ? "PASS" : "FAIL"])
        ok || (nfail[] += 1)
    end
end
report_table("Table 4 — pass/fail with MMTK's improper convention imposed",
    ["System", "Criterion", "Value", "Status"], rows4;
    md=MD, tsv=joinpath(RESULTS_DIR, "validation_forcefield_checks.tsv"))

isempty(failed_cases) || println("\nCases that failed to run: ", join(failed_cases, ", "))
println("\nwritten: validation_forcefield.md and validation_forcefield_*.tsv")
all_ok = nfail[] == 0 && isempty(failed_cases)
println(all_ok ? "\nALL MMTK-CONVENTION CHECKS PASSED" :
        "\n$(nfail[]) CHECK(S) FAILED, $(length(failed_cases)) CASE(S) NOT RUN")
all_ok || exit(1)
