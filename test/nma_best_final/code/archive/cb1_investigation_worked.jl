############################################################
# cb1_investigation.jl — CB1 receptor: WT vs F238L with the elastic network model
#
#   julia test\nma\cb1_investigation.jl WT.pdb MUT.pdb [REP1.pdb REP2.pdb ...]
#
# (needs only Julia's standard library + nma_core.jl, no BiochemicalAlgorithms)
#
# WT, MUT : the two structures prepared with the IDENTICAL protocol
# REPn    : WT pushed through the same protocol again (different repack seed /
#           minimisation settings). They define the noise floor; without >= 3
#           replicates no result is declared significant.
#
# What it does
#   1  pre-flight: numbering (PHE in WT, LEU in the mutant), graph, zero modes
#   2  WT vs mutant on common atoms: dMSF, dDCCM, RMSIP, eigenvalue shifts
#   3  the same metrics for WT vs each replicate -> noise floor, z-scores
#   4  cutoff sensitivity of the WT-vs-mutant result
#   5  stiffness-perturbation scan (mutation-independent): how strongly does each
#      residue's stiffness couple to the region of interest, and where does the
#      target residue rank?
#   6  optional B-factor calibration of gamma from the ORIGINAL experimental PDB
#
# A Calpha network with identical Calpha coordinates CANNOT distinguish WT from
# mutant (the script says so). Use CB1_SEL=heavy for a side-chain-aware model.
#
# Settings (environment variables, defaults in brackets)
#   CB1_SEL [ca|heavy]   CB1_MODEL [dw|anm]   CB1_CHAIN [A]   CB1_TARGET [237]
#   CB1_RANGES [99-306,332-414]   CB1_REGION [332-368]   CB1_CUTOFF [12 ca | 7 heavy]
#   CB1_NMODES [20]   CB1_MASSES [none|residue]   CB1_PRS_FACTOR [0.5]
#   CB1_PRS_STRIDE [1]   CB1_BFAC_PDB [path to the original experimental PDB]
############################################################
using LinearAlgebra, Statistics, Printf, Random
include(joinpath(@__DIR__, "nma_core.jl"))
include(joinpath(@__DIR__, "nma_report.jl"))

env(k, d) = get(ENV, k, d)
function parse_ranges(s)
    out = Tuple{Int,Int}[]
    for r in split(s, ",")
        isempty(strip(r)) && continue
        a, b = split(strip(r), "-")
        push!(out, (parse(Int, a), parse(Int, b)))
    end
    return out
end

const SEL        = Symbol(lowercase(env("CB1_SEL", "ca")))
const MODEL      = Symbol(lowercase(env("CB1_MODEL", "dw")))
const CHAIN      = only(env("CB1_CHAIN", "A"))
const TARGET     = parse(Int, env("CB1_TARGET", "237"))        # PHE in PDB numbering ("F238")
const RANGES     = parse_ranges(env("CB1_RANGES", "99-306,332-414"))
const REGION     = first(parse_ranges(env("CB1_REGION", "332-368")))
const CUTOFF     = parse(Float64, env("CB1_CUTOFF", SEL == :ca ? "12.0" : "7.0"))
const NMODES     = parse(Int, env("CB1_NMODES", "20"))
const MASS_MODE  = Symbol(lowercase(env("CB1_MASSES", "none")))
const PRS_FACTOR = parse(Float64, env("CB1_PRS_FACTOR", "0.5"))
const PRS_STRIDE = parse(Int, env("CB1_PRS_STRIDE", "1"))
const BFAC_PDB   = env("CB1_BFAC_PDB", "")
const KT         = 2.494
const DIGITS     = ['0','1','2','3','4','5','6','7','8','9']
const KFUN       = MODEL == :dw ? kfun_distance_weighted(1.0) : kfun_constant(1.0)

SEL in (:ca, :heavy)   || error("CB1_SEL must be ca or heavy")
MODEL in (:dw, :anm)   || error("CB1_MODEL must be dw or anm")
length(ARGS) >= 2 || error("usage: julia cb1_investigation.jl WT.pdb MUT.pdb [replicate1.pdb ...]")

const OUT = joinpath(@__DIR__, "cb1_results")
mkpath(OUT)
const MD = joinpath(OUT, "cb1_summary.md")
isfile(MD) && rm(MD)

wt_path  = abspath(ARGS[1])
mut_path = abspath(ARGS[2])
rep_paths = abspath.(ARGS[3:end])
for p in vcat(wt_path, mut_path, rep_paths)
    isfile(p) || error("file not found: $p")
end

# ============================================================
# reading
# ============================================================
field(l, a, b) = length(l) < a ? "" : l[a:min(b, lastindex(l))]

function read_structure(path)
    keys_ = Tuple{Char,Int,String}[]; xyz = Vector{Vector{Float64}}(); rn = String[]
    for l in eachline(path)
        startswith(l, "ATOM") || continue
        alt = field(l, 17, 17)
        (alt == " " || alt == "" || alt == "A") || continue
        name = String(strip(field(l, 13, 16)))
        el = String(strip(field(l, 77, 78)))
        isempty(el) && (el = string(first(lstrip(name, DIGITS))))
        el == "H" && continue
        occursin(r"^\d*H", name) && continue
        ch = only(rpad(field(l, 22, 22), 1))
        resi = tryparse(Int, strip(field(l, 23, 26)))
        resi === nothing && continue
        ch == CHAIN || continue
        any(a <= resi <= b for (a, b) in RANGES) || continue
        SEL == :ca && name != "CA" && continue
        push!(keys_, (ch, resi, name)); push!(rn, String(strip(field(l, 18, 20))))
        push!(xyz, [parse(Float64, strip(field(l, a, b))) for (a, b) in ((31,38), (39,46), (47,54))])
    end
    isempty(keys_) && error("no atoms selected from $path (check CB1_CHAIN / CB1_RANGES / CB1_SEL)")
    return (keys=keys_, coords=permutedims(hcat(xyz...)), resn=rn, resi=[k[2] for k in keys_])
end

function masses_for(S)
    MASS_MODE == :none && return ones(length(S.keys))
    if SEL == :ca
        return [get(RESIDUE_MASSES, uppercase(r[1:min(3, end)]), 110.0) for r in S.resn]
    end
    return [get(ATOMIC_MASSES, string(first(lstrip(k[3], DIGITS))), 12.011) for k in S.keys]
end

dofs(idx) = reduce(vcat, [[3i-2, 3i-1, 3i] for i in idx])

# ============================================================
# model
# ============================================================
"""
Mass-weighted ENM solution. Returns eigenvalues, Cartesian mode matrix U, per-atom
MSF (units of kT/γ) and, optionally, the full Cartesian covariance (kT = 1).
"""
function solve(coords, masses, cutoff; pairweight=nothing, want_cov=false)
    H = build_hessian(coords, KFUN; cutoff=cutoff, pairweight=pairweight)
    s = 1.0 ./ sqrt.(repeat(masses, inner=3))
    F = eigen(Symmetric(Diagonal(s) * H * Diagonal(s)))
    nz = count_zero_modes(F.values; tol=1e-6)
    λ = F.values[nz+1:end]
    V = F.vectors[:, nz+1:end]
    U = s .* V
    W = U ./ sqrt.(λ)'
    N = length(masses)
    msf = [sum(abs2, W[3i-2:3i, :]) for i in 1:N]
    C = want_cov ? W * W' : nothing
    return (vals=F.values, nz=nz, U=U, msf=msf, C=C)
end

"Compare two solved models on their common atoms."
function compare(A, B, SA, SB)
    ixb = Dict(k => i for (i, k) in enumerate(B.keys))
    ixa = Dict(k => i for (i, k) in enumerate(A.keys))
    common = [k for k in A.keys if haskey(ixb, k)]
    isempty(common) && error("no common atoms between the two structures")
    ia = [ixa[k] for k in common]; ib = [ixb[k] for k in common]
    res = sort(unique(k[2] for k in common))
    pos = Dict(r => Int[] for r in res)
    for (n, k) in enumerate(common); push!(pos[k[2]], n); end
    rA  = [mean(SA.msf[ia[pos[r]]]) for r in res]
    rB  = [mean(SB.msf[ib[pos[r]]]) for r in res]
    rel = (rB .- rA) ./ rA

    K  = min(NMODES, size(SA.U, 2), size(SB.U, 2))
    QA = Matrix(qr(SA.U[dofs(ia), 1:K]).Q)[:, 1:K]
    QB = Matrix(qr(SB.U[dofs(ib), 1:K]).Q)[:, 1:K]
    rmsp = rmsip(QA, QB, K; orth_warn=false)

    ireg = findall(r -> REGION[1] <= r <= REGION[2], res)
    it   = findfirst(==(TARGET), res)

    ddc, ddt = NaN, NaN
    if SA.C !== nothing && SB.C !== nothing && res == [k[2] for k in common]
        Da = dccm_from_cov(SA.C[dofs(ia), dofs(ia)])
        Db = dccm_from_cov(SB.C[dofs(ib), dofs(ib)])
        ddc = norm(Db .- Da) / length(res)
        it === nothing || (ddt = mean(Db[it, ireg] .- Da[it, ireg]))
    end
    eig = NaN
    if length(SA.vals) == length(SB.vals)
        λa = SA.vals[SA.nz+1:SA.nz+K]; λb = SB.vals[SB.nz+1:SB.nz+K]
        eig = mean(abs.(log.(λb ./ λa)))
    end
    dmax = maximum(norm(A.coords[ia[n], :] .- B.coords[ib[n], :]) for n in eachindex(common))
    return (res=res, rA=rA, rB=rB, rel=rel, rmsip=rmsp, omr=1 - rmsp,
            reg=isempty(ireg) ? NaN : mean(rel[ireg]),
            tgt=it === nothing ? NaN : rel[it],
            maxrel=maximum(abs, rel), ddc=ddc, ddt=ddt, eig=eig,
            ncommon=length(common), identical=dmax < 1e-6, dmax=dmax)
end

# ============================================================
# 1. load and pre-flight
# ============================================================
println("\n=== CB1 investigation ===")
println("  model       : ", MODEL == :dw ? "distance-weighted (k = γ/r²)" : "constant spring (k = γ)",
        ", selection ", SEL, ", cutoff ", CUTOFF, " Å")
println("  target      : residue ", TARGET, "   region of interest: ", REGION[1], "–", REGION[2],
        "   masses: ", MASS_MODE)

WT  = read_structure(wt_path)
MUT = read_structure(mut_path)
REP = [read_structure(p) for p in rep_paths]
M_wt, M_mut = masses_for(WT), masses_for(MUT)

if SEL == :heavy && 3 * length(WT.keys) > 6000
    @warn "heavy-atom model with $(length(WT.keys)) atoms needs several GB of memory; " *
          "consider CB1_SEL=ca or a smaller region"
end

resname_of(S, r) = (i = findfirst(==(r), S.resi); i === nothing ? "—" : S.resn[i])
println("  residues around the target in WT:")
println("    ", join(["$(r):$(resname_of(WT, r))" for r in TARGET-4:TARGET+4], "  "))

print("  solving WT and mutant ... "); flush(stdout)
S_wt  = solve(WT.coords,  M_wt,  CUTOFF; want_cov = SEL == :ca)
S_mut = solve(MUT.coords, M_mut, CUTOFF; want_cov = SEL == :ca)
println("done")

pre = Vector{Any}[
    ["atoms in the model", length(WT.keys), length(MUT.keys)],
    ["residues", length(unique(WT.resi)), length(unique(MUT.resi))],
    ["residue at target $(TARGET)", resname_of(WT, TARGET), resname_of(MUT, TARGET)],
    ["expected at target", "PHE", "LEU"],
    ["zero modes (expect 6)", S_wt.nz, S_mut.nz],
    ["softest non-zero eigenvalue", S_wt.vals[S_wt.nz+1], S_mut.vals[S_mut.nz+1]],
]
report_table("Pre-flight checks", ["item", "WT", "mutant"], pre; md=MD)

preflight_ok = true
if resname_of(WT, TARGET) != "PHE"
    println("  [WARN] WT residue at $TARGET is $(resname_of(WT, TARGET)), not PHE: check the numbering.")
    preflight_ok = false
end
if resname_of(MUT, TARGET) != "LEU"
    println("  [WARN] mutant residue at $TARGET is $(resname_of(MUT, TARGET)), not LEU: wrong file or numbering.")
    preflight_ok = false
end
if S_wt.nz != 6 || S_mut.nz != 6
    println("  [WARN] zero modes ≠ 6: the network is disconnected or strained; results are not valid.")
    preflight_ok = false
end

# ============================================================
# 2. WT vs mutant
# ============================================================
wm = compare(WT, MUT, S_wt, S_mut)
println("\n=== WT vs mutant on $(wm.ncommon) common atoms ===")
if wm.identical
    println("  [NOTE] common-atom coordinates are IDENTICAL (max difference ",
            @sprintf("%.1e", wm.dmax), " Å). With a ",
            SEL == :ca ? "Cα" : "heavy-atom", " model every WT-vs-mutant difference is zero except " *
            "mass effects: the model cannot see this mutation. Use a structure where the mutation " *
            "changed the coordinates, CB1_SEL=heavy, or the stiffness-perturbation scan below.")
end

# ============================================================
# 3. noise floor
# ============================================================
reps = NamedTuple[]
for (p, R) in zip(rep_paths, REP)
    SR = solve(R.coords, masses_for(R), CUTOFF; want_cov = SEL == :ca)
    SR.nz == 6 || println("  [WARN] replicate $(basename(p)): zero modes = $(SR.nz)")
    push!(reps, compare(WT, R, S_wt, SR))
end

specs = [("relative ΔMSF, region mean",         :reg,    :signed),
         ("relative ΔMSF, target residue",      :tgt,    :signed),
         ("max abs(relative ΔMSF)",                :maxrel, :mag),
         ("1 − RMSIP (low-mode subspace)",      :omr,    :mag),
         ("mean abs(log λ ratio) (spectrum shift)", :eig,   :mag),
         ("ΔDCCM, Frobenius / N",               :ddc,    :mag),
         ("ΔDCCM, target × region mean",        :ddt,    :signed)]

rows = Vector{Any}[]
verdicts = String[]
for (label, f, kind) in specs
    v = getfield(wm, f)
    vals = filter(!isnan, [getfield(r, f) for r in reps])
    if isnan(v)
        push!(rows, [label, v, NaN, NaN, NaN, "n/a"])
    elseif length(vals) < 3
        push!(rows, [label, v, isempty(vals) ? NaN : mean(vals), length(vals) < 2 ? NaN : std(vals), NaN,
                     "need ≥ 3 replicates"])
    else
        mu, sd = mean(vals), std(vals)
        z = sd > 0 ? (v - mu) / sd : (v == mu ? 0.0 : Inf)
        out = kind == :signed ? abs(z) > 3 : z > 3
        verdict = out ? "outside noise floor" : "within noise floor"
        push!(rows, [label, v, mu, sd, z, verdict]); push!(verdicts, verdict)
    end
end
report_table("WT vs mutant against the replicate noise floor ($(length(reps)) replicate(s))",
    ["metric", "WT vs mutant", "replicate mean", "replicate sd", "z", "verdict"], rows;
    md=MD, tsv=joinpath(OUT, "cb1_noise_floor.tsv"))
length(reps) < 3 && println("\n  [NOTE] fewer than 3 replicates: no metric can be declared significant.")

# ============================================================
# 4. cutoff sensitivity
# ============================================================
cuts = SEL == :ca ? [10.0, 12.0, 15.0] : [6.0, 7.0, 8.0, 10.0]
crow = Vector{Any}[]
for c in cuts
    a = solve(WT.coords,  M_wt,  c)
    b = solve(MUT.coords, M_mut, c)
    cc = compare(WT, MUT, a, b)
    push!(crow, [c, a.nz, b.nz, cc.reg, cc.tgt, cc.maxrel, cc.rmsip])
end
report_table("Cutoff sensitivity of the WT-vs-mutant result",
    ["cutoff (Å)", "zero modes WT", "zero modes mutant", "ΔMSF region (rel.)", "ΔMSF target (rel.)",
     "max abs(ΔMSF) (rel.)", "RMSIP"], crow; md=MD, tsv=joinpath(OUT, "cb1_cutoff_sensitivity.tsv"))

# ============================================================
# 5. stiffness-perturbation scan on the WT network (mutation-independent)
# ============================================================
println("\n=== Stiffness-perturbation scan (all springs of one residue scaled by $(PRS_FACTOR)) ===")
reslist = sort(unique(WT.resi))
scan    = sort(unique(vcat(reslist[1:PRS_STRIDE:end], TARGET)))
reg_at  = findall(r -> REGION[1] <= r <= REGION[2], WT.resi)
isempty(reg_at) && error("no atoms in the region of interest $(REGION)")
base    = S_wt.msf
infl    = Dict{Int,Float64}(); selfe = Dict{Int,Float64}()
tprofile = Float64[]
for (n, r) in enumerate(scan)
    pw = (i, j) -> (WT.resi[i] == r || WT.resi[j] == r) ? PRS_FACTOR : 1.0
    m  = solve(WT.coords, M_wt, CUTOFF; pairweight=pw).msf
    rel = (m .- base) ./ base
    infl[r]  = mean(abs.(rel[reg_at]))
    selfe[r] = mean(rel[findall(==(r), WT.resi)])
    r == TARGET && (tprofile = rel)
    n % 25 == 0 && (print("  ", n, "/", length(scan), "\r"); flush(stdout))
end
outside  = [r for r in scan if !(REGION[1] <= r <= REGION[2])]
val_t    = infl[TARGET]
pct_all  = 100 * count(r -> infl[r] <= val_t, scan) / length(scan)
pct_out  = TARGET in outside ? 100 * count(r -> infl[r] <= val_t, outside) / length(outside) : NaN
top = sort(outside, by = r -> -infl[r])[1:min(10, length(outside))]
report_table("Residues outside the region whose stiffness most affects the region's fluctuations",
    ["rank", "residue", "name", "mean abs(Δ MSF) in region (rel.)", "own Δ MSF (rel.)"],
    [[i, r, resname_of(WT, r), infl[r], selfe[r]] for (i, r) in enumerate(top)];
    md=MD, tsv=joinpath(OUT, "cb1_perturbation_top.tsv"))
report_table("Target residue $(TARGET) in the perturbation ranking",
    ["quantity", "value"],
    [["mean abs(Δ MSF) in region when residue $(TARGET) is softened", val_t],
     ["percentile among all scanned residues", pct_all],
     ["percentile among scanned residues outside the region", pct_out],
     ["residues scanned", length(scan)]]; md=MD)
save_tsv(joinpath(OUT, "cb1_perturbation_all.tsv"),
         ["residue", "name", "region_influence", "own_dMSF_rel"],
         [[r, resname_of(WT, r), infl[r], selfe[r]] for r in scan])

# ============================================================
# 6. optional B-factor calibration (WT, Cα model only)
# ============================================================
if !isempty(BFAC_PDB) && SEL == :ca
    if !isfile(BFAC_PDB)
        println("\n[WARN] CB1_BFAC_PDB not found: $BFAC_PDB")
    else
        bd = Dict{Int,Float64}()
        for l in eachline(BFAC_PDB)
            (startswith(l, "ATOM") && length(l) >= 66) || continue
            strip(l[13:16]) == "CA" && l[22] == CHAIN && (l[17] == ' ' || l[17] == 'A') || continue
            r = parse(Int, strip(l[23:26])); haskey(bd, r) || (bd[r] = parse(Float64, strip(l[61:66])))
        end
        b = [get(bd, r, NaN) for r in WT.resi]
        ok = isfinite.(b) .& (b .> 0)
        μ, σ = mean(b[ok]), std(b[ok])
        m = ok .& (b .< μ + 2σ)
        m[1:min(5, end)] .= false; m[max(1, end-4):end] .= false
        g = S_wt.msf
        s = sum(g[m] .* b[m]) / sum(abs2, g[m])
        γ = (8π^2 / 3) * KT / s
        println("\n=== B-factor calibration (WT, original experimental PDB) ===")
        @printf("  residues used: %d / %d   Pearson r = %.3f   Spearman ρ = %.3f   γ = %.3f %s\n",
                count(m), length(b), cor(g[m], b[m]), spearman(g[m], b[m]), γ,
                MODEL == :dw ? "kJ/mol" : "kJ/(mol·Å²)")
        println("  caution: crystal B-factors of CB1 include lattice contacts, the fusion partner and " *
                "lipidic-cubic-phase disorder; treat γ as an order-of-magnitude scale.")
        @printf("  region mean MSF = %.3f Å² (physical scale, calibrated γ)\n",
                KT * mean(g[findall(r -> REGION[1] <= r <= REGION[2], WT.resi)]) / γ)
    end
elseif !isempty(BFAC_PDB)
    println("\n(B-factor calibration is only implemented for CB1_SEL=ca)")
end

# ============================================================
# 7. per-residue output
# ============================================================
save_tsv(joinpath(OUT, "cb1_dmsf_wt_vs_mut.tsv"),
         ["residue", "name", "MSF_WT", "MSF_mutant", "dMSF_rel"],
         [[r, resname_of(WT, r), wm.rA[i], wm.rB[i], wm.rel[i]] for (i, r) in enumerate(wm.res)])
save_tsv(joinpath(OUT, "cb1_perturbation_target_profile.tsv"),
         ["atom_index", "residue", "dMSF_rel_when_target_softened"],
         [[i, WT.resi[i], tprofile[i]] for i in eachindex(tprofile)])

println("\nwritten to ", OUT, ":")
println("  cb1_summary.md, cb1_noise_floor.tsv, cb1_cutoff_sensitivity.tsv, cb1_perturbation_*.tsv,")
println("  cb1_dmsf_wt_vs_mut.tsv")

# ============================================================
# 8. factual read-out and limits
# ============================================================
println("\n", "="^66, "\nREAD-OUT\n", "="^66)
preflight_ok || println("• pre-flight problems were reported above: fix them before reading anything below.")
wm.identical && println("• the model cannot distinguish WT and mutant here (identical coordinates).")
if length(reps) < 3
    println("• fewer than 3 replicates: differences cannot be separated from protocol noise.")
elseif any(==("outside noise floor"), verdicts) || any(==("above noise"), verdicts)
    println("• metrics outside the noise floor: ",
            join([r[1] for r in rows if r[6] in ("outside noise floor",)], "; "))
else
    println("• no metric is outside the replicate noise floor.")
end
println("• the stiffness-perturbation scan is independent of the mutant structure: it shows how ",
        "strongly residue $(TARGET)'s environment couples to the region (percentile ",
        @sprintf("%.0f", pct_out), " among residues outside the region).")
println("• limits: harmonic, single-conformation, vacuum model; WT and mutant must come from the same ",
        "protocol; Cα models need changed coordinates, heavy-atom models need the cutoff checked.")