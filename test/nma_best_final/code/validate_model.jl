############################################################
# validate_model.jl — model-level validation of the distance-weighted ENM
# (DWENM, k_ij = γ/r²) on hen egg-white lysozyme (1AKI)
#
# Pillars
#   A  experimental   Cα B-factors, γ calibration, baseline predictor
#   B  independent    ProDy ANM: constant-spring AND distance-weighted (optional files)
#   C  robustness     cutoff sweep, γ-invariance of mode shapes
#   D  falsification  shuffled-B null control
#   E  consistency    zero modes, graph, DCCM, frequencies
#
# ONE model throughout: mass-weighted Hessian (residue masses), Cartesian covariance.
# Writes validation_model.md / .tsv tables; exits non-zero if any check fails.
#
# Run (from the project root):
#   julia --project=. test\nma_best_final\code\validate_model.jl
# Optional ProDy references are read from ../prody; input and model outputs go in ../1aki:
#   prody_const_modes.txt prody_const_eigvals.txt prody_dw_modes.txt prody_dw_eigvals.txt
############################################################
using BiochemicalAlgorithms, LinearAlgebra, Statistics, Printf, Random, Downloads
include(joinpath(@__DIR__, "nma_core.jl"))
include(joinpath(@__DIR__, "nma_report.jl"))

const DATA_ROOT = normpath(joinpath(@__DIR__, ".."))
const ONEAKI_DIR = joinpath(DATA_ROOT, "1aki")
const PRODY_DIR = joinpath(DATA_ROOT, "prody")
const RESULTS_DIR = joinpath(DATA_ROOT, "validation_model")
mkpath(ONEAKI_DIR)
mkpath(RESULTS_DIR)

const KT         = 2.494          # kJ/mol at 300 K
const CUTOFF_REF = 12.0
const N_MODES    = 20
const SWEEP      = [8.0, 10.0, 12.0, 14.0, 15.0, 18.0]
const NZERO      = 6
const TOL_ZERO   = 1e-6
const NSHUF      = 1000
const MD         = joinpath(RESULTS_DIR, "validation_model.md")
isfile(MD) && rm(MD)

Random.seed!(20261028)

# ============================================================
# bookkeeping: every check becomes a row of the summary table
# ============================================================
const VROWS    = Vector{Vector{String}}()     # [pillar, check, value, criterion, status]
const FAILURES = String[]
const NOTES    = String[]

function check(pillar, name, cond::Bool; value="", criterion="")
    push!(VROWS, [pillar, name, string(value), criterion, cond ? "PASS" : "FAIL"])
    println("  [", cond ? "PASS" : "FAIL", "] ", name,
            isempty(string(value)) ? "" : "  ($value)")
    cond || push!(FAILURES, name)
    return cond
end

function info(pillar, name; value="", criterion="")
    push!(VROWS, [pillar, name, string(value), criterion, "INFO"])
    println("  [info] ", name, isempty(string(value)) ? "" : "  ($value)")
end

note(msg) = (println("  [note] ", msg); push!(NOTES, msg); nothing)

"B-factor calibration: B_pred = s·g, γ = (8π²/3) kT / s (γ in kJ/mol for k = γ/r²)."
function calibrate(g::AbstractVector, b::AbstractVector; kT=KT)
    s = sum(g .* b) / sum(abs2, g)
    return (8π^2 / 3) * kT / s, s
end

"Strip terminal variants: \"ALA-M\" -> \"ALA\"."
function base_resname(r::AbstractString)
    m = match(r"^([A-Z]{3})", uppercase(strip(r)))
    return m === nothing ? uppercase(strip(r)) : String(m.captures[1])
end

"Cα B-factors by residue number straight from the PDB text (first alt-loc only)."
function read_ca_bfactors(path)
    d = Dict{Int,Float64}(); chains = Set{Char}()
    for l in eachline(path)
        (startswith(l, "ATOM") && length(l) >= 66) || continue
        strip(l[13:16]) == "CA" || continue
        (l[17] == ' ' || l[17] == 'A') || continue
        resi = parse(Int, strip(l[23:26]))
        push!(chains, l[22])
        haskey(d, resi) && continue
        d[resi] = parse(Float64, strip(l[61:66]))
    end
    return d, chains
end

"Pair count, mean degree and component count for a cutoff."
function graph_stats(coords::AbstractMatrix, cutoff::Real)
    N = size(coords, 1)
    c2  = cutoff^2
    adj = [Int[] for _ in 1:N]
    np  = 0
    for i in 1:N-1, j in i+1:N
        if sum(abs2, coords[i, :] .- coords[j, :]) < c2
            push!(adj[i], j); push!(adj[j], i); np += 1
        end
    end
    seen = falses(N); ncomp = 0
    for s in 1:N
        seen[s] && continue
        ncomp += 1
        stack = [s]; seen[s] = true
        while !isempty(stack)
            u = pop!(stack)
            for v in adj[u]
                seen[v] || (seen[v] = true; push!(stack, v))
            end
        end
    end
    return (npairs=np, mean_degree=2np/N, ncomp=ncomp)
end

# ============================================================
# 1. structure
# ============================================================
pdb_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : joinpath(ONEAKI_DIR, "1aki.pdb")
if !isfile(pdb_path)
    println("Downloading 1AKI.pdb ...")
    Downloads.download("https://files.rcsb.org/download/1AKI.pdb", pdb_path)
end

mol = load_pdb(pdb_path, Float64)
ca = filter(collect(atoms(mol))) do a
    f = parent_fragment(a)
    !isnothing(f) && is_amino_acid(f) && strip(a.name) == "CA"
end
N = length(ca)
println("\n=== Structure ===")
println("  file     : ", basename(pdb_path))
println("  Cα atoms : ", N)

coords    = zeros(Float64, N, 3)
res_names = Vector{String}(undef, N)
res_nums  = Vector{Int}(undef, N)
for (i, a) in enumerate(ca)
    coords[i, :] .= a.r
    f = parent_fragment(a)
    res_names[i] = f.name
    res_nums[i]  = f.number
end
base_names = base_resname.(res_names)
chains     = fill('A', N)
atom_names = fill("CA", N)

bd, chain_set = read_ca_bfactors(pdb_path)
b_exp = [get(bd, r, NaN) for r in res_nums]

check("Input", "Cα count == 129 (1AKI)", N == 129; value=N, criterion="129")
check("Input", "single chain", length(chain_set) == 1; value=join(sort(collect(chain_set))), criterion="1 chain")
check("Input", "no duplicate residue numbers", length(unique(res_nums)) == N)
check("Input", "all B-factors present", !any(isnan, b_exp))
check("Input", "all B-factors positive", all(>(0), filter(!isnan, b_exp));
      value=@sprintf("min %.2f", minimum(filter(!isnan, b_exp))), criterion="> 0")
unknown = setdiff(unique(base_names), collect(keys(RESIDUE_MASSES)))
check("Input", "all residues have a mass", isempty(unknown);
      value=isempty(unknown) ? "$(length(unique(base_names))) types" : "unknown: $unknown")
masses = [RESIDUE_MASSES[r] for r in base_names]

# ============================================================
# 2. Pillar A — B-factors
# ============================================================
println("\n=== Pillar A: B-factor validation (cutoff $CUTOFF_REF Å) ===")
H1 = build_hessian(coords; cutoff=CUTOFF_REF, gamma=1.0)
check("E consistency", "Hessian symmetric", norm(H1 - H1') < 1e-10 * norm(H1))
nz_all = count_zero_modes(eigvals(Symmetric(H1)); tol=TOL_ZERO)
check("E consistency", "exactly 6 zero modes", nz_all == NZERO; value=nz_all, criterion="6")
gs = graph_stats(coords, CUTOFF_REF)
check("E consistency", "Cα graph connected", gs.ncomp == 1; value="$(gs.ncomp) component(s)", criterion="1")
println("  graph    : $(gs.npairs) pairs, ⟨degree⟩ = $(round(gs.mean_degree, digits=2))")

g_unw = msf_from_cov(covariance(H1; nzero=NZERO, kT=1.0))
g     = msf_from_cov(covariance_mw(H1, masses; kT=1.0, nzero=NZERO))
γ_cal, s_fit = calibrate(g, b_exp)
r_all = cor(g, b_exp)
rho   = spearman(g, b_exp)
r_unw = cor(g_unw, b_exp)

@printf("  γ            = %.4f kJ/mol            (k = γ/r²)\n", γ_cal)
@printf("  k_eff(4 Å)   = %.3f   k_eff(12 Å) = %.3f  kJ/(mol·Å²)\n", γ_cal / 4^2, γ_cal / 12^2)
@printf("  Pearson r    = %.4f   (unweighted model: %.4f)\n", r_all, r_unw)
@printf("  Spearman ρ   = %.4f\n", rho)

info("A B-factors", "calibrated γ [kJ/mol]"; value=@sprintf("%.3f", γ_cal))
info("A B-factors", "k_eff(4 Å) / k_eff(12 Å) [kJ/(mol·Å²)]";
     value=@sprintf("%.3f / %.3f", γ_cal/16, γ_cal/144))
info("A B-factors", "Pearson r"; value=@sprintf("%.4f", r_all))
info("A B-factors", "Spearman ρ"; value=@sprintf("%.4f", rho))
# sanity bands (not derived criteria: typical literature range for lysozyme-sized proteins)
check("A B-factors", "Pearson r in sanity band", 0.40 <= r_all <= 0.85;
      value=@sprintf("%.4f", r_all), criterion="0.40–0.85")
check("A B-factors", "Spearman ρ in sanity band", 0.40 <= rho <= 0.90;
      value=@sprintf("%.4f", rho), criterion="0.40–0.90")
check("A B-factors", "Pearson and Spearman consistent", abs(r_all - rho) < 0.20;
      value=@sprintf("Δ = %.3f", abs(r_all - rho)), criterion="< 0.20")
check("A B-factors", "MSF profile insensitive to the mass model", abs(r_all - r_unw) < 0.05;
      value=@sprintf("%.4f vs %.4f", r_all, r_unw), criterion="Δr < 0.05")

μ, σ = mean(b_exp), std(b_exp)
mask_flex = trues(N)
mask_flex[1:min(5, N)]   .= false
mask_flex[max(1, N-4):N] .= false
mask_flex .&= (b_exp .< μ + 2σ)
n_mask   = count(mask_flex)
r_masked = n_mask >= 10 ? cor(g[mask_flex], b_exp[mask_flex]) : NaN
flex_hi  = b_exp .>= μ + σ
n_hi     = count(flex_hi)
r_hi     = n_hi >= 10 ? cor(g[flex_hi], b_exp[flex_hi]) : NaN
@printf("  masked r     = %.4f   (%d / %d residues; termini and B > μ+2σ removed)\n", r_masked, n_mask, N)
@printf("  r (B ≥ μ+σ)  = %.4f   (%d residues)\n", r_hi, n_hi)
check("A B-factors", "masked r ≥ 0.35", !isnan(r_masked) && r_masked >= 0.35;
      value=@sprintf("%.4f", r_masked), criterion="≥ 0.35")
info("A B-factors", "r on the high-B residues (B ≥ μ+σ)"; value=@sprintf("%.4f", r_hi))
if r_masked < r_all
    note(@sprintf("masked r (%.3f) < full r (%.3f): the removed residues carry much of the signal. Report both.", r_masked, r_all))
end
if !isnan(r_hi) && r_hi < r_all
    note("r on the most flexible residues ($(round(r_hi, digits=3))) is below the overall r " *
         "($(round(r_all, digits=3))): high-B regions are dominated by static disorder / " *
         "multi-basin motion that a harmonic single-structure ENM does not describe.")
end

# ---- baseline predictor ------------------------------------------------
cn = [count(j -> j != i && sum(abs2, coords[i, :] .- coords[j, :]) < CUTOFF_REF^2, 1:N) for i in 1:N]
r_cn, rho_cn = cor(1.0 ./ cn, b_exp), spearman(1.0 ./ cn, b_exp)
@printf("  baseline 1/CN: Pearson r = %.4f   Spearman ρ = %.4f\n", r_cn, rho_cn)
info("A B-factors", "baseline: inverse neighbour count (12 Å), Pearson r"; value=@sprintf("%.4f", r_cn))
if r_all > r_cn + 0.02
    note(@sprintf("the DWENM (r = %.3f) beats the packing baseline (r = %.3f).", r_all, r_cn))
else
    note("the DWENM (r = $(round(r_all, digits=3))) does NOT clearly beat simple packing " *
         "(r = $(round(r_cn, digits=3))): B-factor agreement alone does not show that the elastic " *
         "network adds information beyond local density. Say so in the thesis.")
end

# ============================================================
# 3. Pillar D — null control
# ============================================================
println("\n=== Pillar D: null control (shuffled B-factors, $NSHUF permutations) ===")
r_null = [cor(g, shuffle(b_exp)) for _ in 1:NSHUF]
n_ge   = count(>=(r_all), r_null)
p_emp  = (n_ge + 1) / (NSHUF + 1)
@printf("  shuffled r : mean %+.4f | max %+.4f | exceedances %d | p = %.4f\n",
        mean(r_null), maximum(r_null), n_ge, p_emp)
check("D null", "shuffled r centred on 0", abs(mean(r_null)) < 0.15;
      value=@sprintf("mean %+.4f", mean(r_null)), criterion="|mean| < 0.15")
check("D null", "observed r beats the null", p_emp < 0.01;
      value=@sprintf("p = %.4f (floor %.4f)", p_emp, 1/(NSHUF+1)), criterion="p < 0.01")

# ============================================================
# 4. modes, frequencies, DCCM
# ============================================================
println("\n=== Modes and frequencies (calibrated γ) ===")
Cγ     = covariance_mw(H1, masses; kT=KT, nzero=NZERO)
b_calc = (8π^2 / 3) .* msf_from_cov(Cγ) ./ γ_cal      # zero-intercept calibration (as before)
dr = std(b_calc) / std(b_exp)
@printf("  dynamic range std(b_calc)/std(b_exp) = %.3f   (zero-intercept fit)\n", dr)
@printf("  mean(b_calc) / mean(b_exp)            = %.3f\n", mean(b_calc) / mean(b_exp))

# Affine calibration  b_exp ≈ B0 + s_aff·g : B0 is the baseline (lattice / static disorder,
# refinement restraints) that a harmonic Cα network cannot produce. With an intercept the
# std ratio equals r by construction, so the old "0.2–1.2" band tested the wrong thing:
# under a zero-intercept fit, dr = r · s_zero / s_affine exactly.
s_aff = cov(g, b_exp) / var(g)
B0    = mean(b_exp) - s_aff * mean(g)
γ_aff = (8π^2 / 3) * KT / s_aff
@printf("  affine fit   : B0 = %.2f Å²   γ_affine = %.1f kJ/mol   (zero-intercept γ = %.1f, ratio %.2f)\n",
        B0, γ_aff, γ_cal, γ_aff / γ_cal)
dr_id = r_all * s_fit / s_aff
check("E consistency", "scale identity dr = r·s₀/s_affine (code consistency)",
      isapprox(dr, dr_id; rtol=1e-6);
      value=@sprintf("%.4f vs %.4f", dr, dr_id), criterion="rtol 1e-6")
check("E consistency", "affine baseline B0 physical (0 ≤ B0 ≤ ⟨B⟩)", 0 <= B0 <= mean(b_exp);
      value=@sprintf("B0 = %.2f Å², ⟨B⟩ = %.2f Å²", B0, mean(b_exp)), criterion="0 ≤ B0 ≤ ⟨B⟩")
info("E consistency", "std(b_calc)/std(b_exp), zero-intercept fit"; value=@sprintf("%.3f", dr))
info("E consistency", "γ_affine / γ_zero-intercept"; value=@sprintf("%.2f", γ_aff / γ_cal))
rd(x, d=2) = round(x, digits=d)
note("zero-intercept fit gives std ratio $(rd(dr)) = r ($(rd(r_all, 3))) × s₀/s_affine ($(rd(s_fit / s_aff))): " *
     "the one-scalar fit has to reproduce the mean B (baseline B0 = $(rd(B0, 1)) Å² of ⟨B⟩ = $(rd(mean(b_exp), 1)) Å²) " *
     "and therefore over-disperses the profile. γ differs by a factor $(rd(γ_aff / γ_cal, 1)) between the two " *
     "calibrations, so absolute frequencies differ by $(rd(sqrt(γ_aff / γ_cal), 1))×. " *
     "Report γ as a range, or state the calibration used.")

# which calibration drives the mode table:  NMA_CAL=zero (default, as before) | affine
CAL   = lowercase(get(ENV, "NMA_CAL", "zero"))
γ_use = CAL == "affine" ? γ_aff : γ_cal
println("  frequencies below use the ", CAL == "affine" ? "affine" : "zero-intercept",
        " calibration (γ = ", round(γ_use, digits=2), " kJ/mol)")

vals, vecs, invM = compute_modes(H1 .* γ_use, masses)
freqs = frequencies(vals; unit=:cm1, tol=TOL_ZERO)
nz    = count_zero_modes(vals; tol=TOL_ZERO)
check("E consistency", "zero modes at calibrated γ", nz == NZERO; value=nz, criterion="6")
n_imag = count(<(0), freqs[nz+1:end])
check("E consistency", "no imaginary non-zero modes", n_imag == 0; value=n_imag, criterion="0")
f_first, f_last = freqs[nz+1], freqs[end]
@printf("  softest mode %.3f cm⁻¹ | highest mode %.3f cm⁻¹\n", f_first, f_last)
check("E consistency", "softest mode in sanity band", 1.0 <= f_first <= 15.0;
      value=@sprintf("%.2f cm⁻¹", f_first), criterion="1–15 cm⁻¹")
check("E consistency", "highest Cα mode ≤ 250 cm⁻¹", f_last <= 250.0;
      value=@sprintf("%.1f cm⁻¹", f_last), criterion="≤ 250 cm⁻¹")
note("frequencies are of the coarse-grained, B-factor-scaled network: collective-mode " *
     "frequencies, not vibrational frequencies of the protein.")

report_table("First 10 non-zero modes (calibrated γ)",
    ["mode", "eigenvalue (kJ/mol/amu/Å²)", "frequency (cm⁻¹)"],
    [[m, vals[m], freqs[m]] for m in nz+1:min(nz+10, length(freqs))]; md=MD)

println("\n=== DCCM ===")
D = dccm_from_cov(Cγ)
check("E consistency", "DCCM diagonal == 1", maximum(abs, diag(D) .- 1.0) < 1e-10)
check("E consistency", "DCCM symmetric", maximum(abs, D - D') < 1e-10)
check("E consistency", "DCCM entries in [-1, 1]", minimum(D) >= -1 - 1e-10 && maximum(D) <= 1 + 1e-10;
      value=@sprintf("[%.2f, %.2f]", minimum(D), maximum(D)))
off = .!Matrix{Bool}(I, N, N)
info("E consistency", "mean abs(off-diagonal) DCCM"; value=@sprintf("%.3f", mean(abs.(D[off]))))
open(joinpath(ONEAKI_DIR, "1aki_dccm.tsv"), "w") do io
    for i in 1:N
        println(io, join(round.(D[i, :], digits=4), '\t'))
    end
end

# ============================================================
# 5. Pillar C — γ-invariance and cutoff sweep
# ============================================================
println("\n=== Pillar C: γ-invariance (tautological: scaling a matrix keeps its eigenvectors) ===")
_, V1, _ = compute_modes(build_hessian(coords; cutoff=CUTOFF_REF, gamma=1.0), masses)
_, V7, _ = compute_modes(build_hessian(coords; cutoff=CUTOFF_REF, gamma=7.0), masses)
inv_r = rmsip(V1[:, NZERO+1:end], V7[:, NZERO+1:end], N_MODES)
check("C robustness", "mode shapes independent of γ (sanity)", inv_r > 1 - 1e-8;
      value=@sprintf("RMSIP = %.10f", inv_r), criterion="> 1 − 1e-8")

println("\n=== Pillar C: cutoff sweep ===")
sweep = NamedTuple[]
for c in SWEEP
    Hc  = build_hessian(coords; cutoff=c, gamma=1.0)
    nzc = count_zero_modes(eigvals(Symmetric(Hc)); tol=TOL_ZERO)
    gc  = graph_stats(coords, c)
    gmc = msf_from_cov(covariance_mw(Hc, masses; kT=1.0, nzero=NZERO))
    push!(sweep, (cutoff=c, npairs=gc.npairs, degree=gc.mean_degree, ncomp=gc.ncomp, nzero=nzc,
                  r=cor(gmc, b_exp), rho=spearman(gmc, b_exp),
                  r_masked=n_mask >= 10 ? cor(gmc[mask_flex], b_exp[mask_flex]) : NaN,
                  msf_corr=cor(gmc, g)))
end
report_table("Cutoff sweep (distance-weighted ENM)",
    ["cutoff (Å)", "pairs", "⟨degree⟩", "components", "zero modes", "Pearson r", "Spearman ρ",
     "masked r", "MSF corr vs 12 Å"],
    [[s.cutoff, s.npairs, s.degree, s.ncomp, s.nzero, s.r, s.rho, s.r_masked, s.msf_corr] for s in sweep];
    md=MD, tsv=joinpath(ONEAKI_DIR, "1aki_sweep.tsv"))

plateau = [s for s in sweep if 10.0 <= s.cutoff <= 15.0]
check("C robustness", "graph connected for 10–15 Å", all(s.ncomp == 1 for s in plateau))
check("C robustness", "6 zero modes for 10–15 Å", all(s.nzero == NZERO for s in plateau))
r_pl = [s.r for s in plateau]
check("C robustness", "r stable across 10–15 Å", maximum(r_pl) - minimum(r_pl) < 0.10;
      value=@sprintf("spread %.3f", maximum(r_pl) - minimum(r_pl)), criterion="< 0.10")
check("C robustness", "MSF profile stable across 10–15 Å", minimum(s.msf_corr for s in plateau) > 0.90;
      value=@sprintf("min corr %.3f", minimum(s.msf_corr for s in plateau)), criterion="> 0.90")
best = sweep[argmax([s.r for s in sweep])]
if best.cutoff < 10.0 || best.cutoff > 15.0
    note("the r maximum lies outside the 10–15 Å plateau ($(best.cutoff) Å): choose the " *
         "cutoff for graph robustness, not for maximal r, and say so.")
end

# ============================================================
# 6. Pillar B — ProDy (constant-spring AND distance-weighted)
# ============================================================
println("\n=== Pillar B: ProDy ANM cross-check ===")

"Read a ProDy matrix file with 3N lines x k columns -> (3N x k)."
function read_prody_matrix(path, nrows)
    rows = [parse.(Float64, split(l)) for l in eachline(path) if !isempty(strip(l))]
    length(rows) == nrows || error("expected $nrows lines in $path, got $(length(rows))")
    return permutedims(reduce(hcat, rows))
end

function prody_compare(tag, label, kfun, cutoff)
    Wf = joinpath(PRODY_DIR, "prody_$(tag)_modes.txt")
    Ef = joinpath(PRODY_DIR, "prody_$(tag)_eigvals.txt")
    if !(isfile(Wf) && isfile(Ef))
        info("B ProDy", "[$label] ProDy files missing — skipped"; value="run prody_anm_export.py")
        return
    end
    W  = read_prody_matrix(Wf, 3N)
    k  = min(N_MODES, size(W, 2))
    W  = W[:, 1:k]
    λP = [parse(Float64, strip(l)) for l in eachline(Ef) if !isempty(strip(l))][1:k]
    H  = build_hessian(coords, kfun; cutoff=cutoff)
    F  = eigen(Symmetric(H))
    nzj = count_zero_modes(F.values; tol=TOL_ZERO)
    check("B ProDy", "[$label] Julia Hessian has 6 zero modes", nzj == NZERO; value=nzj, criterion="6")
    λJ = F.values[nzj+1:nzj+k]
    VJ = F.vectors[:, nzj+1:nzj+k]
    orth = maximum(abs, W'W - I(k))
    check("B ProDy", "[$label] ProDy modes orthonormal", orth < 1e-6; value=@sprintf("%.1e", orth))
    relλ = maximum(abs.(λJ .- λP) ./ abs.(λP))
    r10  = rmsip(VJ, W, min(10, k); orth_warn=false)
    rk   = rmsip(VJ, W, k;          orth_warn=false)
    check("B ProDy", "[$label] eigenvalues agree", relλ < 1e-6;
          value=@sprintf("max rel. diff %.2e", relλ), criterion="< 1e-6")
    check("B ProDy", "[$label] RMSIP(10)", r10 > 0.999; value=@sprintf("%.6f", r10), criterion="> 0.999")
    check("B ProDy", "[$label] RMSIP($k)", rk > 0.999; value=@sprintf("%.6f", rk), criterion="> 0.999")
    O = overlap_matrix(VJ, W, k)
    open(joinpath(ONEAKI_DIR, "1aki_overlap_prody_$(tag).tsv"), "w") do io
        println(io, "i\tj\toverlap")
        for i in 1:k, j in 1:k
            @printf(io, "%d\t%d\t%.6f\n", i, j, O[i, j])
        end
    end
end
prody_compare("const", "constant spring, 15 Å", kfun_constant(1.0), 15.0)
prody_compare("dw",    "distance-weighted, 12 Å", kfun_distance_weighted(1.0), 12.0)

# ============================================================
# 7. outputs
# ============================================================
open(joinpath(ONEAKI_DIR, "1aki_bfactors.tsv"), "w") do io
    println(io, "ResNum\tResName\tB_exp\tB_calc\tmasked")
    for i in 1:N
        @printf(io, "%d\t%s\t%.3f\t%.3f\t%s\n", res_nums[i], base_names[i], b_exp[i], b_calc[i],
                mask_flex[i] ? "1" : "0")
    end
end
write_mode_pdb(joinpath(ONEAKI_DIR, "1aki_mode_$(nz+1).pdb"), coords, invM * vecs[:, nz+1],
               atom_names, base_names, chains, res_nums; n_frames=30, rms_amp=1.5)
println("\nwritten: 1aki_bfactors.tsv, 1aki_dccm.tsv, 1aki_sweep.tsv, 1aki_mode_$(nz+1).pdb")

# ============================================================
# 8. summary table
# ============================================================
report_table("VALIDATION SUMMARY — distance-weighted ENM on 1AKI",
    ["Pillar", "Check", "Value", "Criterion", "Status"], VROWS;
    md=MD, tsv=joinpath(RESULTS_DIR, "validation_model_summary.tsv"))
if !isempty(NOTES)
    println("\nNotes:")
    foreach(n -> println("  - ", n), NOTES)
    open(MD, "a") do io
        println(io, "\nNotes:")
        foreach(n -> println(io, "- ", n), NOTES)
    end
end
npass = count(r -> r[5] == "PASS", VROWS); nfail = count(r -> r[5] == "FAIL", VROWS)
println("\n", "="^60)
println(isempty(FAILURES) ? "ALL CHECKS PASSED ($npass passed)" :
        "$nfail CHECK(S) FAILED ($npass passed):")
foreach(f -> println("  - ", f), FAILURES)
println("="^60)
isempty(FAILURES) || exit(1)
