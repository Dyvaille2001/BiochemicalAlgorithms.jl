############################################################
# test_nma_correctness.jl
# Physics tests that FAIL if the model or machinery is wrong.
#
#   julia test_nma_correctness.jl
#
# Mutation checks (do these by hand once, confirm the named test fails):
#   * kfun γ/r² -> γ            : "diatomic DWENM" fails
#   * drop mass weighting       : "mass weighting" fails
#   * drop the rrᵀ/r² block     : "finite-difference DWENM" fails
############################################################
using Test, LinearAlgebra, Random, Statistics
include(joinpath(@__DIR__, "nma_core.jl"))
include(joinpath(@__DIR__, "nma_report.jl"))

Random.seed!(1)

@testset "NMA correctness" begin

    @testset "diatomic DWENM: λ = (γ/r²)(1/m₁ + 1/m₂)" begin
        r, γ, m1, m2 = 2.0, 3.0, 12.0, 16.0
        coords = [0.0 0 0; r 0 0]
        H = build_hessian(coords; cutoff=10.0, gamma=γ)
        vals, _, _ = compute_modes(H, [m1, m2])
        k = γ / r^2
        @test vals[end] ≈ k * (1/m1 + 1/m2) rtol=1e-12
        @test count_zero_modes(vals) == 5        # collinear: 5 zero modes
    end

    @testset "diatomic constant-spring ANM: λ = γ(1/m₁ + 1/m₂)" begin
        r, γ, m1, m2 = 2.0, 3.0, 12.0, 16.0
        coords = [0.0 0 0; r 0 0]
        H = build_hessian(coords, kfun_constant(γ); cutoff=10.0)
        vals, _, _ = compute_modes(H, [m1, m2])
        @test vals[end] ≈ γ * (1/m1 + 1/m2) rtol=1e-12
        @test count_zero_modes(vals) == 5
    end

    @testset "spring-law shorthand is consistent" begin
        coords = 6 .* rand(12, 3)
        γ = 2.5
        @test build_hessian(coords; cutoff=14.0, gamma=γ) ==
              build_hessian(coords, kfun_distance_weighted(γ); cutoff=14.0)
    end

    @testset "symmetry and rigid-body invariance" begin
        coords = 6 .* rand(30, 3); N = 30
        H = build_hessian(coords; cutoff=14.0)
        @test H ≈ H'
        for d in 1:3
            t = zeros(3N); t[d:3:end] .= 1.0
            @test norm(H * t) < 1e-10 * norm(H)
        end
        rotz = vec(permutedims(hcat(-coords[:, 2], coords[:, 1], zeros(N))))
        @test norm(H * rotz) < 1e-10 * norm(H) * norm(rotz)
        λ = eigvals(Symmetric(H))
        @test count_zero_modes(λ) == 6
        @test minimum(λ) > -1e-8 * maximum(λ)
    end

    @testset "finite-difference Hessian, DWENM" begin
        N = 6; γ = 1.0; cutoff = 12.0
        coords = 5 .* rand(N, 3)
        pairs = [(i, j) for i in 1:N-1 for j in i+1:N
                 if norm(coords[i,:] - coords[j,:]) < cutoff]
        @test !isempty(pairs)
        r0 = Dict(p => norm(coords[p[1],:] - coords[p[2],:]) for p in pairs)
        energy(x) = sum(0.5 * (γ / r0[p]^2) *
                        (norm(x[3p[1]-2:3p[1]] - x[3p[2]-2:3p[2]]) - r0[p])^2
                        for p in pairs)
        x0 = vec(permutedims(coords)); h = 1e-4
        Hfd = zeros(3N, 3N)
        for a in 1:3N, b in 1:3N
            f(sa, sb) = (x = copy(x0); x[a] += sa*h; x[b] += sb*h; energy(x))
            Hfd[a, b] = (f(1,1) - f(1,-1) - f(-1,1) + f(-1,-1)) / (4h^2)
        end
        @test Hfd ≈ build_hessian(coords; cutoff=cutoff, gamma=γ) atol=1e-5
    end

    @testset "finite-difference Hessian, constant spring" begin
        N = 6; γ = 2.0; cutoff = 12.0
        coords = 5 .* rand(N, 3)
        pairs = [(i, j) for i in 1:N-1 for j in i+1:N
                 if norm(coords[i,:] - coords[j,:]) < cutoff]
        @test !isempty(pairs)
        r0 = Dict(p => norm(coords[p[1],:] - coords[p[2],:]) for p in pairs)
        energy(x) = sum(0.5 * γ *
                        (norm(x[3p[1]-2:3p[1]] - x[3p[2]-2:3p[2]]) - r0[p])^2
                        for p in pairs)
        x0 = vec(permutedims(coords)); h = 1e-4
        Hfd = zeros(3N, 3N)
        for a in 1:3N, b in 1:3N
            f(sa, sb) = (x = copy(x0); x[a] += sa*h; x[b] += sb*h; energy(x))
            Hfd[a, b] = (f(1,1) - f(1,-1) - f(-1,1) + f(-1,-1)) / (4h^2)
        end
        @test Hfd ≈ build_hessian(coords, kfun_constant(γ); cutoff=cutoff) atol=1e-5
    end

    @testset "pair weights scale the stiffness of selected pairs" begin
        coords = 5 .* rand(8, 3)
        H0 = build_hessian(coords, kfun_constant(1.0); cutoff=20.0)
        H2 = build_hessian(coords, kfun_constant(1.0); cutoff=20.0,
                           pairweight=(i, j) -> 2.0)
        @test H2 ≈ 2 .* H0
        # weighting only pairs that contain atom 3 leaves other-only pairs unchanged
        Hw = build_hessian(coords, kfun_constant(1.0); cutoff=20.0,
                           pairweight=(i, j) -> (i == 3 || j == 3) ? 0.0 : 1.0)
        coords2 = coords[[1,2,4,5,6,7,8], :]
        H3 = build_hessian(coords2, kfun_constant(1.0); cutoff=20.0)
        keep = vcat(1:6, 10:24)                     # DOF of atoms 1,2,4..8
        @test Hw[keep, keep] ≈ H3
    end

    @testset "cutoff logic: rigid cluster + far atom" begin
        coords = [0.0 0 0; 4 0 0; 0 4 0; 0 0 4; 40 0 0]
        vals = eigvals(Symmetric(build_hessian(coords; cutoff=12.0)))
        @test count_zero_modes(vals) == 9          # 6 (cluster) + 3 (free atom)
        vals = eigvals(Symmetric(build_hessian(coords; cutoff=100.0)))
        @test count_zero_modes(vals) == 6
    end

    @testset "mass weighting solves H x = λ M x" begin
        coords = 8 .* rand(10, 3); masses = 10 .+ 6 .* rand(10)
        H = build_hessian(coords)
        vals, vecs, invM = compute_modes(H, masses)
        M = Diagonal(repeat(masses, inner=3))
        X = invM * vecs
        @test vecs' * vecs ≈ I atol=1e-10
        @test norm(H * X - M * X * Diagonal(vals)) < 1e-8 * norm(H)
        @test_throws ArgumentError compute_modes(H, vcat(masses[1:9], -1.0))
    end

    @testset "frequency conversion and tolerance" begin
        @test frequencies([1.0]; unit=:arb)   ≈ [1.0]
        @test frequencies([1.0]; unit=:cm1)[1]  ≈ 53.09      rtol=1e-3
        @test frequencies([1.0]; unit=:thz)[1]  ≈ 1.5915     rtol=1e-3
        @test frequencies([1.0]; unit=:radps)[1] ≈ 1.0e13    rtol=1e-12
        @test frequencies([1.0]; unit=:hz)[1]   ≈ 1.5915e12  rtol=1e-3
        @test frequencies([-4.0, 1e-20, 4.0]; unit=:arb) == [-2.0, 0.0, 2.0]
        @test frequencies([-1e-18, 1.0]; unit=:arb)[1] == 0.0
        @test frequencies(Float64[]; unit=:arb) == Float64[]
        @test frequencies([0.0, 0.0]; unit=:cm1) == [0.0, 0.0]
        @test_throws ArgumentError frequencies([1.0]; unit=:furlongs)
    end

    @testset "count_zero_modes edge cases" begin
        @test count_zero_modes(Float64[]) == 0
        @test count_zero_modes(zeros(9)) == 9
        @test count_zero_modes([0.0, 0.0, 1.0]) == 2
    end

    @testset "covariance: MSF ∝ 1/γ" begin
        coords = 8 .* rand(20, 3)
        H1 = build_hessian(coords; gamma=1.0)
        H2 = build_hessian(coords; gamma=2.0)
        @test msf_from_cov(covariance(H1)) ≈ 2 .* msf_from_cov(covariance(H2)) rtol=1e-8
    end

    @testset "covariance_mw equals covariance for unit masses" begin
        coords = 8 .* rand(15, 3)
        H = build_hessian(coords; cutoff=14.0)
        C1 = covariance(H; nzero=6, kT=1.0)
        C2 = covariance_mw(H, ones(15); nzero=6, kT=1.0)
        @test C1 ≈ C2 rtol=1e-8 atol=1e-12
    end

    @testset "DCCM properties" begin
        coords = 8 .* rand(20, 3)
        masses = fill(12.011, 20)
        C = covariance_mw(build_hessian(coords; cutoff=14.0), masses; kT=1.0)
        D = dccm_from_cov(C)
        @test size(D) == (20, 20)
        @test maximum(abs, diag(D) .- 1.0) < 1e-10
        @test maximum(abs, D - D') < 1e-10
        @test minimum(D) >= -1.0 - 1e-12 && maximum(D) <= 1.0 + 1e-12
    end

    @testset "rmsip and overlap sanity" begin
        rng = MersenneTwister(2)
        Q = Matrix(qr(randn(rng, 30, 20)).Q)
        V = Q[:, 1:10]
        W = Q[:, 11:20]
        @test rmsip(V, V, 10) ≈ 1.0 atol=1e-10
        @test rmsip(V, W, 10) ≈ 0.0 atol=1e-10
        O = overlap_matrix(V, W, 10)
        @test size(O) == (10, 10)
        @test maximum(O) < 1e-10
        @test overlap_matrix(V, V, 10) ≈ Matrix(I, 10, 10) atol=1e-10
        @test rmsip(V, V, 100) ≈ 1.0 atol=1e-10    # k is clamped
    end

    @testset "gamma calibration recovers a known gamma" begin
        coords = 10 .* rand(40, 3)
        masses = fill(12.011, 40)
        H1 = build_hessian(coords; gamma=1.0)
        γtrue, kT = 7.0, 2.494

        g = msf_from_cov(covariance(H1; kT=1.0))
        bexp = (8π^2 / 3) * kT .* g ./ γtrue
        γ, r, _ = calibrate_gamma(H1, bexp; kT=kT)
        @test γ ≈ γtrue rtol=1e-8
        @test r ≈ 1.0 atol=1e-10

        gm = msf_from_cov(covariance_mw(H1, masses; kT=1.0))
        bexp2 = (8π^2 / 3) * kT .* gm ./ γtrue
        γ2, r2, _ = calibrate_gamma(H1, bexp2; kT=kT, weighted=true, masses=masses)
        @test γ2 ≈ γtrue rtol=1e-8
        @test r2 ≈ 1.0 atol=1e-10
    end

    @testset "gamma calibration input guards" begin
        coords = 10 .* rand(20, 3)
        H1 = build_hessian(coords; cutoff=14.0)
        @test_throws DimensionMismatch calibrate_gamma(H1, ones(19))
        # non-positive and non-finite B-factors are masked, not fatal
        b = 1.0 .+ rand(20); b[1] = -5.0; b[2] = NaN     # varying B, so r is defined
        γ, r, _ = calibrate_gamma(H1, b)
        @test isfinite(γ) && isfinite(r)
        @test_throws ArgumentError calibrate_gamma(H1, fill(NaN, 20))
    end

    @testset "ranks and spearman handle ties" begin
        @test ranks([10.0, 20.0, 20.0, 30.0]) == [1.0, 2.5, 2.5, 4.0]
        @test spearman([1.0, 2.0, 3.0, 4.0], [1.0, 4.0, 9.0, 16.0]) ≈ 1.0
        @test spearman([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0]) ≈ -1.0
        @test isnan(spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))
    end

    @testset "write_mode_pdb RMS amplitude" begin
        N = 12
        coords = rand(N, 3) .* 10
        disp = randn(3N)
        names  = fill("CA", N)
        resn   = fill("ALA", N)
        chains = fill('A', N)
        resnum = collect(1:N)
        path = joinpath(mktempdir(), "mode.pdb")
        rms_amp = 1.5
        write_mode_pdb(path, coords, disp, names, resn, chains, resnum;
                       n_frames=4, rms_amp=rms_amp)

        lines = readlines(path)
        first_model = findfirst(l -> startswith(l, "ATOM"), lines)
        xyz = zeros(N, 3)
        for i in 1:N
            l = lines[first_model + i - 1]
            xyz[i, 1] = parse(Float64, l[31:38])
            xyz[i, 2] = parse(Float64, l[39:46])
            xyz[i, 3] = parse(Float64, l[47:54])
        end
        rms = sqrt(sum(abs2, xyz .- coords) / N)
        @test rms ≈ rms_amp atol=2e-3          # frame 1 has sin = 1
    end

    @testset "report table formatting" begin
        @test fmt(3) == "3"
        @test fmt(0.5) == "0.5000"
        @test fmt(6.0) == "6"                  # integer-valued float (promoted count)
        @test fmt(-0.0) == "0"
        @test fmt(12.5) == "12.5000"
        @test fmt(1.0e-7) == "1.000e-07"
        @test fmt(NaN) == "nan"
        @test fmt(true) == "yes"
        io = IOBuffer()
        print_table(io, "T", ["a", "b"], [[1, 2.5], [10, 0.25]])
        s = String(take!(io))
        @test occursin("| a  | b      |", s)
        @test count(==('\n'), s) == 6
    end

    @testset "mode shapes are independent of gamma" begin
        coords = 8 .* rand(25, 3)
        masses = fill(12.011, 25)
        _, V1, _ = compute_modes(build_hessian(coords; cutoff=14.0, gamma=1.0), masses)
        _, V7, _ = compute_modes(build_hessian(coords; cutoff=14.0, gamma=7.0), masses)
        @test rmsip(V1[:, 7:end], V7[:, 7:end], 10) ≈ 1.0 atol=1e-8
    end
end