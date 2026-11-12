
First 10 non-zero modes (calibrated γ)
| mode | eigenvalue (kJ/mol/amu/Å²) | frequency (cm⁻¹) |
|------|----------------------------|------------------|
| 7    | 0.0050                     | 3.7648           |
| 8    | 0.0069                     | 4.4150           |
| 9    | 0.0102                     | 5.3611           |
| 10   | 0.0142                     | 6.3190           |
| 11   | 0.0160                     | 6.7191           |
| 12   | 0.0165                     | 6.8121           |
| 13   | 0.0182                     | 7.1590           |
| 14   | 0.0217                     | 7.8285           |
| 15   | 0.0247                     | 8.3371           |
| 16   | 0.0258                     | 8.5323           |

Cutoff sweep (distance-weighted ENM)
| cutoff (Å) | pairs | ⟨degree⟩ | components | zero modes | Pearson r | Spearman ρ | masked r | MSF corr vs 12 Å |
|------------|-------|----------|------------|------------|-----------|------------|----------|------------------|
| 8          | 632   | 9.7984   | 1          | 6          | 0.6122    | 0.6502     | 0.4711   | 0.8086           |
| 10         | 1131  | 17.5349  | 1          | 6          | 0.5649    | 0.6944     | 0.4499   | 0.9883           |
| 12         | 1792  | 27.7829  | 1          | 6          | 0.5653    | 0.7386     | 0.4750   | 1                |
| 14         | 2549  | 39.5194  | 1          | 6          | 0.5633    | 0.7554     | 0.4823   | 0.9917           |
| 15         | 2983  | 46.2481  | 1          | 6          | 0.5628    | 0.7557     | 0.4903   | 0.9856           |
| 18         | 4221  | 65.4419  | 1          | 6          | 0.6255    | 0.7420     | 0.5629   | 0.9540           |

VALIDATION SUMMARY — distance-weighted ENM on 1AKI
| Pillar        | Check                                                    | Value                         | Criterion     | Status |
|---------------|----------------------------------------------------------|-------------------------------|---------------|--------|
| Input         | Cα count == 129 (1AKI)                                   | 129                           | 129           | PASS   |
| Input         | single chain                                             | A                             | 1 chain       | PASS   |
| Input         | no duplicate residue numbers                             |                               |               | PASS   |
| Input         | all B-factors present                                    |                               |               | PASS   |
| Input         | all B-factors positive                                   | min 10.60                     | > 0           | PASS   |
| Input         | all residues have a mass                                 | 20 types                      |               | PASS   |
| E consistency | Hessian symmetric                                        |                               |               | PASS   |
| E consistency | exactly 6 zero modes                                     | 6                             | 6             | PASS   |
| E consistency | Cα graph connected                                       | 1 component(s)                | 1             | PASS   |
| A B-factors   | calibrated γ [kJ/mol]                                    | 173.306                       |               | INFO   |
| A B-factors   | k_eff(4 Å) / k_eff(12 Å) [kJ/(mol·Å²)]                   | 10.832 / 1.204                |               | INFO   |
| A B-factors   | Pearson r                                                | 0.5653                        |               | INFO   |
| A B-factors   | Spearman ρ                                               | 0.7386                        |               | INFO   |
| A B-factors   | Pearson r in sanity band                                 | 0.5653                        | 0.40–0.85     | PASS   |
| A B-factors   | Spearman ρ in sanity band                                | 0.7386                        | 0.40–0.90     | PASS   |
| A B-factors   | Pearson and Spearman consistent                          | Δ = 0.173                     | < 0.20        | PASS   |
| A B-factors   | MSF profile insensitive to the mass model                | 0.5653 vs 0.5649              | Δr < 0.05     | PASS   |
| A B-factors   | masked r ≥ 0.35                                          | 0.4750                        | ≥ 0.35        | PASS   |
| A B-factors   | r on the high-B residues (B ≥ μ+σ)                       | 0.2022                        |               | INFO   |
| A B-factors   | baseline: inverse neighbour count (12 Å), Pearson r      | 0.7352                        |               | INFO   |
| D null        | shuffled r centred on 0                                  | mean -0.0029                  | |mean| < 0.15 | PASS   |
| D null        | observed r beats the null                                | p = 0.0010 (floor 0.0010)     | p < 0.01      | PASS   |
| E consistency | scale identity dr = r·s₀/s_affine (code consistency)     | 2.4261 vs 2.4261              | rtol 1e-6     | PASS   |
| E consistency | affine baseline B0 physical (0 ≤ B0 ≤ ⟨B⟩)               | B0 = 14.51 Å², ⟨B⟩ = 17.58 Å² | 0 ≤ B0 ≤ ⟨B⟩  | PASS   |
| E consistency | std(b_calc)/std(b_exp), zero-intercept fit               | 2.426                         |               | INFO   |
| E consistency | γ_affine / γ_zero-intercept                              | 4.29                          |               | INFO   |
| E consistency | zero modes at calibrated γ                               | 6                             | 6             | PASS   |
| E consistency | no imaginary non-zero modes                              | 0                             | 0             | PASS   |
| E consistency | softest mode in sanity band                              | 3.76 cm⁻¹                     | 1–15 cm⁻¹     | PASS   |
| E consistency | highest Cα mode ≤ 250 cm⁻¹                               | 55.1 cm⁻¹                     | ≤ 250 cm⁻¹    | PASS   |
| E consistency | DCCM diagonal == 1                                       |                               |               | PASS   |
| E consistency | DCCM symmetric                                           |                               |               | PASS   |
| E consistency | DCCM entries in [-1, 1]                                  | [-0.22, 1.00]                 |               | PASS   |
| E consistency | mean abs(off-diagonal) DCCM                              | 0.065                         |               | INFO   |
| C robustness  | mode shapes independent of γ (sanity)                    | RMSIP = 1.0000000000          | > 1 − 1e-8    | PASS   |
| C robustness  | graph connected for 10–15 Å                              |                               |               | PASS   |
| C robustness  | 6 zero modes for 10–15 Å                                 |                               |               | PASS   |
| C robustness  | r stable across 10–15 Å                                  | spread 0.002                  | < 0.10        | PASS   |
| C robustness  | MSF profile stable across 10–15 Å                        | min corr 0.986                | > 0.90        | PASS   |
| B ProDy       | [constant spring, 15 Å] Julia Hessian has 6 zero modes   | 6                             | 6             | PASS   |
| B ProDy       | [constant spring, 15 Å] ProDy modes orthonormal          | 9.8e-12                       |               | PASS   |
| B ProDy       | [constant spring, 15 Å] eigenvalues agree                | max rel. diff 3.15e-13        | < 1e-6        | PASS   |
| B ProDy       | [constant spring, 15 Å] RMSIP(10)                        | 1.000000                      | > 0.999       | PASS   |
| B ProDy       | [constant spring, 15 Å] RMSIP(20)                        | 1.000000                      | > 0.999       | PASS   |
| B ProDy       | [distance-weighted, 12 Å] Julia Hessian has 6 zero modes | 6                             | 6             | PASS   |
| B ProDy       | [distance-weighted, 12 Å] ProDy modes orthonormal        | 1.0e-11                       |               | PASS   |
| B ProDy       | [distance-weighted, 12 Å] eigenvalues agree              | max rel. diff 3.86e-13        | < 1e-6        | PASS   |
| B ProDy       | [distance-weighted, 12 Å] RMSIP(10)                      | 1.000000                      | > 0.999       | PASS   |
| B ProDy       | [distance-weighted, 12 Å] RMSIP(20)                      | 1.000000                      | > 0.999       | PASS   |

Notes:
- masked r (0.475) < full r (0.565): the removed residues carry much of the signal. Report both.
- r on the most flexible residues (0.202) is below the overall r (0.565): high-B regions are dominated by static disorder / multi-basin motion that a harmonic single-structure ENM does not describe.
- the DWENM (r = 0.565) does NOT clearly beat simple packing (r = 0.735): B-factor agreement alone does not show that the elastic network adds information beyond local density. Say so in the thesis.
- zero-intercept fit gives std ratio 2.43 = r (0.565) × s₀/s_affine (4.29): the one-scalar fit has to reproduce the mean B (baseline B0 = 14.5 Å² of ⟨B⟩ = 17.6 Å²) and therefore over-disperses the profile. γ differs by a factor 4.3 between the two calibrations, so absolute frequencies differ by 2.1×. Report γ as a range, or state the calibration used.
- frequencies are of the coarse-grained, B-factor-scaled network: collective-mode frequencies, not vibrational frequencies of the protein.
- the r maximum lies outside the 10–15 Å plateau (18.0 Å): choose the cutoff for graph robustness, not for maximal r, and say so.
