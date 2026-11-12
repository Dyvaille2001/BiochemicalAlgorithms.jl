# CB1 investigation

- WT: [WT_min.pdb](../../WT/WT_min.pdb)
- mutant: [F238L_min.pdb](../../mutants/F238L_min.pdb)
- WT replicates (6): WT_rep2.pdb, WT_rep3.pdb, WT_rep4.pdb, WT_rep5.pdb, WT_rep6.pdb, WT_rep7.pdb
- mutant replicates (4): MUT_rep2.pdb, MUT_rep3.pdb, MUT_rep4.pdb, MUT_rep5.pdb
- model: dw, selection heavy, cutoff 7.0 Å; scan skipped


Pre-flight checks
| item                        | WT     | mutant |
|-----------------------------|--------|--------|
| atoms in the model          | 2305   | 2302   |
| residues                    | 289    | 289    |
| residue at target 237       | PHE    | LEU    |
| expected at target          | PHE    | LEU    |
| zero modes (expect 6)       | 6      | 6      |
| softest non-zero eigenvalue | 0.0016 | 0.0016 |

WT vs mutant against the replicate noise floor (6 replicate(s))
| metric                                 | WT vs mutant | replicate mean | replicate sd | z       | verdict             |
|----------------------------------------|--------------|----------------|--------------|---------|---------------------|
| relative ΔMSF, region mean             | 0.0041       | 5.640e-04      | 0.0023       | 1.5454  | within noise floor  |
| relative ΔMSF, target residue          | 0.1499       | -0.0195        | 0.0018       | 94.7776 | outside noise floor |
| max abs(relative ΔMSF)                 | 0.9997       | 0.7432         | 0.3199       | 0.8018  | within noise floor  |
| 1 − RMSIP (low-mode subspace)          | 0.0264       | 0.0315         | 0.0109       | -0.4691 | within noise floor  |
| mean abs(log λ ratio) (spectrum shift) | 0.0416       | 0.0344         | 0.0115       | 0.6267  | within noise floor  |
| ΔDCCM, Frobenius / N                   | 0.0019       | 0.0019         | 5.090e-04    | 0.0215  | within noise floor  |
| ΔDCCM, target × region mean            | -1.176e-04   | 4.618e-04      | 1.063e-04    | -5.4500 | outside noise floor |

Two-sample comparison: 6 WT replicates vs 5 mutant structures (each compared with WT)
| metric                                 | WT-rep mean | WT-rep sd | mutant mean | mutant sd | Δ mean     | p (permutation) | p floor | Holm p | verdict                       |
|----------------------------------------|-------------|-----------|-------------|-----------|------------|-----------------|---------|--------|-------------------------------|
| relative ΔMSF, region mean             | 5.640e-04   | 0.0023    | 0.0031      | 0.0018    | 0.0026     | 0.0671          | 0.0043  | 0.3355 | no detectable difference      |
| relative ΔMSF, target residue          | -0.0195     | 0.0018    | 0.1397      | 0.0101    | 0.1591     | 0.0022          | 0.0043  | 0.0152 | differs from replicate spread |
| max abs(relative ΔMSF)                 | 0.7432      | 0.3199    | 0.9142      | 0.4804    | 0.1709     | 0.2554          | 0.0022  | 0.7662 | no detectable difference      |
| 1 − RMSIP (low-mode subspace)          | 0.0315      | 0.0109    | 0.0322      | 0.0129    | 6.673e-04  | 0.5541          | 0.0022  | 0.7662 | no detectable difference      |
| mean abs(log λ ratio) (spectrum shift) | 0.0344      | 0.0115    | 0.0376      | 0.0112    | 0.0032     | 0.3268          | 0.0022  | 0.7662 | no detectable difference      |
| ΔDCCM, Frobenius / N                   | 0.0019      | 5.090e-04 | 0.0027      | 0.0012    | 7.540e-04  | 0.1126          | 0.0022  | 0.4502 | no detectable difference      |
| ΔDCCM, target × region mean            | 4.618e-04   | 1.063e-04 | 5.886e-05   | 1.581e-04 | -4.029e-04 | 0.0022          | 0.0043  | 0.0152 | differs from replicate spread |

Cutoff sensitivity of the WT-vs-mutant result
| cutoff (Å) | zero modes WT | zero modes mutant | ΔMSF region (rel.) | ΔMSF target (rel.) | max abs(ΔMSF) (rel.) | RMSIP  |
|------------|---------------|-------------------|--------------------|--------------------|----------------------|--------|
| 6          | 6             | 6                 | 0.0052             | 0.2139             | 1.6103               | 0.9728 |
| 7          | 6             | 6                 | 0.0041             | 0.1499             | 0.9997               | 0.9736 |
| 8          | 6             | 6                 | 0.0044             | 0.1076             | 0.8094               | 0.9738 |
| 10         | 6             | 6                 | 0.0025             | 0.0831             | 0.4123               | 0.9740 |
