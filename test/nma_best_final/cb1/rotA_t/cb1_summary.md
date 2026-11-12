# CB1 investigation

- WT: [WT_min.pdb](../../WT/WT_min.pdb)
- mutant: [F238L_rotA_min.pdb](../../mutants/F238L_rotA_min.pdb)
- WT replicates (6): WT_rep2.pdb, WT_rep3.pdb, WT_rep4.pdb, WT_rep5.pdb, WT_rep6.pdb, WT_rep7.pdb
- mutant replicates (0):
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
| metric                                 | WT vs mutant | replicate mean | replicate sd | z        | verdict             |
|----------------------------------------|--------------|----------------|--------------|----------|---------------------|
| relative ΔMSF, region mean             | -2.750e-04   | 5.640e-04      | 0.0023       | -0.3666  | within noise floor  |
| relative ΔMSF, target residue          | 0.1963       | -0.0195        | 0.0018       | 120.7155 | outside noise floor |
| max abs(relative ΔMSF)                 | 0.4963       | 0.7432         | 0.3199       | -0.7720  | within noise floor  |
| 1 − RMSIP (low-mode subspace)          | 0.0282       | 0.0315         | 0.0109       | -0.3006  | within noise floor  |
| mean abs(log λ ratio) (spectrum shift) | 0.0348       | 0.0344         | 0.0115       | 0.0351   | within noise floor  |
| ΔDCCM, Frobenius / N                   | 0.0031       | 0.0019         | 5.090e-04    | 2.3688   | within noise floor  |
| ΔDCCM, target × region mean            | 1.291e-04    | 4.618e-04      | 1.063e-04    | -3.1293  | outside noise floor |

Cutoff sensitivity of the WT-vs-mutant result
| cutoff (Å) | zero modes WT | zero modes mutant | ΔMSF region (rel.) | ΔMSF target (rel.) | max abs(ΔMSF) (rel.) | RMSIP  |
|------------|---------------|-------------------|--------------------|--------------------|----------------------|--------|
| 6          | 6             | 6                 | 0.0098             | 0.3097             | 0.4992               | 0.9695 |
| 7          | 6             | 6                 | -2.750e-04         | 0.1963             | 0.4963               | 0.9718 |
| 8          | 6             | 6                 | 0.0020             | 0.1876             | 0.3931               | 0.9725 |
| 10         | 6             | 6                 | 2.745e-04          | 0.1262             | 0.3018               | 0.9749 |
