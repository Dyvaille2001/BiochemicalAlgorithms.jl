# NMA validation runbook

Run commands from the `BiochemicalAlgorithms.jl` project root unless stated otherwise. The canonical scripts are in this folder; earlier exploratory variants are in `archive/` and are not part of the final validation pipeline.

## 1. Correctness checks

```powershell
julia --project=. test\nma_best_final\code\test_nma_correctness.jl
```

## 2. MMTK reference data and force-field validation

The reference frequencies, gradients, energies, and minimized structures in `../mmtk/` were produced with the Python 2 MMTK Docker environment. To regenerate them, mount `test/nma_best_final` as `/data`, then run inside the container:

```sh
cd /data/mmtk
python ../code/mmtk_compare.py bala1 bala1
python ../code/mmtk_compare.py AlaAla.pdb alaala
```

Supply the original `AlaAla.pdb` in `/data/mmtk`; it is not part of this result bundle. Keep each minimized PDB and its three MMTK text references from the same run because the PDB coordinates are rounded to three decimal places.

Run the Julia comparison from the project root:

```powershell
julia --project=. test\nma_best_final\code\validate_forcefield.jl
```

It reads references from `test/nma_best_final/mmtk/` and writes the report and TSV tables to `test/nma_best_final/validation_forcefield/`.

## 3. ProDy and 1AKI validation

```powershell
python test\nma_best_final\code\prody_anm_export.py test\nma_best_final\1aki\1aki.pdb
julia --project=. test\nma_best_final\code\validate_model.jl
```

ProDy reference modes are written to `prody/`. The model validator reads the 1AKI structure and ProDy files from their respective folders; it writes the report to `validation_model/` and residue/mode tables and structures to `1aki/`.

## 4. CB1 WT versus F238L

`cb1_investigation.jl` accepts WT and mutant PDB paths, then WT replicates followed by `--mut` and the mutant replicates. For the stored final set, from the project root:

```powershell
julia --project=. test\nma_best_final\code\cb1_investigation.jl `
  test\nma_best_final\WT\WT_min.pdb `
  test\nma_best_final\mutants\F238L_min.pdb `
  test\nma_best_final\WT\WT_rep2.pdb test\nma_best_final\WT\WT_rep3.pdb `
  test\nma_best_final\WT\WT_rep4.pdb test\nma_best_final\WT\WT_rep5.pdb `
  test\nma_best_final\WT\WT_rep6.pdb test\nma_best_final\WT\WT_rep7.pdb `
  --mut test\nma_best_final\mutants\MUT_rep2.pdb test\nma_best_final\mutants\MUT_rep3.pdb `
  test\nma_best_final\mutants\MUT_rep4.pdb test\nma_best_final\mutants\MUT_rep5.pdb
```

The default output directory is `cb1/final/`. Set `$env:CB1_OUT='rotA_t'` before running to write a separate rotamer analysis to `cb1/rotA_t/`.

## 5. Convenience runner

`powershell -ExecutionPolicy Bypass -File test\nma_best_final\code\run_all.ps1` runs the Julia correctness, force-field, and 1AKI validation steps and stores logs under `logs/`. It does not run Docker/ProDy preparation or the CB1 analysis.