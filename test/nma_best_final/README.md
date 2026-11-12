# NMA validation and CB1 investigation

This directory contains the final scripts, reference data, and generated results used for the thesis analyses. The organization keeps scripts and each result family together while keeping the `alpha_sweep*` files out of the GitHub version.

## Directory map

| Directory | Contents |
|---|---|
| `code/` | Canonical Julia, Python, and PowerShell scripts; earlier script variants are in `code/archive/`. |
| `validation_forcefield/` | Force-field versus MMTK report and TSV result tables. |
| `validation_model/` | Lysozyme model-validation report and summary table. |
| `mmtk/` | MMTK reference frequencies, gradients, energies, and minimized bala1/Ala-Ala structures. |
| `prody/` | ProDy reference eigenvalues and modes. |
| `1aki/` | 1AKI input structure, mode files, B-factor/DCCM tables, overlap comparisons, and figures. |
| `WT/` | Wild-type CB1 minimized structure, DMSF structure, and replicate structures. |
| `mutants/` | F238L and mutant replicate structures, including `5TGZ_F238L.pdb`. |
| `cb1/` | CB1 analysis tables and reports (`final/`, `rotA_t/`), the `5TGZ.pdb` structure, and figures. |
| `logs/` | Sanitized historical run transcript. The tables and reports above are the authoritative results. |

`alpha_sweep*` outputs are retained locally but intentionally ignored and not included in this repository publication. Temporary files and Windows `desktop.ini` are also excluded.

## Reproducing the principal validations

Run these commands from the `BiochemicalAlgorithms.jl` project root. They require the project dependencies in `Project.toml`.

```powershell
julia --project=. test\nma_best_final\code\test_nma_correctness.jl
julia --project=. test\nma_best_final\code\validate_forcefield.jl
python test\nma_best_final\code\prody_anm_export.py test\nma_best_final\1aki\1aki.pdb
julia --project=. test\nma_best_final\code\validate_model.jl
```

The model-validation script writes its report and summary under `validation_model/`, its 1AKI tables and mode PDB under `1aki/`, and reads optional ProDy references from `prody/`. The force-field validation script reads the MMTK files from `mmtk/` and writes its report/tables to `validation_forcefield/`.

For CB1, pass the WT and mutant structures and replicate PDBs to `code/cb1_investigation.jl`; results default to `cb1/final/`. Set `CB1_OUT=rotA_t` to write a separate rotamer run to `cb1/rotA_t/`. The stored CB1 summaries use relative structure names so they do not expose local machine paths.

## MMTK reference provenance

The MMTK reference files were generated in the project’s Python 2 Docker environment. The corresponding frequencies, gradients, energies, and minimized structures are grouped in `mmtk/`; each set should come from the same MMTK run because the minimized PDB files are rounded to three decimal places. The original Ala-Ala input PDB is not included in this results folder.
