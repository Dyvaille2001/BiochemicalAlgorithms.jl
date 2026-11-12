
Table 1 — energy terms at the SAME geometry, Julia − MMTK (kJ/mol)
| System | Improper convention | Atoms | Δ bond    | Δ angle    | Δ vdW      | Δ electrostatic | Δ torsion  | Δ total   | MMTK total |
|--------|---------------------|-------|-----------|------------|------------|-----------------|------------|-----------|------------|
| bala1  | amber               | 22    | 2.681e-12 | -3.412e-11 | 1.084e-09  | 2.004e-07       | 1.0649     | 1.0649    | -117.7560  |
| bala1  | mmtk                | 22    | 2.681e-12 | -3.412e-11 | 1.084e-09  | 2.004e-07       | 1.888e-10  | 2.016e-07 | -117.7560  |
| alaala | amber               | 23    | 1.344e-10 | -4.817e-12 | -6.757e-11 | 6.429e-07       | 51.8615    | 51.8615   | -373.2118  |
| alaala | mmtk                | 23    | 1.344e-10 | -4.817e-12 | -6.757e-11 | 6.429e-07       | -2.591e-11 | 6.429e-07 | -373.2118  |

Table 2 — gradients, force/energy consistency and minimum
| System | Improper convention | gradient RMS vs MMTK (kJ/mol/Å) | max force/energy ratio | |g| after polish | Hessian asymmetry | E_min Julia | E_min MMTK | Δ E_min   |
|--------|---------------------|---------------------------------|------------------------|------------------|-------------------|-------------|------------|-----------|
| bala1  | amber               | 10.7585                         | 5.955e-10              | 1.509e-07        | 3.633e-09         | -118.3127   | -117.7640  | -0.5486   |
| bala1  | mmtk                | 2.942e-08                       | 8.440e-10              | 1.915e-08        | 4.357e-09         | -117.7640   | -117.7640  | 2.007e-07 |
| alaala | amber               | 46.2399                         | 4.126e-10              | 1.497e-10        | 3.771e-09         | -336.4484   | -373.2207  | 36.7723   |
| alaala | mmtk                | 8.069e-08                       | 3.764e-10              | 1.406e-11        | 3.304e-09         | -373.2207   | -373.2207  | 6.437e-07 |

Table 3 — normal-mode frequencies vs MMTK (cm⁻¹)
| System | Improper convention | Modes | zero modes | imaginary | mean rel. diff | max rel. diff | low modes (<500) | high modes (>2500) | median ratio | lowest mode Julia | lowest mode MMTK |
|--------|---------------------|-------|------------|-----------|----------------|---------------|------------------|--------------------|--------------|-------------------|------------------|
| bala1  | amber               | 60    | 6          | 0         | 0.0349         | 0.3774        | 0.0970           | 2.154e-04          | 1.0002       | 38.2639           | 27.7788          |
| bala1  | mmtk                | 60    | 6          | 0         | 6.850e-06      | 1.176e-05     | 3.280e-06        | 1.093e-05          | 1.0000       | 27.7789           | 27.7788          |
| alaala | amber               | 63    | 6          | 0         | 0.0110         | 0.1395        | 0.0192           | 3.205e-04          | 1.0000       | 73.2888           | 69.3804          |
| alaala | mmtk                | 63    | 6          | 0         | 6.813e-06      | 1.185e-05     | 3.943e-06        | 1.109e-05          | 1.0000       | 69.3802           | 69.3804          |

Table 4 — pass/fail with MMTK's improper convention imposed
| System | Criterion                    | Value     | Status |
|--------|------------------------------|-----------|--------|
| bala1  | |Δ E_total| < 0.0001         | 2.016e-07 | PASS   |
| bala1  | gradient RMS < 1.0e-5        | 2.942e-08 | PASS   |
| bala1  | |Δ E_min| < 0.0001           | 2.007e-07 | PASS   |
| bala1  | mean rel. freq diff < 0.0001 | 6.850e-06 | PASS   |
| bala1  | max rel. freq diff < 0.001   | 1.176e-05 | PASS   |
| bala1  | 6 zero modes, 0 imaginary    | 6         | PASS   |
| alaala | |Δ E_total| < 0.0001         | 6.429e-07 | PASS   |
| alaala | gradient RMS < 1.0e-5        | 8.069e-08 | PASS   |
| alaala | |Δ E_min| < 0.0001           | 6.437e-07 | PASS   |
| alaala | mean rel. freq diff < 0.0001 | 6.813e-06 | PASS   |
| alaala | max rel. freq diff < 0.001   | 1.185e-05 | PASS   |
| alaala | 6 zero modes, 0 imaginary    | 6         | PASS   |
