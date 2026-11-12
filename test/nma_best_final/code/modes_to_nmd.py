#!/usr/bin/env python3
"""ProDy ANM (constant spring, 15 A; the model whose eigenvalues/modes your Julia code reproduced
to ~1e-13) -> NMD file for VMD's Normal Mode Wizard.
usage: python code/modes_to_nmd.py 1aki/1aki.pdb 1aki/1aki_anm.nmd [n_modes=10]"""
import sys
from prody import parsePDB, ANM, writeNMD
pdb, out = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 10
ca = parsePDB(pdb).select("calpha")
anm = ANM("anm"); anm.buildHessian(ca, cutoff=15.0, gamma=1.0); anm.calcModes(n_modes=n)
writeNMD(out, anm, ca)
print("wrote", out, "with", n, "modes for", ca.numAtoms(), "Calpha atoms")
