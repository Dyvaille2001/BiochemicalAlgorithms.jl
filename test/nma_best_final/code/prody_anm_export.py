# prody_anm_export.py  (ProDy; Python 2 or 3)
#
# Usage:   python code/prody_anm_export.py 1aki/1aki.pdb [n_modes] [output_dir]
#
# Exports the first n_modes (default 20) NON-ZERO modes and eigenvalues of two
# Cα elastic network models, built on the SAME Cα atoms:
#
#   const : constant spring k_ij = 1, cutoff 15 A   (ProDy defaults)
#   dw    : distance-weighted k_ij = 1/r^2, cutoff 12 A   (the thesis DWENM)
#
# Files written (3N lines x n_modes columns; x1 y1 z1 x2 ... ordering):
#   prody_const_modes.txt  prody_const_eigvals.txt
#   prody_dw_modes.txt     prody_dw_eigvals.txt
#
# The distance-weighted model uses ProDy's Gamma class, gamma(self, dist2, i, j).
# This is the signature I expect; check `help(prody.Gamma)` if the DW export fails.
import sys
import os
import numpy as np
from prody import parsePDB, ANM, Gamma, confProDy

confProDy(verbosity='none')

pdb_in = sys.argv[1] if len(sys.argv) > 1 else '1aki.pdb'
nmodes = int(sys.argv[2]) if len(sys.argv) > 2 else 20
outdir = (sys.argv[3] if len(sys.argv) > 3 else
          os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'prody')))
os.makedirs(outdir, exist_ok=True)

prot = parsePDB(pdb_in)
ca = prot.select('name CA and protein')
if ca is None:
    sys.exit('no C-alpha atoms selected')
n = ca.numAtoms()
print('C-alpha atoms: %d' % n)


def export(anm, tag):
    anm.calcModes(n_modes=nmodes)            # zero modes are dropped by default
    V = np.asarray(anm.getEigvecs(), dtype=float)
    if V.shape[0] != 3 * n:
        sys.exit('unexpected eigenvector shape %s' % (V.shape,))
    V = V / np.linalg.norm(V, axis=0)        # explicit column normalisation
    np.savetxt(os.path.join(outdir, 'prody_%s_modes.txt' % tag), V,
               fmt='%.10e', delimiter=' ')
    np.savetxt(os.path.join(outdir, 'prody_%s_eigvals.txt' % tag), np.asarray(anm.getEigvals()),
               fmt='%.12e')
    print('%-5s eigenvalues 1-5: %s' % (tag, anm.getEigvals()[:5]))


const = ANM('const')
const.buildHessian(ca, cutoff=15.0, gamma=1.0)
export(const, 'const')


class DW(Gamma):
    """k_ij = 1 / r_ij^2  (dist2 = r^2)."""
    def gamma(self, dist2, i, j):
        return 1.0 / dist2


dw = ANM('dw')
dw.buildHessian(ca, cutoff=12.0, gamma=DW())
export(dw, 'dw')
print('done')
