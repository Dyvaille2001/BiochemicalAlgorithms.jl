# mmtk_compare.py  (Python 2, MMTK)
#
# Usage:   python mmtk_compare.py <PDB file or MMTK database name> <tag>
# Examples:
#   python mmtk_compare.py bala1 bala1          (blocked alanine from the MMTK database)
#   python mmtk_compare.py AlaAla.pdb alaala    (Ala-Ala zwitterion from your own PDB)
#
# Writes (all from the SAME run; copy them together to the Julia folder):
#   <tag>_min.pdb               minimized structure (the Julia input)
#   mmtk_<tag>_freqs.txt        non-trivial frequencies in 1/cm, one per line
#   mmtk_<tag>_gradients.txt    per-atom gradient (kJ/mol/A) at the re-read geometry
#   mmtk_<tag>_energies.txt     energy terms at the re-read geometry + the minimum energy
#
# MMTK units: frequency/Units.invcm -> 1/cm (MMTK frequencies are in THz).

import sys
from MMTK import *
from MMTK import Units
from MMTK.Proteins import Protein
from MMTK.ForceFields import Amber94ForceField
from MMTK.NormalModes import VibrationalModes
from MMTK.Minimization import ConjugateGradientMinimizer

if len(sys.argv) < 3:
    print __doc__
    sys.exit(1)
spec, tag = sys.argv[1], sys.argv[2]

universe = InfiniteUniverse(Amber94ForceField())
universe.protein = Protein(spec)
natoms = universe.numberOfAtoms()
print 'atoms:', natoms, '  (non-trivial modes expected: %d)' % (3 * natoms - 6)

minimizer = ConjugateGradientMinimizer(universe)
minimizer(convergence=1.e-5, steps=50000)
e_min = universe.energy()
print 'minimized energy (kJ/mol):', e_min

universe.protein.writeToFile(tag + '_min.pdb')


def print_terms(u, label):
    print label
    terms = u.energyTerms()
    for k in sorted(terms.keys()):
        print '   %-30s %14.6f' % (k, terms[k])
    print '   %-30s %14.6f' % ('TOTAL', u.energy())

print_terms(universe, 'MMTK energy terms at the minimum (kJ/mol):')

# Same geometry Julia will read: PDB coordinates have only 3 decimals.
u2 = InfiniteUniverse(Amber94ForceField())
u2.protein = Protein(tag + '_min.pdb')
print_terms(u2, 'MMTK energy terms after re-reading %s_min.pdb (compare THIS with Julia):' % tag)

ef = open('mmtk_%s_energies.txt' % tag, 'w')
ef.write('minimum_energy %.9f\n' % e_min)
terms = u2.energyTerms()
for k in sorted(terms.keys()):
    ef.write('reread %s %.9f\n' % (k.replace(' ', '_'), terms[k]))
ef.write('reread TOTAL %.9f\n' % u2.energy())
ef.close()

e2, g2 = u2.energyAndGradients()
gf = open('mmtk_%s_gradients.txt' % tag, 'w')
for atom in u2.atomList():
    p = atom.position()
    g = g2[atom]
    gf.write('%s %.6f %.6f %.6f %.8f %.8f %.8f\n' % (
        atom.fullName().replace(' ', '_'),
        p[0] * 10., p[1] * 10., p[2] * 10., g[0] / 10., g[1] / 10., g[2] / 10.))
gf.close()

modes = VibrationalModes(universe)
out = open('mmtk_%s_freqs.txt' % tag, 'w')
for mode in modes[6:]:                       # the first 6 are rigid-body modes
    out.write('%.6f\n' % (mode.frequency / Units.invcm))
out.close()
print 'check: print of one mode ->', modes[6]
print 'check: frequency/invcm   ->', modes[6].frequency / Units.invcm
print 'wrote: %s_min.pdb  mmtk_%s_freqs.txt  mmtk_%s_gradients.txt  mmtk_%s_energies.txt' % (
    tag, tag, tag, tag)
