#!/usr/bin/env python3
"""Write a PDB whose B-factor column holds the per-residue relative dMSF (mutant vs WT)
so VMD can colour the structure by it.
usage: python code/dmsf_to_pdb.py WT/WT_min.pdb cb1/final/cb1_dmsf_wt_vs_mut.tsv out.pdb"""
import sys
pdb, tsv, out = sys.argv[1:4]
val = {}
for l in open(tsv, encoding="utf-8"):
    f = l.rstrip("\n").split("\t")
    try: val[int(float(f[0]))] = float(f[4])
    except (ValueError, IndexError): continue            # header / blank lines
n = miss = 0
with open(out, "w") as o:
    for l in open(pdb):
        if l.startswith(("ATOM", "HETATM")):
            r = int(l[22:26]); v = val.get(r)
            if v is None: miss += 1; v = 0.0
            l = l[:60] + "%6.3f" % max(min(v, 99.0), -9.0) + l[66:]
            n += 1
        o.write(l)
print("residues with dMSF: %d, atoms written: %d, atoms without a value (set to 0): %d" % (len(val), n, miss))
