# mutate.py   usage: python mutate.py 5TGZ.pdb 5TGZ_F238L.pdb
import sys
inp, out = sys.argv[1], sys.argv[2]
keep = {"N", "CA", "C", "O", "CB"}          # Leu atoms that share positions with Phe
seen = False
with open(out, "w") as o:
    for l in open(inp):
        if l.startswith(("ATOM", "ANISOU")) and l[21] == "A" and int(l[22:26]) == 237:
            assert l[17:20] == "PHE", "residue 237 is not PHE, check numbering"
            if l[12:16].strip() not in keep:
                continue                     # drop CG, CD1, CD2, CE1, CE2, CZ
            if l.startswith("ANISOU"):
                continue
            l = l[:17] + "LEU" + l[20:]
            seen = True
        o.write(l)
assert seen, "residue 237 not found"