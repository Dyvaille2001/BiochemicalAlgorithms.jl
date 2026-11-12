#!/usr/bin/env python3
"""Build explicit Leu rotamers at one residue of a PDB, without PyMOL.

usage: python make_leu_rotamers.py 5TGZ.pdb PREFIX [residue=237] [chain=A]

Writes PREFIX_<tag>.pdb for the best clash-free rotamer of each chi1 class
(g-, t, g+), i.e. up to three structurally distinct F238L models, and prints a
table of all nine (chi1, chi2) candidates. CB, N, CA, C, O of the residue stay
exactly as in the input; CG, CD1, CD2 are placed with ideal Leu geometry.
CD2 = CD1 + 120 deg: the nine (chi1, H-slot) combinations cover every distinct
staggered conformer, whichever way CD1/CD2 are labelled.
"""
import sys, math
import numpy as np

CHI = {"g-": -65.0, "t": 175.0, "g+": 60.0}

def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.degrees(math.atan2(y, x))

def place(a, b, c, bond, angle, dih):
    """atom d with |c-d| = bond, angle(b,c,d) = angle, dihedral(a,b,c,d) = dih (degrees)"""
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc); n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    th, ph = math.radians(angle), math.radians(dih)
    d2 = np.array([-bond * math.cos(th), bond * math.sin(th) * math.cos(ph), bond * math.sin(th) * math.sin(ph)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n

def read_atoms(path):
    out = []
    for l in open(path):
        if l.startswith("ATOM"):
            out.append((l, l[12:16].strip(), l[17:20], l[21], int(l[22:26]),
                        np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out

def candidates(N, CA, CB):
    cands = []
    for t1, chi1 in CHI.items():
        CG = place(N, CA, CB, 1.53, 116.0, chi1)
        for t2, chi2 in CHI.items():
            CD1 = place(CA, CB, CG, 1.52, 110.5, chi2)
            CD2 = place(CA, CB, CG, 1.52, 110.5, chi2 + 120.0)
            cands.append((t1, t2, chi1, chi2, CG, CD1, CD2))
    return cands

def main():
    if len(sys.argv) < 3: sys.exit(__doc__)
    inp, prefix = sys.argv[1], sys.argv[2]
    resi = int(sys.argv[3]) if len(sys.argv) > 3 else 237
    chain = sys.argv[4] if len(sys.argv) > 4 else "A"
    atoms = read_atoms(inp)
    own = {a[1]: a for a in atoms if a[3] == chain and a[4] == resi}
    assert own, "residue %s %d not found" % (chain, resi)
    for k in ("N", "CA", "CB"): assert k in own, "missing " + k
    N, CA, CB = own["N"][5], own["CA"][5], own["CB"][5]
    others = np.array([a[5] for a in atoms if not (a[3] == chain and a[4] == resi)])
    rows = []
    for t1, t2, c1, c2, CG, CD1, CD2 in candidates(N, CA, CB):
        # self-test of the placement
        assert abs(dihedral(N, CA, CB, CG) - c1) < 0.1 or abs(abs(dihedral(N, CA, CB, CG) - c1) - 360) < 0.1
        sc = np.array([CG, CD1, CD2])
        d = np.linalg.norm(sc[:, None, :] - others[None, :, :], axis=2)
        rows.append(dict(t1=t1, t2=t2, c1=c1, c2=c2, atoms=(CG, CD1, CD2),
                         dmin=d.min(), n_close=int((d < 4.0).sum())))
    print("%-4s %-4s %7s %7s %12s %14s" % ("chi1", "chi2", "chi1°", "chi2°", "min dist (Å)", "contacts<4 Å"))
    for r in rows:
        print("%-4s %-4s %7.1f %7.1f %12.2f %14d" % (r["t1"], r["t2"], r["c1"], r["c2"], r["dmin"], r["n_close"]))
    chosen = []
    for t1 in CHI:
        ok = [r for r in rows if r["t1"] == t1 and r["dmin"] >= 2.9]
        if ok: chosen.append(max(ok, key=lambda r: r["dmin"]))      # least strained of that chi1 class
        else:  print("no clash-free rotamer for chi1 =", t1)
    if not chosen: sys.exit("no usable rotamer")
    bfac = float(own["CB"][0][60:66])
    for k, r in enumerate(chosen):
        tag = "rot%s" % "ABC"[k]
        out = []
        for l, name, rn, ch, ri, xyz in atoms:
            if ch == chain and ri == resi:
                if name not in ("N", "CA", "C", "O", "CB"): continue
                l = l[:17] + "LEU" + l[20:]
            out.append(l)
            if ch == chain and ri == resi and name == "CB":
                for nm, xyz2, el in zip(("CG", "CD1", "CD2"), r["atoms"], "CCC"):
                    out.append("ATOM  %5d  %-3s LEU %s%4d    %8.3f%8.3f%8.3f  1.00%6.2f           %s  \n"
                               % (0, nm, chain, resi, xyz2[0], xyz2[1], xyz2[2], bfac, el))
        n = 0; fin = []
        for l in out:
            n += 1; fin.append(l[:6] + "%5d" % n + l[11:])
        path = "%s_%s.pdb" % (prefix, tag)
        open(path, "w").writelines(fin + ["END\n"])
        print("wrote %s   (chi1 %s %.0f°, chi2 %.0f°, min dist %.2f Å)" % (path, r["t1"], r["c1"], r["c2"], r["dmin"]))

if __name__ == "__main__":
    main()
