#!/usr/bin/env python3
"""s46: verify the 'Region I' epitope on OUR instrument, before any epitope decision.

WHY THIS EXISTS. A literature survey proposed Region I (Chen et al., Commun Biol 2025,
doi:10.1038/s42003-025-09030-7) as an alternative to our R108 receptor-site epitope, and
supplied residue labels V74/L75/I97/W114/Y115 plus cations K98/K65/K112/H78 in "mature"
numbering. Those labels do NOT land on those residues in our sequence -- but the MOTIF
WINDOWS the same source quoted (VLLTHT, KPWYEPIYL) do. Five numbering schemes are live in
this project and s2 says never infer a position. So this locates Region I by motif and
reports every number in all three of our schemes, computed, not carried over.

positional = mature - 5   (our renumbered trimer is mature 6-157 as 1-152)
canonical  = mature + 76  (UNIPROT_START 77, asserted by species_and_histidines.py T2)

Run: python3 analysis/02-tnf/region1_verify.py
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from species_and_histidines import HUMAN, MOUSE, align, UNIPROT_START  # the validated aligner

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
       "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
       "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
TRIMER = os.path.join(HERE, "..", "..", "targets", "tnf", "tnf_trimer_renum.pdb")

def mature_to_pos(m): return m - 5
def mature_to_canon(m): return m + 76

def human_mouse_map():
    """mature human index (1-based) -> mouse residue letter or None if deleted in mouse."""
    A, B = align(HUMAN, MOUSE)
    out, hm = {}, 0
    for x, y in zip(A, B):
        if x != "-":
            hm += 1
            out[hm] = (x, None if y == "-" else y)
    return out

def load_trimer():
    """chain -> {positional resnum: (resname, [(x,y,z) heavy atoms])}"""
    ch = {}
    with open(TRIMER) as f:
        for l in f:
            if not l.startswith("ATOM"):
                continue
            el = l[76:78].strip() or l[12:16].strip()[0]
            if el == "H":
                continue
            c, rn, rs = l[21], int(l[22:26]), l[17:20].strip()
            d = ch.setdefault(c, {}).setdefault(rn, (rs, []))
            d[1].append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return ch

def dmin(a, b):
    return min(((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ** 0.5
               for x1, y1, z1 in a for x2, y2, z2 in b)

def main():
    hm = human_mouse_map()

    # --- T1: our numbering chain must reproduce the project anchor R108 = mature 32 = pos 27
    assert hm[32][0] == "R", f"T1 FAIL: mature 32 is {hm[32][0]}, not R"
    assert mature_to_pos(32) == 27 and mature_to_canon(32) == 108, "T1 FAIL: numbering chain"
    print("T1  numbering chain OK: R108 canonical = mature 32 = positional 27")

    # --- locate Region I by MOTIF, never by the supplied labels
    motifs = {"VLLTHT": None, "KPWYEPIYL": None}
    for m in motifs:
        i = HUMAN.find(m)
        assert i >= 0 and HUMAN.find(m, i+1) == -1, f"motif {m} absent or not unique"
        motifs[m] = i + 1                      # 1-based mature start
    print(f"T2  motifs unique: VLLTHT at mature {motifs['VLLTHT']}, "
          f"KPWYEPIYL at mature {motifs['KPWYEPIYL']}")

    # the five hotspots Chen et al. name, resolved by their chemistry inside those windows
    want = [("V", "VLLTHT", 0), ("L", "VLLTHT", 1), ("W", "KPWYEPIYL", 2),
            ("Y", "KPWYEPIYL", 3)]
    hotspots = []
    for letter, motif, off in want:
        m = motifs[motif] + off
        assert hm[m][0] == letter, f"expected {letter} at mature {m}, got {hm[m][0]}"
        hotspots.append(m)
    # the isoleucine hotspot: the only I in the second cluster's upstream strand
    ile = [m for m in range(90, 100) if hm[m][0] == "I"]
    assert len(ile) == 1, f"expected one Ile in mature 90-99, found {ile}"
    hotspots.append(ile[0])
    hotspots.sort()

    print("\n=== REGION I HOTSPOTS, located by motif ===")
    print(f"{'mature':>7}{'pos':>6}{'canon':>7}  human  mouse  conserved")
    cons = 0
    for m in hotspots:
        h, mo = hm[m]
        ok = (mo == h)
        cons += ok
        print(f"{m:>7}{mature_to_pos(m):>6}{mature_to_canon(m):>7}    {h}      "
              f"{mo or '-'}     {'YES' if ok else 'no'}")
    print(f"\nconserved human->mouse: {cons}/{len(hotspots)}")

    # --- geometry on OUR corrected trimer
    ch = load_trimer()
    print(f"\n=== GEOMETRY, {os.path.basename(TRIMER)} chains {sorted(ch)} ===")
    c1 = [m for m in hotspots if m < 100]
    c2 = [m for m in hotspots if m >= 100]
    print(f"cluster 1 (mature {c1}) vs cluster 2 (mature {c2})")
    best = None
    for a in c1:
        for b in c2:
            for ca, cb in (("A","B"),("B","C"),("C","A"),("B","A"),("C","B"),("A","C")):
                pa, pb = mature_to_pos(a), mature_to_pos(b)
                if pa in ch[ca] and pb in ch[cb]:
                    d = dmin(ch[ca][pa][1], ch[cb][pb][1])
                    if best is None or d < best[0]:
                        best = (d, a, ca, b, cb)
    print(f"closest inter-protomer hotspot pair: mature {best[1]}({best[2]}) <-> "
          f"{best[3]}({best[4]}) = {best[0]:.2f} A")

    # --- cations near Region I, and whether mouse keeps them
    patch = [(mature_to_pos(m), "A" if m < 100 else "B") for m in hotspots]
    atoms = [a for p, c in patch for a in ch[c][p][1]]
    print("\n=== CATIONS within 12 A of the Region I patch ===")
    print(f"{'mature':>7}{'pos':>6}{'canon':>7}  res  chain  dist_A  mouse  usable")
    found = []
    for c in sorted(ch):
        for p, (rs, at) in sorted(ch[c].items()):
            if rs not in ("ARG", "LYS", "HIS"):
                continue
            d = dmin(at, atoms)
            if d > 12.0:
                continue
            m = p + 5
            h, mo = hm[m]
            # usable for a His-cation mechanism only if mouse keeps a cation there
            usable = mo in ("R", "K", "H")
            found.append((d, m, p, rs, c, mo, usable))
    for d, m, p, rs, c, mo, usable in sorted(found)[:14]:
        print(f"{m:>7}{p:>6}{mature_to_canon(m):>7}  {AA3[rs]}    {c}   {d:6.2f}    "
              f"{mo or '-'}     {'YES' if usable else 'NO'}")
    nmouse = sum(1 for f in found if f[6])
    print(f"\ncations within 12 A: {len(found)} total, {nmouse} retained as a cation in mouse")

    # --- how far is our current anchor?
    r108 = ch["B"][27][1]
    print(f"\n=== OUR CURRENT ANCHOR ===")
    print(f"R108 (positional 27, chain B) to nearest Region I hotspot atom: "
          f"{dmin(r108, atoms):.2f} A  -> {'same site' if dmin(r108, atoms) < 10 else 'DIFFERENT SITE'}")
    print("\nregion1_verify OK")

if __name__ == "__main__":
    main()
