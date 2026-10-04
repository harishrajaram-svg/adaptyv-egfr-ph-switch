#!/usr/bin/env python3
"""Free geometry filter: find designs where ANY carboxylate reaches a target histidine.

THE STRATEGY THIS IMPLEMENTS. Measured 2026-10-04 on our own pools:
    unpinned designs with an Asp/Glu within 4.0 A of H433:  46/240 = 19.2%
    pinned carboxylate, full ECD:                            2/160 =  1.3%
    pinned, tight patch:                                     0/90
    pinned, cropped target:                                  0/60
Pinning fixes a sequence position before the backbone exists, so it fixes which FACE
of the binder the acid occupies -- usually the wrong one. Letting inverse folding place
acids freely and SCREENING afterwards is 15x better. rank-1 rimA01_r15's GLU92 (2.74 A)
and ASP67 were never pinned; the only verified switch in the project came this way.

This screen costs nothing -- a distance measurement on structures that already exist --
and it belongs BEFORE the GPU scoring step, not after. Running the pH gate last is what
cost this project a week.

Usage:
    python3 bin/geometry_filter.py <glob> --his N [--chain-binder A --chain-target B]
    python3 bin/geometry_filter.py 'runs/egfr-crop/cr_*/final_ranked_designs/final_30_designs/rank*.cif' --his 123
"""
import argparse, glob, sys
import gemmi

TIP = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
BAR = 4.0      # conventional salt-bridge / H-bond distance, NOT fitted to our positives


def nearest_acid(path, his_num, cb="A", ct="B"):
    st = gemmi.read_structure(path); st.setup_entities()
    C = {c.name: c for c in st[0]}
    A, B = C.get(cb), C.get(ct)
    if A is None or B is None:
        return None
    h = [r for r in B if r.seqid.num == his_num and r.name == "HIS"]
    if not h:
        return None
    ring = [a.pos for a in h[0] if a.name in ("ND1", "NE2")]
    if not ring:
        return None
    best, who = 999.0, None
    for r in A:
        if r.name not in TIP:
            continue
        for a in r:
            if a.name in TIP[r.name]:
                for rn in ring:
                    d = a.pos.dist(rn)
                    if d < best:
                        best, who = d, f"{r.name}{r.seqid.num}"
    return best, who


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pattern")
    ap.add_argument("--his", type=int, required=True,
                    help="target histidine residue number AS NUMBERED IN THESE FILES "
                         "(123 for the domain-III crop output, 99 for the d3 complex "
                         "PDBs, 409 for full-ECD BoltzGen output -- they differ)")
    ap.add_argument("--chain-binder", default="A")
    ap.add_argument("--chain-target", default="B")
    ap.add_argument("--bar", type=float, default=BAR)
    a = ap.parse_args()

    files = sorted(glob.glob(a.pattern))
    if not files:
        print(f"no files matched {a.pattern}", file=sys.stderr); return 2
    rows, skipped = [], 0
    for f in files:
        r = nearest_acid(f, a.his, a.chain_binder, a.chain_target)
        if r is None:
            skipped += 1; continue
        rows.append((r[0], r[1], f))
    rows.sort()
    hits = [r for r in rows if r[0] < a.bar]
    for d, who, f in hits:
        print(f"{d:6.2f} A  {who:>9}  {f.split('/')[-1]}")
    print(f"\n{len(hits)}/{len(rows)} within {a.bar} A"
          f" = {100*len(hits)/max(len(rows),1):.1f}%"
          + (f"   ({skipped} skipped: no HIS {a.his} in chain {a.chain_target})" if skipped else ""))
    print(f"reference rates: unpinned 19.2% | pinned full-ECD 1.3% | pinned cropped 0%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
