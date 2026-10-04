#!/usr/bin/env python3
"""Which target residues does the binder actually touch? Declared epitope vs delivered.

Window A measured on 2026-10-03 that BoltzGen's `binding_types` is a WEAK HINT, not a
constraint: across 7 arms only 4-12 of 20 designs per arm touched the declared epitope and
roughly half bound domain I instead. Narrowing the declared patch from 2 to 7 residues did
not help. That makes "did the binder go where we asked" a measurement worth having rather
than an assumption, for any generator.

This reads complexes written by modal_mosaic.py step 4: binder = chain A, target = chain B
in POSITIONAL numbering (egfr_ecd_6aru_renum.pdb, contiguous 1..609, positional =
canonical - 27, so H433 = 406). It does NOT apply BoltzGen's +3 output shift -- pass
--offset 3 for files that carry it.

usage:
  python3 bin/epitope_contact.py runs/mosaic/tune01/complex/*.pdb
  python3 bin/epitope_contact.py --epitope 403-409 --anchor 406 runs/mosaic/*/complex/*.pdb
"""
import argparse
import glob
import sys

import gemmi

# EGFR ectodomain boundaries in POSITIONAL numbering (canonical - 27).
DOMAINS = (("I", 1, 138), ("II", 139, 283), ("III", 284, 453), ("IV", 454, 609))


def parse_residues(spec):
    out = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = (int(v) for v in part.split("-"))
            out.extend(range(lo, hi + 1))
        elif part:
            out.append(int(part))
    return sorted(set(out))


def domain_of(resnum):
    for name, lo, hi in DOMAINS:
        if lo <= resnum <= hi:
            return name
    return "?"


def contacts(path, cutoff, offset):
    st = gemmi.read_structure(path)
    st.remove_hydrogens()
    model = st[0]
    if len(model) < 2:
        return None
    binder, target = model[0], model[1]
    # heavy-atom distance, every binder atom against every target atom. 76 x 170 residues
    # is small enough that a NeighborSearch buys nothing and costs clarity.
    hit = {}
    for tres in target:
        num = tres.seqid.num - offset
        best = 1e9
        for tatom in tres:
            for bres in binder:
                for batom in bres:
                    d = tatom.pos.dist(batom.pos)
                    if d < best:
                        best = d
        if best <= cutoff:
            hit[num] = round(best, 2)
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdbs", nargs="+")
    ap.add_argument("--epitope", default="403-409",
                    help="declared target residues, positional (default the H433 patch)")
    ap.add_argument("--anchor", type=int, default=406, help="H433 in positional numbering")
    ap.add_argument("--cutoff", type=float, default=5.0, help="heavy-atom contact cutoff A")
    ap.add_argument("--offset", type=int, default=0,
                    help="subtract this from target resnums (use 3 for BoltzGen output)")
    a = ap.parse_args()

    declared = set(parse_residues(a.epitope))
    paths = [p for pat in a.pdbs for p in sorted(glob.glob(pat))]
    if not paths:
        sys.exit("no PDBs matched")

    print(f"declared epitope: {sorted(declared)}  anchor {a.anchor}  "
          f"cutoff {a.cutoff} A  offset {a.offset}")
    print(f"{'design':34} {'n_contact':>9} {'in_patch':>8} {'frac':>6} "
          f"{'anchor':>6} {'domains':>14}  nearest_patch_residue")
    n_on_patch = 0
    for path in paths:
        hit = contacts(path, a.cutoff, a.offset)
        name = path.split("/")[-1].replace(".pdb", "")
        if hit is None:
            print(f"{name:34} -- single-chain file, skipped")
            continue
        in_patch = sorted(declared & set(hit))
        doms = sorted({domain_of(r) for r in hit})
        frac = len(in_patch) / len(hit) if hit else 0.0
        anchor_hit = "YES" if a.anchor in hit else "no"
        nearest = min(((hit[r], r) for r in in_patch), default=(None, None))
        if in_patch:
            n_on_patch += 1
        print(f"{name:34} {len(hit):9d} {len(in_patch):8d} {frac:6.2f} {anchor_hit:>6} "
              f"{','.join(doms):>14}  "
              + (f"{nearest[1]} at {nearest[0]} A" if nearest[1] else "none"))
    print(f"\n{n_on_patch}/{len(paths)} designs touch the declared patch at all")
    print("NOTE: a domain-III crop makes domain I/II/IV contact impossible by construction "
          "-- that escape route is closed by the crop, not by the objective.")


if __name__ == "__main__":
    main()
