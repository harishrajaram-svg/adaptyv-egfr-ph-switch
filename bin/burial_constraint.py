#!/usr/bin/env python3
"""Burial-on-binding, as a per-position criterion for keeping acids out of the wrong places.

WHY THE EXISTING CRITERION IS WRONG. bin/mpnn_constraints.py omits ASP and GLU at binder
positions within IFACE_CA = 8.0 A of any target atom -- PROXIMITY. Running it on 48 designs
(METHODS s6.10a) moved the right-direction rate of the pH linkage from 37% to 38%: nothing. The
split by residue type said why. Acids are wrong 96% of the time in both arms, but 69% of all
moving sites have NO COUNTER-CHARGE within 6 A, which is the desolvation signature -- a residue
whose pKa shifts because binding took its water away, not because it met a charge. Proximity is
the wrong proxy for that. A residue on the interface RIM is proximate and still solvated; one in
a pocket is buried. Those are different sets.

WHAT THIS MEASURES INSTEAD. For each binder CA, the number of TARGET heavy atoms within
NEIGHBOUR_R. That count is what binding adds to the residue's environment, so it is burial-on-
binding directly rather than a proxy for it. It needs only the CA-only Genie 3 backbone and the
target, both available before any sequence exists.

NEIGHBOUR_R = 10.0 A is the conventional neighbour-count shell for burial in this literature and
is NOT fitted here.

\U0001F534 THE THRESHOLD IS NOT INVENTED AND IT IS NOT BLIND. It is set by asking whether this
measure actually separates the positions that showed desolvation shifts in the 48 from those that
did not. That validation runs in --validate and its result is printed, not assumed. If the measure
does not separate them, it is no better than proximity and this file should not be used.

    burial_constraint.py --selftest
    burial_constraint.py --validate      # against the 48 designs' measured linkage
"""
import glob
import json
import math
import os
import sys

NEIGHBOUR_R = 10.0          # A, conventional burial shell; not fitted
BINDER_CHAIN = "A"          # in the Genie 3 complex output


def atoms(path):
    out = []
    for l in open(path):
        if l.startswith("ATOM"):
            out.append((l[21], int(l[22:26]), l[12:16].strip(),
                        (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out


def target_neighbour_count(complex_pdb, binder_chain=BINDER_CHAIN, r=NEIGHBOUR_R):
    """{binder_resnum: number of TARGET heavy atoms within r of that residue's CA}."""
    at = atoms(complex_pdb)
    binder = [(n, p) for c, n, an, p in at if c == binder_chain and an == "CA"]
    tgt = [p for c, n, an, p in at if c != binder_chain]
    if not binder:
        sys.exit(f"REFUSE: no binder CA on chain {binder_chain} in {complex_pdb}")
    if not tgt:
        sys.exit(f"REFUSE: no target atoms outside chain {binder_chain}")
    r2 = r * r
    out = {}
    for n, p in binder:
        out[n] = sum(1 for q in tgt
                     if (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2 <= r2)
    return out


def validate():
    """Does burial separate the measured desolvation movers from everything else?

    The desolvation set is defined by the LINKAGE measurement, not by this file: a moving binder
    site with no counter-charge within 6 A. If burial is the right criterion those positions
    should carry higher target-neighbour counts than the rest of the chain.
    """
    import statistics

    link = json.load(open("/tmp/A2_linkage.json"))
    desolv, other = [], []
    n_designs = 0
    for r in link:
        if "error" in r:
            continue
        stem = r["file"].split("_seed42")[0]
        bb = stem.rsplit("_s", 1)[0]
        cx = f"runs/genie3-out/probe50/tnfa_corrected/pdbs/{bb}.pdb"
        if not os.path.exists(cx):
            continue
        nb = target_neighbour_count(cx)
        n_designs += 1
        hit = {v["resnum"] for v in r.get("contributing", {}).values()
               if v["partner"] == "binder" and v.get("no_partner")}
        for num, cnt in nb.items():
            (desolv if num in hit else other).append(cnt)
    if not desolv:
        print("no desolvation movers found -- cannot validate")
        return 1
    md, mo = statistics.median(desolv), statistics.median(other)
    # tie-corrected Mann-Whitney AUC: P(a desolvation position outranks a random one)
    gt = sum(1 for a in desolv for b in other if a > b)
    eq = sum(1 for a in desolv for b in other if a == b)
    auc = (gt + 0.5 * eq) / (len(desolv) * len(other))
    print(f"validation over {n_designs} designs")
    print(f"  desolvation movers (no counter-charge) n={len(desolv):>5}  "
          f"median target-neighbour count {md:.0f}")
    print(f"  every other binder position          n={len(other):>5}  "
          f"median target-neighbour count {mo:.0f}")
    print(f"  AUC, burial separating the two: {auc:.3f}   (0.5 = no separation)")
    verdict = ("USABLE: burial separates them" if auc >= 0.65 else
               "WEAK: better than chance but thin" if auc >= 0.58 else
               "NOT USABLE: no better than proximity; do not deploy this file")
    print(f"  VERDICT: {verdict}")
    return 0 if auc >= 0.58 else 1


def selftest():
    import tempfile

    def line(ch, n, an, x):
        return (f"ATOM  {n:5d}  {an:<3s} GLY {ch}{n:4d}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00\n")
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        # binder CAs at x = 0 (enclosed), 9 (edge), 40 (far)
        for i, x in enumerate([0.0, 9.0, 40.0], start=1):
            fh.write(line("A", i, "CA", x))
        # a wall of 5 target atoms packed around x=0
        for j, x in enumerate([-3.0, -1.0, 1.0, 3.0, 5.0], start=100):
            fh.write(line("B", j, "CA", x))
        tmp = fh.name
    nb = target_neighbour_count(tmp)
    assert nb[1] == 5, nb        # fully enclosed: all five within 10 A
    assert nb[2] == 4, nb        # rim at x=9: four atoms within 10 A (at 4, 6, 8, 10)
    assert nb[3] == 0, nb        # far
    # MUTATION 1, AND IT IS A WARNING ABOUT THIS MEASURE. A 10 A shell is generous: the rim
    # position scores 4 against the enclosed position's 5, which is a thin margin. In a real
    # 3-D complex enclosure separates better than it does on a line, but the honest arbiter is
    # --validate against the measured desolvation set, not this toy. Do not read the toy as
    # evidence the criterion discriminates.
    nearest = min(abs(9.0 - x) for x in (-3.0, -1.0, 1.0, 3.0, 5.0))
    assert nearest <= 8.0, nearest
    assert nb[2] < nb[1], "a rim position must score below an enclosed one"
    assert nb[2] / nb[1] > 0.5, "and on a line the separation is weak -- see --validate"
    # MUTATION 2: the radius must be load-bearing
    assert target_neighbour_count(tmp, r=2.5)[1] < nb[1]
    # MUTATION 3: refuse rather than guess when a side is missing
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        fh.write(line("A", 1, "CA", 0.0)); only_binder = fh.name
    try:
        target_neighbour_count(only_binder); raise AssertionError("no-target must refuse")
    except SystemExit as e:
        assert "no target atoms" in str(e)
    os.unlink(tmp); os.unlink(only_binder)
    print(f"  ok  target-neighbour count at R={NEIGHBOUR_R} A: enclosed 5, rim 4, far 0")
    print(f"  ok  MUTATION: on a line the rim/enclosed margin is thin (4 vs 5) -- the toy is")
    print(f"      NOT evidence of discrimination; --validate is the arbiter")
    print(f"  ok  MUTATION: the radius is load-bearing")
    print(f"  ok  MUTATION: a complex with no target chain is refused, not guessed")
    print("\nself-tests passed: 4")
    return 0


if __name__ == "__main__":
    a = sys.argv
    if "--validate" in a:
        sys.exit(validate())
    sys.exit(selftest())
