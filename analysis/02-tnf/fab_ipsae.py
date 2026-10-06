#!/usr/bin/env python3
"""ipSAE for a MULTI-CHAIN binder against a multi-chain target, with the guard §30 lacked.

WHY THIS EXISTS. `bin/ipsae_min.py` refuses any complex with more than one inter-chain pair,
on purpose: ipSAE_min is a min over two ALIGNMENT DIRECTIONS of ONE interface, and letting
min() roam across chain pairs reports an intra-binder VH/VL packing score as a binding score.
A Fab against the TNF trimer has 5 chains and 10 pairs, so it has no defined answer there.
The three aggregation rules below are the ones PRE-REGISTERED in challenge §29, committed
before any number was computed; §30 reports what they said. Do not add a fourth rule after
seeing results -- that is the move §13 exists to prevent.

THE GUARD. §30's G4 leg read exactly 0.0000 on all nine inter-chain pairs and was nearly
recorded as a spectacular selectivity result. It was a dead fold: the TARGET's own protomers
(A:B, A:C, B:C) also read 0.0000, where human TNF reads 0.71-0.74, because the construct
carried LT-alpha's signal peptide -- the exact error §27 had already diagnosed and fixed, in a
file this run never pointed at. A target that has not assembled cannot report a meaningful
interface, so `--min-target-packing` asserts the target's internal packing before any
binder:target number is believed. It FAILS CLOSED: no score, not a zero.

usage:
  python3 analysis/02-tnf/fab_ipsae.py --dir runs/gate-g-fab/g-fab \
      --target A,B,C --binder D,E
  python3 analysis/02-tnf/fab_ipsae.py --selftest
"""
import glob
import itertools
import json
import math
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RULES = ("max_pair", "sum_pairs", "top_contact")
MIN_TARGET_PACKING = 0.10


def _ipsae_paths():
    """Same resolution bin/ipsae_min.py uses, and for the same clone-safety reasons."""
    vend = ROOT / "outbox" / "ipsae-fixtures" / "vendor" / "ipsae" / "ipsae.py"
    ips = ROOT / "ipsae" / "ipsae.py"
    if not ips.exists() and vend.exists():
        ips = vend
    venv = ROOT / ".venv" / "bin" / "python"
    return (venv if venv.exists() else Path(sys.executable)), ips


def pair_minima(pae, struct, pae_cut=10, dist_cut=10):
    """{(c1,c2): min over the two alignment directions} for every inter-chain pair.

    Only pairs carrying BOTH finite reciprocal directions are returned -- a single
    direction reported as if it were a minimum is the bug bin/ipsae_min.py guards (4).
    """
    pae, struct = Path(pae).resolve(), Path(struct).resolve()
    out = struct.with_name(f"{struct.stem}_{pae_cut}_{dist_cut}.txt")
    if not out.exists():
        py, ips = _ipsae_paths()
        r = subprocess.run([str(py), str(ips), str(pae), str(struct),
                            str(pae_cut), str(dist_cut)],
                           capture_output=True, text=True, cwd=struct.parent)
        if r.returncode != 0 or not out.exists():
            return None
    asym = defaultdict(list)
    for line in out.read_text().splitlines():
        row = line.split()
        if line.startswith("Chn1") or len(row) <= 5 or row[4] != "asym":
            continue
        try:
            v = float(row[5])
        except ValueError:
            continue
        if math.isfinite(v):
            asym[tuple(sorted((row[0], row[1])))].append(v)
    return {k: min(v) for k, v in asym.items() if len(v) == 2}


def check_target_packing(pairs, target, floor=MIN_TARGET_PACKING):
    """Return a complaint string if the TARGET did not assemble, else ''. Fails closed."""
    internal = {k: v for k, v in pairs.items() if k[0] in target and k[1] in target}
    if len(target) > 1 and not internal:
        return f"target chains {','.join(target)} report NO mutual pairs -- cannot trust any interface"
    dead = {":".join(k): v for k, v in internal.items() if v < floor}
    if dead:
        return (f"target did not assemble: internal packing {dead} below {floor} "
                f"-- a target that has not folded cannot report a binder interface "
                f"(see challenge §30, the voided LT-alpha G4 leg)")
    return ""


def contacts(struct, a, b, cutoff=5.0):
    """Heavy-atom contact count between two chains. STRUCTURE, not score -- so the
    top_contact rule's pair selection cannot see the answer it is about to report."""
    import gemmi
    st = gemmi.read_structure(str(struct))
    st.setup_entities()
    st.remove_ligands_and_waters()
    ns = gemmi.NeighborSearch(st, cutoff).populate()
    n = 0
    for res in st[0][a]:
        for atom in res:
            if any(m.to_cra(st[0]).chain.name == b
                   for m in ns.find_atoms(atom.pos, "\0", radius=cutoff)):
                n += 1
    return n


def aggregate(pairs, target, binder, struct=None, pick=None):
    """The three §29 rules. Returns (dict_of_rule_values, chosen_pair) or (None, complaint)."""
    bad = check_target_packing(pairs, target)
    if bad:
        return None, bad
    bt = {k: v for k, v in pairs.items() if (k[0] in binder) != (k[1] in binder)}
    want = len(target) * len(binder)
    if len(bt) != want:
        return None, f"{len(bt)} binder:target pairs with two directions, expected {want}"
    if pick is None and struct is not None:
        pick = tuple(sorted(max(((b, t) for b in binder for t in target),
                                key=lambda p: contacts(struct, *p))))
    vals = {"max_pair": max(bt.values()), "sum_pairs": sum(bt.values())}
    if pick is not None:
        vals["top_contact"] = bt[pick]
    return vals, pick


def rank(vals):
    order = sorted(range(len(vals)), key=lambda i: -vals[i])
    r = [0] * len(vals)
    for pos, i in enumerate(order):
        r[i] = pos + 1
    return r


def spearman(a, b):
    ra, rb, n = rank(a), rank(b), len(a)
    return 1 - 6 * sum((x - y) ** 2 for x, y in zip(ra, rb)) / (n * (n * n - 1))


def exact_p(a, b):
    """One-sided permutation p on small n -- no scipy, and exact beats asymptotic at n=4."""
    obs = spearman(a, b)
    hits = tot = 0
    for perm in itertools.permutations(b):
        tot += 1
        hits += spearman(a, list(perm)) >= obs - 1e-12
    return hits / tot


def selftest():
    T, B = ("A", "B", "C"), ("D", "E")
    packed = {("A", "B"): 0.73, ("A", "C"): 0.71, ("B", "C"): 0.72}
    bt = {(t, b): v for (t, b), v in zip(
        [(t, b) for t in T for b in B], [0.26, 0.19, 0.45, 0.35, 0.46, 0.43])}

    # 1. the guard catches the G4 failure: EVERY pair dead, including the target's own
    dead = {k: 0.0 for k in list(packed) + list(bt)}
    vals, why = aggregate(dead, T, B, pick=("C", "E"))
    assert vals is None and "did not assemble" in why, why

    # 2. MUTATION TEST -- the guard must not be satisfiable by the binder pairs alone.
    #    A target that failed to fold while the binder:target numbers look fine is the
    #    case that nearly got recorded as selectivity.
    sneaky = dict(bt)
    sneaky.update({k: 0.0 for k in packed})
    vals, why = aggregate(sneaky, T, B, pick=("C", "E"))
    assert vals is None, "guard passed on an unassembled target -- it is blind"

    # 3. a properly packed target scores, and the rules compute as §29 defines them
    good = dict(packed)
    good.update(bt)
    vals, pick = aggregate(good, T, B, pick=("C", "E"))
    assert vals is not None, pick
    assert abs(vals["max_pair"] - 0.46) < 1e-9, vals
    assert abs(vals["sum_pairs"] - sum(bt.values())) < 1e-9, vals
    assert abs(vals["top_contact"] - good[("C", "E")]) < 1e-9, vals
    # max_pair must be a MAX over binder:target only -- never the 0.73 target packing
    assert vals["max_pair"] < 0.5, "max_pair leaked an intra-target pair"

    # 4. a missing pair fails closed rather than scoring on five of six
    short = dict(good)
    del short[("C", "E")]
    vals, why = aggregate(short, T, B, pick=("C", "E"))
    assert vals is None and "expected 6" in why, why

    # 5. the statistics, against hand-checked values from §30
    kd = [-4.6, -46.3, -77.3, -112.0]
    mx = [0.5757, 0.4963, 0.4197, 0.5405]
    assert abs(spearman(mx, kd) - 0.4) < 1e-9, spearman(mx, kd)
    assert abs(spearman([0.4412, 0.3794, 0.5405, 0.4939],
                        [-46.3, -77.3, -112.0, -4.6]) - (-0.2)) < 1e-9
    assert abs(spearman(kd, kd) - 1.0) < 1e-9
    assert abs(exact_p(mx, kd) - 0.375) < 1e-9, exact_p(mx, kd)
    assert exact_p(kd, kd) < 0.05

    print("guard         all-dead complex rejected; target-only zeros rejected (mutation);")
    print("              missing pair fails closed; max_pair cannot leak target packing")
    print(f"stats         rho(max_pair) = {spearman(mx, kd):.3f}  p = {exact_p(mx, kd):.3f}  "
          f"rho(top_contact) = -0.200  rho(self) = 1.000")
    print("selftest OK")


def main():
    if "--selftest" in sys.argv:
        return selftest()
    args = sys.argv[1:]
    g = lambda f, d: args[args.index(f) + 1] if f in args else d
    root = g("--dir", None)
    if not root:
        sys.exit(__doc__)
    target = tuple(g("--target", "A,B,C").split(","))
    binder = tuple(g("--binder", "D,E").split(","))

    by, pick, void = defaultdict(lambda: defaultdict(list)), {}, {}
    for j in sorted(glob.glob(os.path.join(root, "**", "*_ipsae.json"), recursive=True)):
        cif = Path(j.replace("_ipsae.json", ".cif"))
        if not cif.exists():
            continue
        design = re.sub(r"(_seed\d+)?_sample_\d+$", "", cif.stem)
        pairs = pair_minima(j, cif)
        if pairs is None:
            void.setdefault(design, "ipsae produced no output")
            continue
        vals, info = aggregate(pairs, target, binder, struct=cif, pick=pick.get(design))
        if vals is None:
            void.setdefault(design, info)
            continue
        pick.setdefault(design, info)
        for k, v in vals.items():
            by[design][k].append(v)

    print(f"{'design':26s} " + " ".join(f"{r:>12s}" for r in RULES) + "   pick   n")
    rows = {}
    for d, m in sorted(by.items()):
        rows[d] = {k: max(v) for k, v in m.items()}
        print(f"{d:26s} " + " ".join(f"{rows[d][r]:12.4f}" for r in RULES)
              + f"   {':'.join(pick[d])}   {len(m['max_pair'])}")
    for d, why in sorted(void.items()):
        print(f"🔴 VOID {d}: {why}")
    print(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
