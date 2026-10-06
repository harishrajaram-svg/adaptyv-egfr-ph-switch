#!/usr/bin/env python3
"""Score the problem-2 validation gate. The frozen scorer could not, and that is a finding.

\U0001F534 bin/ipsae_min.py REFUSES this target, correctly. It asserts that exactly ONE
inter-chain pair is present, because ipSAE_min is a min over the two ALIGNMENT DIRECTIONS of one
interface, and a flat min() across chain PAIRS would let an intra-binder interface (the VH:VL of
a badly packed Fv) be reported as the binding score. That guard earned its keep on problem 1's
g532 ladder.

Problem 2's target is a TRIMER. Binder D against protomers A, B, C gives six inter-chain pairs,
so the frozen scorer stops. THE INSTRUMENT WAS BUILT FOR TWO-CHAIN COMPLEXES AND THE PROBLEM-2
TARGET IS NOT ONE. Nothing in frozen-decisions.md anticipated that, and it would have surfaced at
ranking time on a deadline rather than on day 2.

This does NOT change the metric. ipSAE_min per interface is untouched: min over the two
directions of one binder:protomer pair, from the same vendored ipsae.py. What this file adds is
the one thing the trimer forces -- a stated rule for combining three binder:protomer interfaces
into one number per design -- and it reports every candidate rule side by side rather than
picking one quietly. frozen-decisions.md freezes the instrument; it leaves "target trimming"
and per-week application open, and this is that.

THE CANDIDATE RULES, and why the choice is not obvious:
  max        the best-engaged protomer. s1 says the epitope STRADDLES two protomers, so a
             correct binder should light up two interfaces; max reports the stronger.
  second     the SECOND-best protomer. A binder that straddles must have two real interfaces,
             so this is the one that distinguishes straddling from a single-protomer graze --
             which is exactly the failure mode s1's geometry predicts.
  sum        total engagement across the trimer.
  min_engaged  min over protomers scoring above zero; conservative, and undefined when only
             one engages.

The PASS RULE was pre-registered in s7 before any score existed: the positive must score clearly
above all four composition-matched negatives. It is evaluated below under EVERY rule, because a
pass that depends on the aggregation choice is not a pass.
"""
import glob, json, math, os, re, subprocess, sys, collections
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(ROOT / "bin"))
import ipsae_min as im

BINDER_DEFAULT = "D"


def pairs_for(pae_file: Path, struct: Path, pae_cut=10, dist_cut=10):
    """{frozenset(chain pair): [both directional ipSAE values]} from the vendored ipsae.py."""
    pae_file, struct = pae_file.resolve(), struct.resolve()
    out = struct.with_name(f"{struct.stem}_{pae_cut}_{dist_cut}.txt")
    if out.exists():
        out.unlink()
    r = subprocess.run([str(im.PY), str(im.IPSAE), str(pae_file), str(struct),
                        str(pae_cut), str(dist_cut)],
                       capture_output=True, text=True, cwd=struct.parent)
    if r.returncode != 0 or not out.exists():
        return None
    asym = {}
    for row in (l.split() for l in out.read_text().splitlines()
                if l.strip() and not l.startswith("Chn1")):
        if len(row) > 5 and row[4] == "asym":
            try:
                v = float(row[5])
            except ValueError:
                continue
            if math.isfinite(v):
                asym[(row[0], row[1])] = v
    pairs = {}
    for (c1, c2), v in asym.items():
        pairs.setdefault(frozenset((c1, c2)), []).append(v)
    return pairs


def per_protomer(pairs, binder_chains):
    """ipSAE_min for each binder:TARGET interface. Min over the two directions, unchanged.

    \U0001F534 `binder_chains` is a SET, and INTRA-BINDER pairs are excluded. The first version
    took a single chain and the Fab leg immediately reproduced the exact failure that
    bin/ipsae_min.py's guard exists to prevent: with binder="D" the D:E pair -- the Fab's own
    VH:VL interface -- was reported as a binder:target interface at 0.8334, HIGHER than any real
    contact with TNF (0.52-0.55). A badly packed Fv would have been scored as excellent binding.
    The guard's docstring names this case on problem 1's g532 ladder; it took eleven minutes to
    walk into it on problem 2."""
    out = {}
    for p, vals in pairs.items():
        if len(vals) != 2:
            continue
        inside = p & binder_chains
        if len(inside) != 1:
            continue                    # 0 = target:target, 2 = INTRA-BINDER. Neither is ours.
        other = next(c for c in p if c not in binder_chains)
        out.setdefault(other, []).append(min(vals))
    # a two-chain binder touches each protomer twice (VH and VL); take the stronger
    return {c: max(v) for c, v in out.items()}


def aggregates(d):
    v = sorted(d.values(), reverse=True)
    eng = [x for x in v if x > 0]
    return {
        "max": v[0] if v else 0.0,
        "second": v[1] if len(v) > 1 else 0.0,
        "sum": float(sum(v)),
        "min_engaged": min(eng) if eng else 0.0,
        "n_engaged": len(eng),
    }


def main():
    binder = sys.argv[1] if len(sys.argv) > 1 else BINDER_DEFAULT
    base = ROOT / "runs" / "gate-tnf" / "tnf_gate_full"
    per_case = collections.defaultdict(list)
    for j in sorted(glob.glob(str(base / "**" / "*_ipsae.json"), recursive=True)):
        cif = Path(j.replace("_ipsae.json", ".cif"))
        if not cif.exists():
            continue
        case = re.sub(r"(_seed\d+)?_sample_\d+$", "", cif.stem)
        b = {"D", "E"} if case.startswith("pos2") else {"D"}
        pp = pairs_for(Path(j), cif)
        if pp is None:
            print(f"FAILED {cif.name}")
            continue
        d = per_protomer(pp, b)
        if not d:
            print(f"  {cif.name}: no interface involving {sorted(b)}; pairs "
                  f"{[':'.join(sorted(p)) for p in pp]}")
            continue
        per_case[case].append(d)

    print(f"\n{'case':<36}{'n':>3}  per-protomer ipSAE_min, median over seeds")
    print("-" * 92)
    summary = {}
    for case, lst in sorted(per_case.items()):
        chains = sorted({c for d in lst for c in d})
        med = {c: float(np.median([d.get(c, 0.0) for d in lst])) for c in chains}
        agg = aggregates(med)
        summary[case] = agg
        cells = "  ".join(f"{c}={med[c]:.4f}" for c in chains)
        print(f"{case:<36}{len(lst):>3}  {cells}")
    print("-" * 92)
    print(f"\n{'case':<36}" + "".join(f"{k:>14}" for k in ("max", "second", "sum", "min_engaged"))
          + "   engaged")
    print("-" * 92)
    for case, a in sorted(summary.items(), key=lambda kv: -kv[1]["max"]):
        print(f"{case:<36}" + "".join(f"{a[k]:>14.4f}" for k in
                                      ("max", "second", "sum", "min_engaged"))
              + f"{a['n_engaged']:>10}")

    pos = summary.get("pos_tnf_tnfr2")
    negs = {k: v for k, v in summary.items() if k.startswith("neg_")}
    print("\n" + "=" * 92)
    print("PRE-REGISTERED PASS RULE (s7, written before any score existed):")
    print("  PASS = pos_tnf_tnfr2 scores CLEARLY ABOVE all four neg_shuffled_tnfr2_*.")
    print("=" * 92)
    if not pos or not negs:
        print("  cannot evaluate: positive or negatives missing")
        return
    for rule in ("max", "second", "sum", "min_engaged"):
        worst = max(n[rule] for n in negs.values())
        ok = pos[rule] > worst
        margin = pos[rule] - worst
        print(f"  {rule:<13} pos {pos[rule]:.4f}  vs  best negative {worst:.4f}   "
              f"margin {margin:+.4f}   {'PASS' if ok else 'FAIL'}")
    print("\nA pass that holds under every rule is a pass. One that depends on the aggregation")
    print("choice is an artifact of that choice and must be reported as such.")


if __name__ == "__main__":
    main()
