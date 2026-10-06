#!/usr/bin/env python3
"""Run the repacked dddG_elec leg over every fixture that already has a rigid answer.

The order is deliberate: the Schroter series FIRST, because it is the only set where the
repacked leg can overturn a conclusion already in the challenge file (s6b's 0 of 5). The anchor
sets test P2 and P3 and cannot overturn anything -- they can only resize s9's effect.
"""
import argparse
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dddg_elec_repacked as dd

STRUCT = os.path.join(HERE, "structures")
SETS = [
    ("schroter", "H,L", "A", [
        ("adalimumab_WT", os.path.join(STRUCT, "3WD5_adalimumab_WT.pdb")),
        ("PSV1", os.path.join(STRUCT, "3WD5_PSV1.pdb")),
        ("PSV2", os.path.join(STRUCT, "3WD5_PSV2.pdb")),
        ("PSV3", os.path.join(STRUCT, "3WD5_PSV3.pdb")),
        ("scramble_neg", os.path.join(STRUCT, "3WD5_scramble.pdb")),
    ]),
    ("anchor_K166", "Z", "A,B,C",
     [(f"K166#{i}", os.path.join(HERE, f"pose_anchor_K166_rebuilt_{i}.pdb")) for i in range(1, 6)]),
    ("anchor_K166Q", "Z", "A,B,C",
     [(f"K166Q#{i}", os.path.join(HERE, f"pose_anchor_K166Q_{i}.pdb")) for i in range(1, 6)]),
    ("far_control", "Z", "A,B,C",
     [(f"far#{i}", os.path.join(HERE, f"pose_far_control_{i}.pdb")) for i in range(1, 6)]),
]
MEASURED = {"adalimumab_WT": "9x", "PSV1": "231x", "PSV2": "785x", "PSV3": "505x",
            "scramble_neg": "-"}
LEGS = ("rigid", "repacked_local")   # overridable with --legs
# repacked_all is dropped from the n=5 trial sweep on cost, not on principle: it repacks ~230-270
# residues per trajectory against repacked_local's 9-13, and its n=1 numbers are already on
# record in the previous run. The decisive leg is repacked_local, which is the one P2 is about.


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", help="comma-separated subset of set names")
    ap.add_argument("--legs", help="comma-separated subset of legs")
    ap.add_argument("--out", default="repacked_sweep.json")
    a = ap.parse_args()
    sets = SETS if not a.sets else [x for x in SETS if x[0] in a.sets.split(",")]
    legs = LEGS if not a.legs else tuple(a.legs.split(","))
    globals()["LEGS"] = legs
    results = {}
    for setname, binder, target, items in sets:
        print(f"\n### {setname}   binder {binder} / target {target}")
        print(f"{'case':<18}{'measured':>10}"
              + "".join(f"{l + ' med [min..max]':>32}" for l in LEGS)
              + f"{'n_his':>7}{'repack n/rmsd':>18}")
        print("-" * 128)
        for label, path in items:
            if not os.path.exists(path):
                print(f"{label:<18}{'':>10}  MISSING {path}")
                continue
            try:
                r = dd.score(path, binder.replace(",", ""), target.replace(",", ""),
                             legs_wanted=LEGS)
            except SystemExit as e:
                print(f"{label:<18}{'':>10}  REFUSED: {e}")
                results.setdefault(setname, {})[label] = {"refused": str(e)}
                continue
            results.setdefault(setname, {})[label] = r
            cells = "".join(
                f"{r[l]['dddG']:>+9.3f} [{r[l]['dddG_min']:+.2f}..{r[l]['dddG_max']:+.2f}]"
                + ("!" if r[l]["sign_unstable"] else " ") for l in LEGS)
            loc = next((r[l] for l in legs if l != "rigid" and l in r), None) or {}
            print(f"{label:<18}{MEASURED.get(label, ''):>10}{cells}{r['n_his']:>7}"
                  f"{str(loc.get('repacked_residues')) + '/' + str(loc.get('sidechain_rmsd')):>18}")

    print("\n\n### MEDIANS and the `>= 0` pass rate per leg")
    print(f"{'set':<18}" + "".join(f"{l + ' med':>20}" for l in LEGS)
          + "".join(f"{l + ' >=0':>20}" for l in LEGS))
    print("-" * 100)
    summary = {}
    for setname, cases in results.items():
        meds, rates = [], []
        for leg in LEGS:
            v = [c[leg]["dddG"] for c in cases.values() if leg in c]
            meds.append(float(np.median(v)) if v else float("nan"))
            rates.append(f"{sum(1 for x in v if x >= 0)}/{len(v)}" if v else "-")
        summary[setname] = {"medians": dict(zip(LEGS, meds)),
                            "pass_ge_zero": dict(zip(LEGS, rates))}
        print(f"{setname:<18}" + "".join(f"{m:>+20.4f}" for m in meds)
              + "".join(f"{r:>20}" for r in rates))

    out = os.path.join(HERE, a.out)
    json.dump({"legs": list(legs), "shell_A": dd.SHELL, "summary": summary,
               "results": results}, open(out, "w"), indent=1, default=str)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
