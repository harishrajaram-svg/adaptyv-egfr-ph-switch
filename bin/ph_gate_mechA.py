#!/usr/bin/env python3
"""Mechanism A pH gate: a histidine on the BINDER, paired with a carboxylate on the TARGET.

Mechanism B (bin/ph_gate.py) exploits the target's own H433 and is the only mechanism
tested so far -- one histidine, one location. Measured 2026-10-02: burying H433 drops its
pKa (binder contacts vs bound pKa r = -0.519), so B fights desolvation by construction.

Mechanism A inverts the geometry. The titratable group lives on the binder, where we
control its environment, and the stabilising negative charge is a target Asp/Glu. Same
thermodynamic linkage, different sign of control:

    pKa_bound > pKa_free  ->  protonation favours binding  ->  acid STRENGTHENS  (wanted)

The difference that matters for scoring: pKa_free is NOT a constant here. It is measured
on the unbound binder, per design, per histidine -- so this needs PROPKA twice.

No GPU. Structures already exist; PROPKA is CPU.

Usage
-----
    python3 bin/ph_gate_mechA.py --complex-dir <dir> --binder-dir <dir> [--tsv out.tsv]
"""
import argparse, csv, os, re, shutil, subprocess, sys, tempfile

PH_LO, PH_HI = 6.5, 7.4
RATIO_BAR = 1.20     # below this the "switch" is inside PROPKA's own noise (matches ph_gate.py)
BINDER_CHAIN = "A"   # default: BoltzGen writes the design first
TARGET_CHAIN = "B"


def ratio(pka_free, pka_bound, ph_lo=PH_LO, ph_hi=PH_HI):
    """K(ph_lo)/K(ph_hi). >1 = stronger in acid."""
    K = lambda ph: (1 + 10 ** (pka_bound - ph)) / (1 + 10 ** (pka_free - ph))
    return K(ph_lo) / K(ph_hi)


def propka(pdb_path):
    """Run PROPKA3, return {(resname, resnum, chain): pKa}."""
    with tempfile.TemporaryDirectory() as wd:
        local = os.path.join(wd, os.path.basename(pdb_path))
        shutil.copy(pdb_path, local)
        subprocess.run([sys.executable, "-m", "propka", local],
                       capture_output=True, text=True, cwd=wd)
        pka_file = os.path.splitext(local)[0] + ".pka"
        if not os.path.exists(pka_file):
            return None
        out, in_summary = {}, False
        for line in open(pka_file):
            if line.startswith("SUMMARY"):
                in_summary = True
                continue
            if in_summary:
                if line.startswith("-----"):
                    break
                f = line.split()
                if len(f) >= 4:
                    try:
                        out[(f[0], int(f[1]), f[2])] = float(f[3])
                    except ValueError:
                        pass
        return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--complex-dir", required=True)
    ap.add_argument("--binder-dir", required=True)
    ap.add_argument("--tsv")
    # BoltzGen writes the design as chain A; BindCraft writes the target as A and the
    # binder as B. Scoring the wrong chain silently returns "no histidines found".
    ap.add_argument("--binder-chain", default=BINDER_CHAIN)
    a = ap.parse_args()
    binder_chain = a.binder_chain

    rows = []
    names = sorted(x[:-4] for x in os.listdir(a.complex_dir) if x.endswith(".pdb"))
    for n in names:
        cpx = propka(os.path.join(a.complex_dir, n + ".pdb"))
        bnd = propka(os.path.join(a.binder_dir, n + ".pdb"))
        if not cpx or not bnd:
            print(f"{n}: PROPKA failed", file=sys.stderr)
            continue
        # every histidine on the binder chain, scored free vs bound
        for (rn, num, ch), free in sorted(bnd.items()):
            if rn != "HIS" or ch != binder_chain:
                continue
            bound = cpx.get((rn, num, binder_chain))
            if bound is None:
                continue
            r = ratio(free, bound)
            rows.append(dict(design=n, his=num, pka_free=round(free, 2),
                             pka_bound=round(bound, 2), delta=round(bound - free, 2),
                             ratio=round(r, 3),
                             passes=r >= RATIO_BAR))
        print(f"  scored {n}", file=sys.stderr)

    rows.sort(key=lambda r: -r["ratio"])
    print(f"{'design':<34} {'His':>5} {'pKa_free':>9} {'pKa_bound':>10} {'delta':>7} {'ratio':>7}  verdict")
    print("-" * 92)
    for r in rows:
        if r["ratio"] >= RATIO_BAR:
            v = "PASS - real switch"
        elif r["ratio"] > 1.0:
            v = f"fail - rises but under {RATIO_BAR:.2f}x (untouched or marginal)"
        else:
            v = "fail - acid weakens"
        print(f'{r["design"]:<34} {r["his"]:>5} {r["pka_free"]:>9.2f} {r["pka_bound"]:>10.2f} '
              f'{r["delta"]:>+7.2f} {r["ratio"]:>6.2f}x  {v}')
    if rows:
        n = sum(1 for r in rows if r["passes"])
        print("-" * 92)
        weak = sum(1 for r in rows if 1.0 < r["ratio"] < RATIO_BAR)
        print(f"{n}/{len(rows)} binder histidines are real switches (ratio >= {RATIO_BAR:.2f}x) "
              f"across {len({r['design'] for r in rows})} designs; "
              f"{weak} rise but stay under the bar.")
    if a.tsv and rows:
        with open(a.tsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"wrote {a.tsv}")


if __name__ == "__main__":
    main()
