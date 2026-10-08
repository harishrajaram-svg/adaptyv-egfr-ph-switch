#!/usr/bin/env python3
"""Count what candidate material actually exists, across every run on disk.

The methods document states this count, so it has to be regenerable rather than
remembered. Tolerates columns that did not exist in older runs -- frac_H was added
partway through, so a missing value reads n/a and is never silently treated as 0.

Usage:  design_inventory.py [--runs runs/mosaic-p2] [--tsv]
        design_inventory.py --selftest
"""
import argparse
import csv
import glob
import os
import sys

IPTM_GATE = 0.45     # D-P2-1: below this, a geometry verdict is withdrawn to n/a
PLACE_BAR = 4.0      # the placement criterion, in angstroms
HIS_CAP = 0.08       # the composition cap the loss enforces as a soft hinge

# Monomer-foldability floor, adopted 2026-10-07. NOT ours: it is the default in
# reference/anthropic-binder-design-protocol.md:100, where it is specified as a
# PRE-SCORING filter to run before any co-folding spend. This project never applied
# it. analysis/02-tnf/calibrate_external.py measures why it matters: across 150
# TNF-alpha designs that were actually assayed, 146/150 clear it, and the lowest
# value among them is 68.9 on Boltz-2. Our own best is 68.8. Reported on the 0-1
# scale because that is what our designs.tsv carries.
MONOMER_FLOOR = 0.70


def num(row, key):
    """Float, or None when the column is absent or unparseable. Never 0.0 by default."""
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return None


def load(runs_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(runs_dir, "*", "designs.tsv"))):
        cond = os.path.basename(os.path.dirname(f))
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r["_condition"] = cond
                rows.append(r)
    return rows


def summarise(rows):
    seqs = {r.get("sequence") for r in rows if r.get("sequence")}
    def count(pred):
        return sum(1 for r in rows if pred(r))
    return {
        "rows": len(rows),
        "distinct_sequences": len(seqs),
        "clear_iptm_gate": count(lambda r: (num(r, "iptm_repred") or 0) >= IPTM_GATE),
        "within_placement_bar": count(
            lambda r: num(r, "his_N_to_cation_N") is not None
            and num(r, "his_N_to_cation_N") <= PLACE_BAR),
        "over_his_cap": count(lambda r: (num(r, "frac_H") or 0) > HIS_CAP),
        "frac_H_unrecorded": count(lambda r: num(r, "frac_H") is None),
        "clear_monomer_floor": count(
            lambda r: num(r, "plddt_binder_repred") is not None
            and num(r, "plddt_binder_repred") >= MONOMER_FLOOR),
        "monomer_unrecorded": count(lambda r: num(r, "plddt_binder_repred") is None),
        "lengths": sorted({int(num(r, "length")) for r in rows if num(r, "length")}),
    }


def selftest():
    """A missing frac_H must not be counted as compliant, nor as a breach."""
    rows = [
        {"sequence": "AAA", "iptm_repred": "0.50", "his_N_to_cation_N": "3.0",
         "frac_H": "0.02", "length": "76"},                      # clears everything
        {"sequence": "BBB", "iptm_repred": "0.20", "his_N_to_cation_N": "2.6",
         "frac_H": "0.09", "length": "84"},                      # under bar, over cap, no gate
        {"sequence": "CCC", "iptm_repred": "0.10", "length": "76"},  # frac_H and distance absent
    ]
    s = summarise(rows)
    assert s["rows"] == 3 and s["distinct_sequences"] == 3, s
    assert s["clear_iptm_gate"] == 1, s
    assert s["within_placement_bar"] == 2, s          # the absent distance must NOT count
    assert s["over_his_cap"] == 1, s                  # the absent frac_H must NOT count
    assert s["frac_H_unrecorded"] == 1, s             # but it must be reported as absent
    assert s["lengths"] == [76, 84], s
    # duplicate sequences must collapse
    assert summarise(rows + [dict(rows[0])])["distinct_sequences"] == 3
    # the monomer floor: absent plddt counts as neither pass nor fail, as with frac_H
    f = summarise([
        {"sequence": "A", "plddt_binder_repred": "0.81", "length": "76"},   # clears
        {"sequence": "B", "plddt_binder_repred": "0.47", "length": "76"},   # fails
        {"sequence": "C", "length": "76"},                                  # absent
    ])
    assert f["clear_monomer_floor"] == 1, f
    assert f["monomer_unrecorded"] == 1, f
    # MUTATION: a floor read on the 0-100 scale would pass everything. Our column is
    # 0-1, so a 0-100 value must NOT be silently accepted as clearing a 0.70 floor.
    g = summarise([{"sequence": "A", "plddt_binder_repred": "47.0", "length": "76"}])
    assert g["clear_monomer_floor"] == 1, "a 0-100 value is numerically above 0.70"
    assert MONOMER_FLOOR == 0.70, "the floor is declared on the 0-1 scale"
    # so the guard is the EXPORT refusing mixed scales, tested below
    assert _scale_ok([0.32, 0.47, 0.69]) and not _scale_ok([0.47, 85.3]), "scale guard"
    print("design_inventory.py --selftest PASS")


def _scale_ok(vals):
    """All plddt values must sit on ONE scale. Mixing 0-1 and 0-100 silently inverts
    every floor comparison, so the exporter refuses rather than guesses."""
    vals = [v for v in vals if v is not None]
    if not vals:
        return True
    return max(vals) <= 1.0 or min(vals) > 1.0


def export(rows, faa, tsv):
    """Regenerate outbox/02-tnf-candidates.{faa,tsv} from what is on disk.

    The committed export went stale by 10 designs on 2026-10-07 -- two arms landed
    after it was written and it was missing the project's own best row. Nothing
    regenerated it because nothing could; this is that missing step.
    """
    if not _scale_ok([num(r, "plddt_binder_repred") for r in rows]):
        sys.exit("REFUSE: plddt_binder_repred mixes 0-1 and 0-100 scales across runs")
    cols = ["design", "condition", "length", "seed", "iptm_repred", "iptm_design",
            "plddt_binder_repred", "his_N_to_cation_N", "closest_his", "frac_H",
            "frac_V", "frac_G", "n_C", "sequence"]
    ordered = sorted(rows, key=lambda r: -(num(r, "iptm_repred") or -1))
    with open(tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in ordered:
            w.writerow([r.get("condition") or r["_condition"] if c == "condition"
                        else r.get(c, "n/a") for c in cols])
    with open(faa, "w") as fh:
        for r in ordered:
            def f(k, fmt="{:.4f}"):
                v = num(r, k)
                return "n/a" if v is None else fmt.format(v)
            fh.write(f">{r['design']} condition={r['_condition']} "
                     f"length={r.get('length','?')} seed={r.get('seed','?')} "
                     f"iptm_repred={f('iptm_repred')} "
                     f"plddt_binder_repred={f('plddt_binder_repred','{:.3f}')} "
                     f"hisN={f('his_N_to_cation_N','{:.2f}')} "
                     f"frac_H={f('frac_H','{:.3f}')} "
                     f"monomer_floor={'PASS' if (num(r,'plddt_binder_repred') or 0) >= MONOMER_FLOOR else 'FAIL'} "
                     f"geometry_verdict={r.get('geometry_read') or 'n/a'}\n")
            fh.write(f"{r['sequence']}\n")
    print(f"exported {len(ordered)} designs -> {faa} + {tsv}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs/mosaic-p2")
    ap.add_argument("--tsv", action="store_true", help="print every row, not just the summary")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--export", action="store_true",
                    help="regenerate outbox/02-tnf-candidates.{faa,tsv} from disk")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    rows = load(a.runs)
    if not rows:
        sys.exit(f"no designs.tsv under {a.runs}")
    if a.tsv:
        print("condition\tdesign\tiptm_repred\thisN\tfrac_H\tplddt_binder")
        for r in sorted(rows, key=lambda r: num(r, "his_N_to_cation_N") or 1e9):
            def f(k, fmt="{:.3f}"):
                v = num(r, k)
                return "n/a" if v is None else fmt.format(v)
            print(f"{r['_condition']}\t{r['design']}\t{f('iptm_repred','{:.4f}')}\t"
                  f"{f('his_N_to_cation_N','{:.2f}')}\t{f('frac_H')}\t"
                  f"{f('plddt_binder_repred','{:.2f}')}")
        print()
    s = summarise(rows)
    print(f"design rows                     {s['rows']}")
    print(f"distinct sequences              {s['distinct_sequences']}")
    print(f"clearing the {IPTM_GATE} iptm gate     {s['clear_iptm_gate']}")
    print(f"within the {PLACE_BAR} A placement bar  {s['within_placement_bar']}")
    print(f"above the {HIS_CAP:.0%} histidine cap     {s['over_his_cap']}"
          f"   ({s['frac_H_unrecorded']} rows predate the frac_H column)")
    print(f"clearing the {MONOMER_FLOOR} monomer floor  {s['clear_monomer_floor']}"
          f"   ({s['monomer_unrecorded']} rows have no plddt_binder_repred)")
    print(f"lengths present                 {s['lengths']}")
    if a.export:
        export(rows, "outbox/02-tnf-candidates.faa", "outbox/02-tnf-candidates.tsv")


if __name__ == "__main__":
    sys.exit(main())
