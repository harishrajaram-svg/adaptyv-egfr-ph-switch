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
    print("design_inventory.py --selftest PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs/mosaic-p2")
    ap.add_argument("--tsv", action="store_true", help="print every row, not just the summary")
    ap.add_argument("--selftest", action="store_true")
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
    print(f"lengths present                 {s['lengths']}")


if __name__ == "__main__":
    sys.exit(main())
