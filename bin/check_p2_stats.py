#!/usr/bin/env python3
"""Gate every condition statistic in the problem-2 methods document against the run tables.

Exists because of a real error: two cells in the condition comparison reported a ROW
VALUE where a median belonged (0.153 and 0.183 instead of 0.136 and 0.179), and the
derived claim moved with them, from +31% to a stated +20%. Both wrong figures were real
numbers from the same tables, so every spot-check of a cell against the data succeeded.
Only recomputation catches that class, and check_claims.py reads problem-1 docs only.

Usage:  check_p2_stats.py            # exit 1 on any disagreement
        check_p2_stats.py --selftest
"""
import csv
import glob
import os
import re
import statistics
import sys

DOC = "submissions/02-tnf-METHODS.md"
RUNS = "runs/mosaic-p2"

# condition number in the s6.3 table -> the run directories that produced it
CONDITIONS = {
    1: ["p2dry04"],
    2: ["p2probe-a", "p2probe-b", "p2probe-c", "p2probe-d"],
    3: ["p2probe2-e", "p2probe2-f", "p2probe2-g", "p2probe2-h"],
    4: ["p2deep-p", "p2deep-q", "p2deep-r", "p2deep-s"],
    5: ["p2free-p", "p2free-q", "p2free-r", "p2free-s"],
}


def rows_for(dirs):
    out = []
    for d in dirs:
        f = os.path.join(RUNS, d, "designs.tsv")
        if not os.path.isfile(f):
            continue
        with open(f) as fh:
            out += list(csv.DictReader(fh, delimiter="\t"))
    return out


def stats(rows):
    v = [float(r["iptm_repred"]) for r in rows if r.get("iptm_repred")]
    h = [float(r["his_N_to_cation_N"]) for r in rows if r.get("his_N_to_cation_N")]
    if not v:
        return None
    return {"median": statistics.median(v), "max": max(v),
            "closest_his": min(h) if h else None, "n": len(v)}


def parse_table(text):
    """Pull {condition number: {median, max, closest_his}} from the s6.3 table."""
    got = {}
    for line in text.splitlines():
        m = re.match(r"\|\s*(\d) — ", line)
        if not m:
            continue
        cells = [c.strip().replace("**", "") for c in line.strip("|").split("|")]
        if len(cells) < 7:
            continue
        try:
            # keep the RAW text too: the number of decimals shown sets the tolerance,
            # so a doc writing 18.0 is held to +/-0.05 and one writing 18.04 to +/-0.005
            got[int(m.group(1))] = {
                "median": cells[4],
                "max": cells[5],
                "closest_his": cells[6].replace("Å", "").strip(),
            }
            {k: float(v) for k, v in got[int(m.group(1))].items()}
        except ValueError:
            got.pop(int(m.group(1)), None)
            continue
    return got


def decimals(text):
    """How many decimal places the document actually displays."""
    return len(text.split(".")[1]) if "." in text else 0


def close(stated_text, actual):
    """Is the displayed figure a correct rounding of the actual one?

    Tolerance is half a unit of the LAST DISPLAYED PLACE. Not round(actual, n):
    0.1555 displays correctly as 0.156 by half-up rounding, while Python's round
    returns 0.155, and holding the document to Python's tie-breaking would fail a
    cell that is right.
    """
    d = decimals(stated_text)
    return abs(float(stated_text) - actual) <= 0.5 * 10 ** -d + 1e-12


INVENTORY_LABELS = {
    "sequences, all distinct": "distinct_sequences",
    "clearing the 0.45 interface gate": "clear_iptm_gate",
    "histidine fraction above the 8% cap": "over_his_cap",
    "nearest histidine within the 4.0 \u00c5 criterion": "within_placement_bar",
}


def parse_inventory(text):
    """Pull {label: stated count} from the s10 inventory table."""
    got = {}
    for line in text.splitlines():
        cells = [c.strip().replace("**", "") for c in line.strip().strip("|").split("|")]
        if len(cells) != 2:
            continue
        label = cells[0]
        if label in INVENTORY_LABELS:
            m = re.match(r"^(\d+)", cells[1])
            if m:
                got[label] = int(m.group(1))
    return got


def check_inventory(doc_text, all_rows):
    """The s10 counts must equal design_inventory's, which reads the same tables."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "di", os.path.join("analysis", "02-tnf", "design_inventory.py"))
    di = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(di)
    actual = di.summarise(all_rows)
    fails = []
    stated = parse_inventory(doc_text)
    if not stated:
        return ["could not parse the s10 inventory table at all"], {}
    for label, key in INVENTORY_LABELS.items():
        if label not in stated:
            continue
        if stated[label] != actual[key]:
            fails.append(f"s10 inventory '{label}': doc says {stated[label]}, "
                         f"recomputed {actual[key]}")
    return fails, stated


def check(doc_text, condition_rows):
    fails = []
    stated = parse_table(doc_text)
    if not stated:
        fails.append("could not parse the condition table at all")
        return fails
    for cond, want in sorted(stated.items()):
        rows = condition_rows.get(cond, [])
        got = stats(rows)
        if got is None:
            fails.append(f"condition {cond}: stated in the doc, no run data on disk")
            continue
        for key in ("median", "max", "closest_his"):
            a = got[key]
            if a is None:
                continue
            if not close(want[key], a):
                fails.append(f"condition {cond} {key}: doc says {want[key]}, "
                             f"recomputed {a:.4f} from n={got['n']}")
    return fails


def selftest():
    doc = (
        "| condition | pH weight | binding weights | steps | median | max | closest His |\n"
        "|---|---|---|---|---|---|---|\n"
        "| 1 — baseline | 2.0 | 1.0 | 50 | 0.136 | 0.156 | 18.0 Å |\n"
    )
    rows = {1: [
        {"iptm_repred": "0.1527", "his_N_to_cation_N": "18.04"},
        {"iptm_repred": "0.1198", "his_N_to_cation_N": "24.82"},
        {"iptm_repred": "0.1107", "his_N_to_cation_N": "47.93"},
        {"iptm_repred": "0.1555", "his_N_to_cation_N": "39.20"},
    ]}
    assert check(doc, rows) == [], check(doc, rows)
    # MUTATION TEST: the exact historical error must be caught. 0.1527 is a REAL row
    # value from this very table, which is why reading the prose never found it.
    bad = doc.replace("| 50 | 0.136 |", "| 50 | 0.153 |")
    f = check(bad, rows)
    assert any("median" in x for x in f), f
    # a correct half-up rounding must PASS: max is 0.1555, displayed as 0.156,
    # where Python's round() would say 0.155 and fail a cell that is right
    assert close("0.156", 0.1555) and close("0.155", 0.1555)
    assert not close("0.153", 0.13625)
    assert close("18.0", 18.04) and not close("18.04", 18.1)
    # and a wrong max, and a wrong distance
    assert any("max" in x for x in check(doc.replace("0.156 |", "0.200 |"), rows))
    assert any("closest_his" in x for x in check(doc.replace("18.0 Å", "2.62 Å"), rows))
    # a condition with no data on disk must fail loudly, not pass silently
    assert any("no run data" in x for x in check(doc, {}))
    # the s10 inventory parser: a stated count must be read, and a wrong one rejected
    inv = ("| | count of 3 |\n|---|---|\n| sequences, all distinct | 3 |\n"
           "| clearing the 0.45 interface gate | 0 |\n")
    assert parse_inventory(inv) == {"sequences, all distinct": 3,
                                    "clearing the 0.45 interface gate": 0}, parse_inventory(inv)
    assert parse_inventory("| nothing | here |") == {}
    print("check_p2_stats.py --selftest PASS")


def main():
    if "--selftest" in sys.argv:
        return selftest()
    if not os.path.isfile(DOC):
        sys.exit(f"run me from the repo root: {DOC} not found")
    rows = {c: rows_for(d) for c, d in CONDITIONS.items()}
    doc = open(DOC).read()
    fails = check(doc, rows)
    # the s10 inventory counts EVERY run on disk, not only the five conditions
    all_rows = []
    for f in sorted(glob.glob(os.path.join(RUNS, "*", "designs.tsv"))):
        with open(f) as fh:
            all_rows += list(csv.DictReader(fh, delimiter="\t"))
    inv_fails, inv_stated = check_inventory(doc, all_rows)
    fails += inv_fails
    stated = parse_table(doc)
    for c in sorted(stated):
        s = stats(rows.get(c, []))
        n = s["n"] if s else 0
        print(f"  condition {c}: {n} run(s) on disk, "
              f"median {s['median']:.4f} max {s['max']:.4f} closest {s['closest_his']:.2f} Å"
              if s else f"  condition {c}: NO DATA")
    if fails:
        print("\nFAIL:")
        for f in fails:
            print(f"  {f}")
        sys.exit(1)
    print(f"  s10 inventory: {len(inv_stated)} stated count(s), all recomputed from "
          f"{len(all_rows)} rows on disk")
    print(f"\nPASS: {len(stated)} condition(s) and the s10 inventory, every median, max, "
          f"distance and count recomputed from designs.tsv.")


if __name__ == "__main__":
    sys.exit(main())
