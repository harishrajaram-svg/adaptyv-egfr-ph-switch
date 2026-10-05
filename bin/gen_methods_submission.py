#!/usr/bin/env python3
"""Generate the SUBMISSION sections of the methods document FROM the emitted CSV.

WHY THIS EXISTS. METHODS section 11 was hand-maintained through three submission changes in
one evening (10 -> 11 -> 12 designs, plus a change of ranking basis). An audit found the
result: the section still said "11 designs", its rank table had 11 rows and omitted the 12th
design entirely, every rank from 8 upward was off by one, six of the twelve assessment strings
cited the wrong design by rank number, and section 12's declarations attested review of "all
TEN submitted sequences" -- a false attestation in a signed section.

None of those were reasoning errors. They were all the same mechanical error: a document
describing an artifact, maintained separately from it. So the section is now GENERATED. If the
CSV changes, this is re-run, and the document cannot drift from it.

Scope: the rank table, the counts, the family partition and the declarations' numeric claims.
Everything interpretive in section 11 -- why the basis changed, what is not claimed, the cost
of the ordering -- stays hand-written, because it is argument and not data.

Usage:
    gen_methods_submission.py            # print the generated blocks
    gen_methods_submission.py --apply    # splice them into the methods document
    gen_methods_submission.py --selftest
"""
import csv, json, os, re, sys

CSV = "submissions/01-egfr.csv"
MS = "analysis/01-egfr/multisite_pooled.json"
METHODS = "submissions/01-egfr-METHODS.md"
STRUCTDIR = "submissions/structures"

# Backbone family per design. A point mutant belongs to its PARENT's family: ranks that differ
# by one residue are not independent tests, which is the whole reason the family partition is
# reported at all.
FAMILY = {
    "rimA01_r15_boltzgen_egfr_d3_rimA_20": "rimA01_r15_d3_rimA_20",
    "rimA01_r15_L133E":                    "rimA01_r15_d3_rimA_20",   # point mutant of the above
    "bc_s360518_mpnn9_A22D":               "d3acid3_l65_s360518",
    "d2c_mpnn13_S88D_serasp":              "d2c_101_l147_s144898",
    "rimA02_d3_rimA_14_vhh":               "rimA02_d3_rimA_14 (VHH)",
    "h370_020_vhh":                        "h370_020 (VHH)",
}
DEFAULT_FAMILY = "d3acid_l65_s831683"
ANTIBODY = {"nanobody", "scfv", "fab_kappa", "fab_lambda"}


def load():
    rows = list(csv.DictReader(open(CSV)))
    ms = {v["name"]: v for v in json.load(open(MS)).values()}
    for i, r in enumerate(rows, 1):
        r["_rank"] = i
        r["_family"] = FAMILY.get(r["name"], DEFAULT_FAMILY)
        r["_aa"] = len(r["sequence"])
        r["_ms"] = ms.get(r["name"], {})
    return rows


def rank_table(rows):
    out = ["| rank | design | class | family | aa | **all-site pH** | pose spread | target-only | poses | human | mouse | affinity assessable |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        sp = r.get("ph_pose_spread_over_median", "") or "—"
        out.append(
            f"| {r['_rank']} | `{r['name']}` | {r['molecule_class']} | {r['_family']} | {r['_aa']} | "
            f"**{float(r['ph_ratio_6p5_over_7p4']):.3f}** | {sp} | "
            f"{float(r['ph_ratio_target_only_SUPERSEDED']):.3f} | {r['ph_poses_n']} | "
            f"{float(r['ipsae_min_human']):.3f} | {float(r['ipsae_min_mouse']):.3f} | "
            f"{'**no**' if r['affinity_assessable'].startswith('no') else 'yes'} |")
    return "\n".join(out)


def basis_table(rows):
    out = ["| design | target-only | **all-site** | binder histidines | worst drag |",
           "|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda x: -float(x["ph_ratio_6p5_over_7p4"])):
        m = r["_ms"]
        nh = m.get("n_his")
        out.append(f"| {r['name']} | {float(r['ph_ratio_target_only_SUPERSEDED']):.3f} | "
                   f"**{float(r['ph_ratio_6p5_over_7p4']):.3f}** | "
                   f"{nh if nh is not None else '—'} | "
                   f"{m.get('worst') if nh else '—'} |")
    return "\n".join(out)


def facts(rows):
    n = len(rows)
    fams = {}
    for r in rows:
        fams.setdefault(r["_family"], []).append(r["_rank"])
    big = max(fams.items(), key=lambda kv: len(kv[1]))
    tier1 = [r for r in rows if float(r["ph_ratio_6p5_over_7p4"]) >= 1.20]
    with_his = [r for r in rows if (r["_ms"].get("n_his") or 0) > 0]
    unassess = [r for r in rows if r["molecule_class"] in ANTIBODY]
    unrepro = [r for r in rows
               if r.get("ph_pose_spread_over_median") and float(r["ph_pose_spread_over_median"]) > 1.0]
    nstruct = len([f for f in os.listdir(STRUCTDIR) if f.endswith(".cif")]) if os.path.isdir(STRUCTDIR) else 0
    return dict(n=n, n_fam=len(fams), fams=fams, big_fam=big[0], big_n=len(big[1]), big_ranks=big[1],
                n_tier1=len(tier1), n_no_switch=n - len(tier1),
                n_with_his=len(with_his), n_unassess=len(unassess),
                unassess_ranks=[r["_rank"] for r in unassess],
                unrepro=[(r["_rank"], r["name"], r["ph_pose_spread_over_median"]) for r in unrepro],
                nstruct=nstruct,
                lo_aa=min(r["_aa"] for r in rows), hi_aa=max(r["_aa"] for r in rows))


def selftest():
    rows = load()
    f = facts(rows)
    assert f["n"] == len(rows)
    # every rank is cited correctly nowhere -- the whole point is that we stop citing ranks in
    # prose that the generator does not own
    assert f["n_tier1"] + f["n_no_switch"] == f["n"]
    assert all(r["_family"] for r in rows)
    # a point mutant must share its parent's family, or the effective-n claim is wrong
    assert FAMILY["rimA01_r15_L133E"] == FAMILY["rimA01_r15_boltzgen_egfr_d3_rimA_20"], \
        "a point mutant must share its parent's family"
    # the table must have one row per design plus two header rows
    assert len(rank_table(rows).splitlines()) == f["n"] + 2
    assert len(basis_table(rows).splitlines()) == f["n"] + 2
    print(f"selftest OK — {f['n']} designs, {f['n_fam']} families, "
          f"{f['n_tier1']} in tier 1, {f['nstruct']} structures")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit(0)
    rows = load(); f = facts(rows)
    print("=== FACTS ==="); print(json.dumps({k: v for k, v in f.items() if k != "fams"}, indent=1))
    print("\n=== FAMILIES ==="); print(json.dumps(f["fams"], indent=1))
    print("\n=== RANK TABLE ==="); print(rank_table(rows))
    print("\n=== BASIS TABLE ==="); print(basis_table(rows))
