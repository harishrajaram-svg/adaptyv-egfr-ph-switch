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
from pathlib import Path
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
    # ADDED 2026-10-05 from the reopened exclusion pool (METHODS 10). These must be listed
    # explicitly: DEFAULT_FAMILY is the s831683 backbone, so an unmapped name was silently
    # reported as belonging to the family that already holds six designs -- which would
    # have made the submission look MORE concentrated than it is and misstated the one
    # thing these five were selected to improve.
    "c5_cf_short__boltzgen_egfr_cropfree_short_48":  "cf_cropfree_short (c5)",
    "c5_cr_crop_patch__boltzgen_egfr_crop_patch_05": "cr_crop_patch (c5)",
    "cons_gap_h370_only__boltzgen_egfr_h370_018":    "h370_018 (gap)",
    # SWAPPED IN 2026-10-05, replacing sd_d2c_..._m_T65D and
    # ss_bc_s831683_mpnn6_S15D_S62H_routeA. Those two were 0.986 and 0.985 identical to
    # designs already shipped -- "nearly identical variants", which the review explicitly
    # said not to fill slots with. These two are at most 0.467 identical to anything else
    # submitted and open a backbone family that had no representation.
    "bcr_d3acid3_l60_s647537_mpnn3":                 "d3acid3_l60_s647537",
    "bcr_d3acid3_l60_s647537_mpnn11":                "d3acid3_l60_s647537",
    "ss_bc_s831683_mpnn6_S15D_S62H_routeA":          "d3acid_l65_s831683",
    # The six that previously relied on DEFAULT_FAMILY. Listing them is the point:
    # the default was correct for exactly these and silently wrong for anything new.
    "bc_s831683_mpnn6_S15D":                           "d3acid_l65_s831683",
    "bc_s831683_mpnn19_S15D":                          "d3acid_l65_s831683",
    "bc_s831683_mpnn9_S15D":                           "d3acid_l65_s831683",
    "bc_s831683_mpnn9_WT":                             "d3acid_l65_s831683",
    "bc_d3acid_l65_s831683_mpnn11":                    "d3acid_l65_s831683",
    "bc_s831683_mpnn8_S15D":                           "d3acid_l65_s831683",
}
# No silent default. An unmapped design used to be reported as the s831683 backbone, so a
# newly added design would quietly inherit the wrong family in the graded methods table.
DEFAULT_FAMILY = None
ANTIBODY = {"nanobody", "scfv", "fab_kappa", "fab_lambda"}


HEAD_PREFIX = "ph_ratio_6p5_over_7p4"


def headcol(rows):
    """The headline pH column, resolved by PREFIX rather than by exact name.

    The column was renamed to ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE when the
    three-basis sensitivity analysis landed, because the old name implied a single
    settled quantity. This generator exists so the methods document cannot drift from
    the CSV, so it resolves the column instead of hardcoding a spelling -- and raises
    if it cannot, rather than falling back to something that looks similar."""
    cands = [c for c in rows[0] if c.startswith(HEAD_PREFIX)]
    if len(cands) != 1:
        raise SystemExit(f"expected exactly one column starting {HEAD_PREFIX!r}, "
                         f"found {cands!r} in {list(rows[0])!r}")
    return cands[0]


def load():
    rows = list(csv.DictReader(open(CSV)))
    ms = {v["name"]: v for v in json.load(open(MS)).values()}
    for i, r in enumerate(rows, 1):
        r["_rank"] = i
        fam = FAMILY.get(r["name"], DEFAULT_FAMILY)
        if fam is None:
            raise SystemExit(
                f"{r['name']} has no FAMILY entry in bin/gen_methods_submission.py.\n"
                "  Add it. There is deliberately no default: the old one silently assigned\n"
                "  the s831683 backbone, so an unmapped design misreported the submission's\n"
                "  family concentration in the graded methods table.")
        r["_family"] = fam
        r["_aa"] = len(r["sequence"])
        r["_ms"] = ms.get(r["name"], {})
    return rows


def rank_table(rows):
    out = ["| rank | design | class | family | aa | **pH his-only (ranked)** | all-site | partnered | rank range | pose spread | target-only | poses | human | mouse | affinity assessable |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        sp = r.get("ph_pose_spread_over_median", "") or "—"
        out.append(
            f"| {r['_rank']} | `{r['name']}` | {r['molecule_class']} | {r['_family']} | {r['_aa']} | "
            f"**{float(r[headcol(rows)]):.3f}** | "
            f"{r.get('ph_ratio_allsite_SENSITIVITY') or '—'} | "
            f"{r.get('ph_ratio_partnered_SENSITIVITY') or '—'} | "
            f"{r.get('ph_rank_range_across_bases') or '—'} | {sp} | "
            f"{float(r['ph_ratio_target_only_SUPERSEDED']):.3f} | {r['ph_poses_n']} | "
            f"{float(r['ipsae_min_human']):.3f} | {float(r['ipsae_min_mouse']):.3f} | "
            f"{'**no**' if r['affinity_assessable'].startswith('no') else 'yes'} |")
    return "\n".join(out)


def binder_his_sentence(rows):
    """The binder-histidine census of SS 11.1, generated.

    This count has been wrong FIVE times by hand: "six of the eleven", "seven of twelve",
    "eight of the seventeen" (a pre-addition count carried forward), "ten of seventeen"
    (correct for five additions, two of which were then swapped out), and "eight of the
    seventeen" again after the swap. It is derived here and never typed.

    The target construct contributes its own histidines to n_his, so a binder carries
    histidines when n_his exceeds the target's count; that baseline is read from the
    designs that have none rather than hard-coded.
    """
    import json as _j
    sens = _j.load(open('analysis/01-egfr/ph_sensitivity.json'))
    ship = {r["sequence"].strip().upper(): r["name"] for r in rows}
    per = {}
    for v in sens.values():
        q = v.get("seq", "").strip().upper()
        if q in ship and "n_his" in v:
            per[ship[q]] = v["n_his"]
    if not per:
        raise SystemExit("no n_his data in ph_sensitivity.json; re-run ph_sensitivity_multisite.py")
    base = min(per.values())                      # target-only histidine count
    carry = {k: v - base for k, v in per.items() if v > base}
    WORDS = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
             7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve",
             13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen",
             17: "seventeen", 18: "eighteen", 19: "nineteen", 20: "twenty"}
    n, tot = len(carry), len(rows)
    groups = {}
    for k, v in sorted(carry.items(), key=lambda kv: (-kv[1], kv[0])):
        groups.setdefault(v, []).append(k)
    parts = []
    for cnt, names in sorted(groups.items(), reverse=True):
        lst = ", ".join(f"`{x}`" for x in names)
        parts.append(f"{lst} ({cnt} each)" if len(names) > 1 else f"{lst} ({cnt})")
    return (f"**{WORDS.get(n, n)} of the {WORDS.get(tot, tot)} submitted designs carry at "
            f"least one histidine of their own**: " + "; ".join(parts) +
            f". The other {WORDS.get(tot - n, tot - n)} carry none.")


def family_list(rows):
    """The family breakdown of SS 11.3, generated.

    This paragraph was hand-maintained and went stale the moment the submission changed:
    after five designs were added it still read "Six families, twelve designs" and listed
    ranks from the previous build. The family map, the counts and the ranks are all
    derivable from the CSV, so they are derived."""
    from collections import defaultdict
    fams = defaultdict(list)
    for r in rows:
        fams[r["_family"]].append(r["_rank"])
    out = []
    for fam, ranks in sorted(fams.items(), key=lambda kv: (-len(kv[1]), min(kv[1]))):
        rs = ", ".join(str(x) for x in sorted(ranks))
        if len(ranks) > 1:
            out.append(f"`{fam}` **x{len(ranks)}** (ranks {rs})")
        else:
            out.append(f"`{fam}` (rank {rs})")
    body = " - ".join(out)
    big = max(fams.items(), key=lambda kv: len(kv[1]))
    singles = sum(1 for v in fams.values() if len(v) == 1)
    return (body + "\n\n**Effective n is " + str(len(fams)) + " clusters, not " +
            str(len(rows)) + " designs.** The largest cluster, `" + big[0] + "`, holds " +
            str(len(big[1])) + " designs at ranks " +
            ", ".join(str(x) for x in sorted(big[1])) + "; " + str(singles) +
            " families contribute a single design each. Any interval must be computed on "
            "families, not designs.")


def basis_table(rows):
    out = ["| design | target-only | **all-site** | binder histidines | worst drag |",
           "|---|---|---|---|---|"]
    HC = headcol(rows)
    for r in sorted(rows, key=lambda x: -float(x[HC])):
        m = r["_ms"]
        nh = m.get("n_his")
        out.append(f"| {r['name']} | {float(r['ph_ratio_target_only_SUPERSEDED']):.3f} | "
                   f"**{float(r[HC]):.3f}** | "
                   f"{nh if nh is not None else '—'} | "
                   f"{m.get('worst') if nh else '—'} |")
    return "\n".join(out)


def facts(rows):
    n = len(rows)
    fams = {}
    for r in rows:
        fams.setdefault(r["_family"], []).append(r["_rank"])
    big = max(fams.items(), key=lambda kv: len(kv[1]))
    tier1 = [r for r in rows if float(r[headcol(rows)]) >= 1.20]
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


METHODS = "submissions/01-egfr-METHODS.md"


def write_into_methods(rows):
    """Replace the marked generated blocks in METHODS, in place.

    Before this, the generator PRINTED the tables and a human pasted them. That is the
    drift channel it was built to close -- §11 had already gone stale once that way,
    describing eleven designs with a rank table missing the twelfth. Now the tables
    cannot disagree with the CSV unless someone edits between the markers, and
    --check reports it if they do."""
    import re
    src = Path(METHODS).read_text()
    blocks = {"BASIS-TABLE": basis_table(rows), "RANK-TABLE": rank_table(rows),
              "FAMILY-LIST": family_list(rows), "BINDER-HIS": binder_his_sentence(rows)}
    for tag, body in blocks.items():
        pat = re.compile(rf"(<!-- GENERATED:{tag}[^>]*-->\n).*?(<!-- /GENERATED:{tag} -->)",
                         re.S)
        if not pat.search(src):
            raise SystemExit(f"no GENERATED:{tag} block found in {METHODS}")
        src = pat.sub(lambda m: m.group(1) + body + "\n" + m.group(2), src)
    Path(METHODS).write_text(src)
    print(f"wrote {len(blocks)} generated blocks into {METHODS}")


def check_methods(rows):
    """Exit non-zero if a generated block in METHODS differs from what the CSV implies."""
    import re
    src = Path(METHODS).read_text()
    bad = []
    for tag, body in {"BASIS-TABLE": basis_table(rows), "RANK-TABLE": rank_table(rows),
                      "FAMILY-LIST": family_list(rows),
                      "BINDER-HIS": binder_his_sentence(rows)}.items():
        m = re.search(rf"<!-- GENERATED:{tag}[^>]*-->\n(.*?)<!-- /GENERATED:{tag} -->", src, re.S)
        if not m:
            bad.append(f"{tag}: block missing")
        elif m.group(1).strip() != body.strip():
            bad.append(f"{tag}: differs from the CSV")
    if bad:
        print("METHODS is STALE: " + "; ".join(bad))
        print("  fix with: python3 bin/gen_methods_submission.py --write")
        sys.exit(1)
    print("METHODS generated blocks match the CSV")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit(0)
    if "--write" in sys.argv:
        write_into_methods(load()); sys.exit(0)
    if "--check" in sys.argv:
        check_methods(load()); sys.exit(0)
    rows = load(); f = facts(rows)
    print("=== FACTS ==="); print(json.dumps({k: v for k, v in f.items() if k != "fams"}, indent=1))
    print("\n=== FAMILIES ==="); print(json.dumps(f["fams"], indent=1))
    print("\n=== RANK TABLE ==="); print(rank_table(rows))
    print("\n=== BASIS TABLE ==="); print(basis_table(rows))
    print("\n=== FAMILY LIST ==="); print(family_list(rows))
    print("\n=== BINDER-HIS ==="); print(binder_his_sentence(rows))
