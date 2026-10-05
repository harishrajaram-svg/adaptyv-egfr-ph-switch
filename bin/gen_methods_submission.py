import json
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
import csv, statistics, json, os, re, sys

CSV = "submissions/01-egfr.csv"
MS = "analysis/01-egfr/ph_sensitivity.json"   # superseded multisite_pooled.json
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


SENS = 'analysis/01-egfr/ph_sensitivity.json'


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
    # ONE SOURCE, JOINED ON SEQUENCE.
    #
    # Two bugs met here. (1) `_ms` was keyed on NAME, and the records carry the internal
    # run name, not the submission name. (2) It read multisite_pooled.json, which
    # emit_submission_csv.py had ALREADY superseded for the CSV (see its multisite()
    # docstring) and which holds 12 records against 18 shipped designs. Together they
    # printed "--" in the published basis table's binder-histidine and worst-drag columns
    # for six of the eighteen rows -- not because the quantity was unmeasured, but because
    # the join missed a stale file.
    #
    # Both columns are now derived from ph_sensitivity.json, the same single source the
    # CSV reads. Verified against the old artifact on all 12 of its records: n_his
    # reproduces exactly, and worst reproduces to three decimals on ten, moving in the
    # third on two where the larger pose set sees a drag n=5 could not.
    ms = {}
    for v in json.load(open(MS)).values():
        q = (v.get("seq") or "").strip().upper()
        poses = v.get("sites_all_poses") or {}
        if not q or not poses:
            continue
        first = next(iter(poses.values()))
        nb = sum(1 for st_ in first.values()
                 if st_["resname"] == "HIS" and st_["partner"] == "binder")
        pw = []
        for pose in poses.values():
            rr = [st_["ratio"] for st_ in pose.values()
                  if st_["resname"] == "HIS" and st_["partner"] == "binder"]
            if rr:
                pw.append(min(rr))
        # No binder histidine means no binder-borne drag to report. The old artifact
        # wrote 1.0 here, which reads as "measured, none found"; it was never measured.
        ms[q] = dict(name=v.get("name"), n_his=nb,
                     worst=round(statistics.median(pw), 3) if pw else None)
    unjoined = [r["name"] for r in rows
                if r["sequence"].strip().upper() not in ms]
    if unjoined:
        raise SystemExit(
            "no pH-sensitivity record for " + ", ".join(unjoined) +
            "\n  The join is on SEQUENCE. Re-run bin/ph_sensitivity_multisite.py so every\n"
            "  shipped design has a record; the basis table must not print a blank cell\n"
            "  for a design whose record merely failed to join.")
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
        r["_ms"] = ms[r["sequence"].strip().upper()]
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
    # COUNT THE BINDER'S OWN HISTIDINES DIRECTLY.
    #
    # This used to compute `n_his - min(n_his)`, treating the global minimum as "the
    # target's contribution". That is invalid, because the target's contribution is NOT
    # constant: designs docked against the domain-III crop see 5 target histidines, but
    # `d2c_mpnn13_S88D_serasp` is scored against the FULL ectodomain and sees 17. So the
    # block published "d2c_mpnn13_S88D_serasp (14)" for a 147 aa design whose sequence
    # contains exactly 2 histidines -- 19 total minus a 5 that did not apply to it.
    #
    # `_ms['n_his']` is the per-design count of binder-partner HIS sites read straight from
    # sites_all_poses, so there is no baseline to get wrong.
    carry = {r["name"]: r["_ms"]["n_his"] for r in rows if r["_ms"].get("n_his")}
    if not carry and rows:
        raise SystemExit("no binder histidine data; re-run bin/ph_sensitivity_multisite.py")
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


def limit_family(rows):
    """Limitation 19's body, generated.

    It was hand-written and duplicated the counts that SS 11.3 already derives, so it
    drifted: it still read "seventeen submitted designs" after the eighteenth was added,
    and "overstates n by up to six-fold" was arithmetic on the old partition. The family
    map, the collapse factor and the totals all come from the CSV, so they are derived.
    SS 11.3 remains the single source of the partition itself; this block only restates
    its consequence."""
    fams = {}
    for r in rows:
        fams.setdefault(r["_family"], []).append(r["_rank"])
    big = max(fams.items(), key=lambda kv: len(kv[1]))
    pairs = sorted((f, v) for f, v in fams.items() if len(v) == 2)
    pd = "; ".join(f"`{f}` (ranks {', '.join(str(x) for x in sorted(v))})"
                  for f, v in pairs)
    return "    " + (f"**Effective n is {len(fams)}, not {len(rows)}.** {len(big[1])} of the "
            f"{len(rows)} submitted designs sit on one backbone (`{big[0]}`, ranks "
            f"{', '.join(str(x) for x in sorted(big[1]))}), and "
            f"{len(pairs)} further families are two-design clusters: {pd}. Any hit rate "
            f"or interval computed over designs rather than sequence families overstates "
            f"n by up to {len(big[1])}-fold on the arm carrying our only causal claim. "
            f"See SS 11.3 for the partition.").replace("SS ", "\u00a7")


def limit_affinity(rows):
    """Limitation 22's body, generated.

    Hand-written as "two of our ten rows ... a fifth of the submission" when the
    submission held ten designs. The denominator changed twice afterwards and the
    fraction did not, so both the count and the share are derived now."""
    bad = [r for r in rows if r["affinity_assessable"].strip().lower()
           not in ("yes", "true", "1")]
    n = len(rows)
    share = (f"{len(bad)} of {n}" if not bad else
             f"{len(bad)} of the {n}")
    frac = f"{100.0 * len(bad) / n:.0f}%"
    names = ", ".join(f"`{r['name']}`" for r in bad)
    return "    " + (f"**The organisers rank outcomes partly on affinity at pH 6.5, and {share} "
            f"submitted rows have no usable affinity reading at all** (SS 4.5): {names}. "
            f"We submitted them anyway, because excluding them would mean scoring them at "
            f"0.0000, which is the error SS 4.2 documents -- but it means {frac} of the "
            f"submission cannot compete on one of the stated criteria."
            ).replace("SS ", "\u00a7").replace(" -- ", " \u2014 ")


def switch_site(rows):
    """Which EGFR histidine each shipped design actually switches on, generated.

    Hand-written as "All ten submitted designs switch on H433". That was stale on the
    count AND wrong on the substance: once the submission grew, one design stopped
    agreeing -- and it is rank 1. A blanket "all" sentence is exactly the shape that hides
    a single disagreeing row, so the distribution is derived and the exception named."""
    import collections
    sens = json.load(open(SENS))
    ship = {r["sequence"].strip().upper(): r for r in rows}
    dom = {}
    for v in sens.values():
        r = ship.get(v.get("seq", "").strip().upper())
        if not r:
            continue
        poses = v.get("sites_all_poses") or {}
        if not poses:
            continue
        best = None
        for k, st_ in next(iter(poses.values())).items():
            if st_["resname"] == "HIS" and st_.get("moved"):
                if best is None or abs(st_["ratio"] - 1) > abs(best[1] - 1):
                    best = (k, st_["ratio"])
        if best:
            dom.setdefault(best[0], []).append((r["name"], r["_rank"]))
    if not dom:
        raise SystemExit("switch_site: no moved histidine on any shipped design")
    order = sorted(dom.items(), key=lambda kv: -len(kv[1]))
    top, rest = order[0], order[1:]
    txt = (f"{len(top[1])} of the {len(rows)} submitted designs switch on "
           f"**{top[0].split(':')[1]}**")
    if rest:
        bits = []
        for site, ms in rest:
            who = ", ".join(f"`{n}` (rank {k})" for n, k in sorted(ms, key=lambda t: t[1]))
            bits.append(f"**{site.split(':')[1]}** -- {who}")
        txt += ". The remaining " + str(sum(len(m) for _, m in rest)) + ": " + "; ".join(bits)
    return txt + "."


# The submission as it stood on 2026-10-04, recovered from git commit d53be33
# ("Fold in the L133 triad..."), the last 10-04 commit before "12 -> 17". Recorded as a
# constant rather than shelled out to git so that the attestation does not depend on the
# repository's history being intact. Verified with:
#   git show d53be33:submissions/01-egfr.csv | tail -n +2 | cut -d, -f1
BASELINE_1004 = [
    'rimA01_r15_boltzgen_egfr_d3_rimA_20', 'bc_s360518_mpnn9_A22D',
    'd2c_mpnn13_S88D_serasp', 'bc_s831683_mpnn6_S15D', 'bc_s831683_mpnn19_S15D',
    'rimA02_d3_rimA_14_vhh', 'h370_020_vhh', 'rimA01_r15_L133E',
    'bc_s831683_mpnn9_S15D', 'bc_s831683_mpnn9_WT', 'bc_d3acid_l65_s831683_mpnn11',
    'bc_s831683_mpnn8_S15D',
]
STRUCTDIR_REL = 'submissions/structures'


def declaration_review(rows):
    """SS 12's human-review attestation, generated.

    It was hand-written and went stale in the worst possible place for a hand-written
    number: an ATTESTATION. It read "the twelve sequences submitted on 2026-10-04" plus
    "five sequences were added on 2026-10-05", naming `sd_d2c_101_l147_s144898_m_T65D` --
    which was added on 10-05 and then SWAPPED OUT the same day -- while omitting designs
    that did ship. So the attestation covered 17 names, one of them not in the submission.
    An attestation that does not match the submission is worse than no attestation.
    """
    names = [r["name"] for r in rows]
    base = [n for n in names if n in BASELINE_1004]
    added = [n for n in names if n not in BASELINE_1004]
    dropped = [n for n in BASELINE_1004 if n not in names]
    out = []
    out.append(f"**Human review.** The submitting researcher has reviewed all "
               f"**{len(names)}** submitted sequences -- their `molecule_class` labels, their "
               f"lengths, and the claims made about them in this document and in the CSV.")
    out.append("")
    out.append(f"Of these, **{len(base)}** were in the submission as it stood on 2026-10-04 "
               f"and **{len(added)}** were added on 2026-10-05. The additions are "
               + ", ".join(f"`{n}`" for n in added) + ".")
    if dropped:
        out.append("")
        out.append("Designs that were in the 2026-10-04 set and are **no longer submitted**: "
                   + ", ".join(f"`{n}`" for n in dropped) + ".")
    out.append("")
    out.append(f"What has been verified for the {len(added)} additions by code, and is "
               f"reproducible from the repository: each comes from this project's own "
               f"generation runs (SS 10); each was re-scored on the same three pH bases over "
               f"its own human-leg poses; and the provenance audit below covers them. What has "
               f"**not** been done for them: expression QC, which was only ever run on the "
               f"original candidate set. Novelty is **not** uniformly established -- see "
               f"`bin/check_novelty_coverage.py`, which is RED.")
    return "\n".join(out).replace("SS ", "\u00a7")


def declaration_structures(rows):
    """SS 12's published-structure coverage, generated.

    Read "Coverage is 10 of 12" naming two missing designs. The directory holds 10 poses
    and the submission holds 18, so 8 are missing. The previous version of the same line
    claimed full coverage of "all ten designs" when the submission held twelve -- the
    identical failure, one revision earlier, which is why it is derived now."""
    d = STRUCTDIR_REL
    files = os.listdir(d) if os.path.isdir(d) else []
    miss = [r["name"] for r in rows
            if not any(r["name"] in f for f in files)]
    have = len(rows) - len(miss)
    txt = (f"**Structures.** Predicted complexes are published at `{d}/`, one median-ipSAE "
           f"pose each -- not the best pose, which would be selection on the outcome. "
           f"**Coverage is {have} of {len(rows)}.**")
    if miss:
        txt += (" Without a published structure: " + ", ".join(f"`{n}`" for n in miss) +
                ". Their poses exist and are scored; they are simply not exported. Stated "
                "rather than implied.")
    return txt


PERT = 'analysis/01-egfr/ph_pka_perturbation.json'


def perturbation_findings(rows):
    """SS 11.7's perturbation findings, generated from the artifact.

    This block was a restatement of a 17-DESIGN run and every claim in it had drifted
    from the 18-design artifact beside it. Worst: it asserted "No design holds a
    top-three position in more than 50% of draws" when `rimA01_r15_L133E` holds one in
    **60%** -- and that design is rank 1, so the sentence was both false and false in the
    submission's own favour. It also named `sd_d2c..._T65D` in its "top set", a design
    that is not submitted, and understated the maximum rank span by a third.
    """
    d = json.load(open(PERT))
    ds = d['designs']
    ship = {r['name'] for r in rows}
    ds = [x for x in ds if x['name'] in ship]
    if len(ds) != len(rows):
        raise SystemExit(f"perturbation artifact covers {len(ds)} of {len(rows)} shipped "
                         f"designs; re-run bin/ph_pka_perturbation.py")
    span = lambda x: x['p95_rank'] - x['p5_rank']
    widest = sorted(ds, key=lambda x: -span(x))[:3]
    best = max(ds, key=lambda x: x['top3_fraction'])
    top = sorted([x for x in ds if x['top3_fraction'] >= 0.25],
                 key=lambda x: -x['top3_fraction'])
    rest = [x for x in ds if x['top3_fraction'] < 0.25]
    floor_ = sorted([x for x in ds if x['top3_fraction'] == 0.0],
                    key=lambda x: x['p5_rank'])
    n = len(rows)
    o = []
    o.append(f"**At PROPKA's own stated accuracy the ordering is not identifiable.** "
             + "; ".join(f"`{x['name'][:34]}` spans ranks {x['p5_rank']}-{x['p95_rank']}"
                         for x in widest)
             + f". The widest span is {span(widest[0])} of {n} ranks.")
    o.append("")
    o.append(f"**The single most-stable design holds a top-three slot in "
             f"{best['top3_fraction'] * 100:.0f}% of draws** (`{best['name']}`, base rank "
             f"{best['base_rank']}, median {best['median_rank']:.0f}). An earlier version of "
             f"this section claimed no design exceeded 50%; it did, and it is the top-ranked "
             f"design, so the error ran in the submission's favour. The conclusion does not "
             f"depend on sigma: it already holds at the optimistic 0.4.")
    o.append("")
    o.append(f"**What does survive.** Two things. First, the **bottom group is robustly at "
             f"the bottom**: "
             + ", ".join(f"`{x['name'][:30]}` never rises above {x['p5_rank']}"
                         for x in floor_[:3])
             + f" -- {len(floor_)} designs take a top-three slot in 0% of draws. "
             f"\"These are not switches\" is stable under the noise. Second, a **top set "
             f"exists even though its order does not**: {len(top)} designs "
             + ", ".join(f"`{x['name'][:30]}` ({x['top3_fraction'] * 100:.0f}%)" for x in top)
             + f" hold a top-three slot in at least 25% of draws, against 0-"
             f"{max((x['top3_fraction'] for x in rest), default=0) * 100:.0f}% for the "
             f"other {len(rest)}.")
    return "\n".join(o)


FOOT = 'analysis/01-egfr/finalist_footprints.json'


def footprint_table(rows):
    """SS 10b's four reviewer-requested footprint checks, generated.

    Hand-written against a 17-design run and wrong on the count that matters: it read
    "1 of 17 touches an N-glycosylation sequon -- and it is the top-ranked design". Three
    of eighteen touch it. Naming only the top-ranked one made the exposure look like a
    single unlucky row rather than a shared property of three designs at ranks 1, 8 and 9.
    """
    d = json.load(open(FOOT))
    ds = [x for x in d['designs'] if x['name'] in {r['name'] for r in rows}]
    if len(ds) != len(rows):
        raise SystemExit(f"footprints cover {len(ds)} of {len(rows)} shipped designs; "
                         f"re-run bin/finalist_footprints.py")
    rank = {r['name']: r['_rank'] for r in rows}
    n = len(ds)
    out = [x for x in ds if x.get('n_outside_d3_crop')]
    gl = [x for x in ds if x.get('glycan_sequon_hits')]
    ids = [x['hu_mo_identity_at_epitope'] for x in ds
           if x.get('hu_mo_identity_at_epitope') is not None]
    lo, hi = d['crop_mature_range']
    rowsout = [
        "| **domain** | **every design, 100% of contacts, in domain III (L2)** "
        "-- no domain-II contact anywhere |",
        f"| **full-ECD** | **{len(out)} of {n}** have any contact outside the "
        f"{hi - lo + 1} aa domain-III crop (mature {lo}-{hi}), so the crop is adequate and "
        f"no footprint required the full ECD to assess |",
    ]
    if gl:
        sites = sorted({f"Asn{p}" for x in gl for p in x['glycan_sequon_hits']})
        who = ", ".join(f"`{x['name'][:36]}` (rank {rank[x['name']]})"
                        for x in sorted(gl, key=lambda y: rank[y['name']]))
        rowsout.append(
            f"| **glycan** | **{len(gl)} of {n}** touch an N-glycosylation sequon, all of "
            f"them {'/'.join(sites)}: {who}. Of {d['n_sequons']} sequons in the construct, "
            f"only {len(sites)} is contacted |")
    else:
        rowsout.append(f"| **glycan** | **0 of {n}** touch an N-glycosylation sequon |")
    rowsout.append(
        f"| **human/mouse** | median identity **at the contacted positions** is "
        f"**{statistics.median(ids):.2f}**; range {min(ids):.2f}-{max(ids):.2f} over "
        f"{len(ids)} designs |")
    return "\n".join(rowsout)


CHAI = 'analysis/01-egfr/chai_interface_summary.json'


def chai_table(rows):
    """Chai-1's verdict on the shipped finalists plus its calibration set, generated.

    The Chai arm ran ELEVEN complexes -- five calibration/reference complexes and SIX
    shipped finalists -- and the document reported one of them (SS 11.2's note on
    `ss_bc_s831683_mpnn6_S15D_S62H_routeA` at 0.838). Reporting one of six finalists from
    an independent predictor, when that predictor disagrees sharply across them, is
    selective in the way SS 13's own arm-accountability limitation warns about. So the
    whole set is tabulated and the comparison against the calibration complexes is stated.
    """
    d = json.load(open(CHAI))
    cx = d['complexes']
    byname = {r['name']: r for r in rows}

    def match(tag):
        t = tag[4:] if tag.startswith('fin_') else tag
        best = None
        for n in byname:
            if n.startswith(t) or t.startswith(n) or n.replace('__', '_').startswith(t):
                if best is None or len(n) > len(best):
                    best = n
        return best

    fin, ref = [], []
    for c in cx:
        n = match(c['tag']) if c['tag'].startswith('fin_') else None
        (fin if n else ref).append((n, c))
    fin.sort(key=lambda t: -t[1]['iptm_median'])
    ref.sort(key=lambda t: -t[1]['iptm_median'])
    out = ["| complex | Chai-1 ipTM (median of 5) | interface residues | clashing models |",
           "|---|---|---|---|"]
    for n, c in fin:
        out.append(f"| **`{n}`** (rank {byname[n]['_rank']}) | **{c['iptm_median']:.3f}** | "
                   f"{c['iface_residues_median']} | {c['clash_models']} |")
    out.append("| *— calibration and reference complexes —* | | | |")
    for _, c in ref:
        out.append(f"| `{c['tag']}` | {c['iptm_median']:.3f} | "
                   f"{c['iface_residues_median']} | {c['clash_models']} |")
    hi = fin[0]
    lo = fin[-1]
    egf = next((c for _, c in ref if c['tag'].startswith('egf')), None)
    txt = "\n".join(out)
    txt += (f"\n\nSix of the {len(rows)} shipped designs were folded by Chai-1, and it does "
            f"**not** rate them alike: `{hi[0]}` reads {hi[1]['iptm_median']:.3f} against "
            f"`{lo[0]}` at {lo[1]['iptm_median']:.3f}, a spread of "
            f"{hi[1]['iptm_median'] - lo[1]['iptm_median']:.3f} ipTM across designs our own "
            f"pH objective orders quite differently.")
    if egf:
        below = [n for n, c in fin if c['iptm_median'] < egf['iptm_median']]
        txt += (f" Three of the six sit at or above the cetuximab scFv positive control "
                f"(0.793); {len(below)} sit **below human EGF** "
                f"({egf['iptm_median']:.3f}): " + ", ".join(f"`{n}`" for n in below) + ".")
    txt += (" Chai emits no residue-level PAE, so ipSAE cannot be computed on these and "
            "ipTM is not comparable to our ranking metric. It is a second opinion on whether "
            "an interface forms at all, not a second measurement of the objective.")
    g532 = next((c for _, c in ref if c['tag'].startswith('g532')), None)
    nano = next((c for _, c in ref if c['tag'].startswith('nano')), None)
    if g532:
        txt += (f"\n\n**And the calibration set says not to over-read it.** `g532_ecd` is a "
                f"PUBLISHED, experimentally-confirmed pH-switchable EGFR binder, and Chai-1 "
                f"scores it **{g532['iptm_median']:.3f}** -- below three of our six designs "
                f"and well below human EGF. ")
        if nano:
            txt += (f"`nano2_ecd` reads {nano['iptm_median']:.3f} on "
                    f"{nano['iface_residues_median']} interface residues, the largest "
                    f"interface in the set and the lowest score. ")
        txt += ("So a low Chai ipTM is **not** evidence that a design does not bind: on the "
                "one molecule here with a real measured answer, this metric is wrong. The "
                "table supports the positive direction only -- three designs form an "
                "interface an independent predictor rates at the level of the cetuximab "
                "control -- and it cannot be used to argue against the designs at the bottom, "
                "including rank 1. Reporting it the other way round would be the "
                "single most tempting over-read available in this submission.")
    return txt


def sigma_table(rows):
    """SS 11.7's sigma sweep, generated from the three persisted runs.

    The table was a 17-design run ("4 of 17", "16 of 17") and limitation 17 cited it for a
    figure it did not contain. The three sigmas now live in
    analysis/01-egfr/ph_pka_perturbation_sigma{0.4,0.8,1.2}.json so the sweep is an
    artifact rather than a memory of three separate invocations."""
    import glob as _g
    out = ["| sigma (pKa units) | keep baseline rank | span >= 5 ranks |", "|---|---|---|"]
    found = 0
    for sg in ('0.4', '0.8', '1.2'):
        f = f'analysis/01-egfr/ph_pka_perturbation_sigma{sg}.json'
        if not os.path.exists(f):
            continue
        d = json.load(open(f))
        ds = [x for x in d['designs'] if x['name'] in {r['name'] for r in rows}]
        if len(ds) != len(rows):
            raise SystemExit(f"{f} covers {len(ds)} of {len(rows)} shipped designs")
        keep = sum(1 for x in ds if x['median_rank'] == x['base_rank'])
        wide = sum(1 for x in ds if x['p95_rank'] - x['p5_rank'] >= 5)
        label = (f"**{sg} (PROPKA's own RMSD)**" if sg == '0.8' else sg)
        bold = '**' if sg == '0.8' else ''
        out.append(f"| {label} | {bold}{keep} of {len(ds)}{bold} | "
                   f"{bold}{wide} of {len(ds)}{bold} |")
        found += 1
    if found < 3:
        raise SystemExit("sigma sweep incomplete; re-run bin/ph_pka_perturbation.py at "
                         "--sigma 0.4, 0.8 and 1.2 and persist each")
    return "\n".join(out)


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
              "FAMILY-LIST": family_list(rows), "BINDER-HIS": binder_his_sentence(rows),
              "LIMIT-FAMILY": limit_family(rows),
              "LIMIT-AFFINITY": limit_affinity(rows),
              "SWITCH-SITE": switch_site(rows),
              "DECL-REVIEW": declaration_review(rows),
              "DECL-STRUCT": declaration_structures(rows),
              "PERT-FINDINGS": perturbation_findings(rows),
              "FOOTPRINT-TABLE": footprint_table(rows),
              "CHAI-TABLE": chai_table(rows),
              "SIGMA-TABLE": sigma_table(rows)}
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
                      "BINDER-HIS": binder_his_sentence(rows),
                      "LIMIT-FAMILY": limit_family(rows),
              "LIMIT-AFFINITY": limit_affinity(rows),
              "SWITCH-SITE": switch_site(rows),
              "DECL-REVIEW": declaration_review(rows),
              "DECL-STRUCT": declaration_structures(rows),
              "PERT-FINDINGS": perturbation_findings(rows),
              "FOOTPRINT-TABLE": footprint_table(rows),
              "CHAI-TABLE": chai_table(rows),
              "SIGMA-TABLE": sigma_table(rows)}.items():
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
