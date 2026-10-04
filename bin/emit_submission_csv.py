#!/usr/bin/env python3
"""Emit the Proteinbase submission CSV, ranked by the challenge's own stated hierarchy.

COMPLIANCE, from the challenge page and FAQ (read 2026-10-04):
  * Track 3 allows AT MOST 20 designs per challenge. We were at 21.
  * Required columns: `name`, `sequence`, `molecule_class`
    (one of protein / nanobody / scfv / fab_kappa / fab_lambda).
    Our JSON carried `seq`, not `sequence`, and no molecule_class at all.
  * "Submit your designs as a CSV ordered by how you would rank your molecules
    (top row higher)."
  * Single chains must be 10-250 aa.
  * Extra metrics are explicitly encouraged and are fed to the selector, so the
    optional columns are included rather than stripped.

RANKING = the challenge page's stated priority order, not our old iptm order:
    (1) pH-selective binding  (2) mouse cross-reactivity  (3) human affinity
  tier 1: designs with a real switch (ratio >= RATIO_BAR), by ratio descending
  tier 2: everything else, by MOUSE then HUMAN -- because the page ranks mouse
          cross-reactivity ABOVE affinity to human EGFR, which our previous
          iptm-descending order did not reflect.

THE CUT. Three designs score zero on all three objectives. We drop the weakest to
reach 20. Note PK's advice -- "I would not fill the allocation simply to reach
twenty" -- so the data supports cutting further; that is a judgment call left open.

Usage:
    python3 bin/emit_submission_csv.py            # writes submissions/01-egfr.csv
    python3 bin/emit_submission_csv.py --selftest
"""
import csv, json, os, sys, importlib.util

SUB   = "analysis/01-egfr/submission_final.json"
OUT   = "submissions/01-egfr.csv"
LIMIT = 20
RATIO_BAR = 1.20
MIN_AA, MAX_AA = 10, 250
MOLECULE_CLASS = "protein"     # DEFAULT only -- per-design `molecule_class` overrides it.
# Was hardcoded for every row. That is a compliance error the moment a non-protein format
# enters the file: `rimA02_d3_rimA_14_vhh` is a VHH and must ship as `nanobody`, both because
# the label must be true and because Adaptyv score ANTIBODY novelty by a different rule
# (CDRH3 < 70% AND global >= 70% = Level 3). Under the general-protein rule that design reads
# 77.5% identity and looks rejected; under the correct rule it is Level 3 and eligible.
VALID_CLASSES = ("protein", "nanobody", "scfv", "fab_kappa", "fab_lambda")


def load_scores():
    spec = importlib.util.spec_from_file_location("ms", "bin/make_submission.py")
    ms = importlib.util.module_from_spec(spec); spec.loader.exec_module(ms)
    return ms, ms.build_index()


MIN_N = 5          # poses required before a ratio may put a design in tier 1

# AFFINITY PRECONDITION — RETIRED AS A HARD GATE on PK's review, 2026-10-04.
#
# PK: "I would not replace the compromised bar with another universal number ... If adequate
# calibration is infeasible within the existing budget, retain continuous scores and flag
# uncertainty rather than inventing a hard gate."
#
# The 0.1493 / 0.2218 bars came from `NEG_nonbinder`, which is 81% mature human EGF. PK's
# correction to our own reading of that: the sequence similarity is enough to STOP treating it
# as an established nonbinder, but it does NOT establish that the altered sequence retains
# EGF's affinity or agonist activity, and the higher domain-III score is not proof of binding
# either — changing the target crop changes prediction behaviour. Both EGF-derived controls are
# therefore **activity-unknown**, not "binders". Claims resting on their negative status are
# withdrawn, and the original numbers are preserved for the audit trail.
#
# So: no design is EXCLUDED on affinity. Affinity becomes a reported, continuous column with an
# explicit uncertainty flag, and the tier-1 test is the pH hypothesis plus pose count only.
MIN_AFFINITY = None        # retired; see above. Set a float only to re-enable a hard gate.
AFFINITY_FLAG = 0.2218     # REPORTING threshold: "above the (compromised) computational null"


def rank_key(r):
    """Challenge hierarchy: pH tier first by ratio, then mouse, then human.

    TIER 1 NOW REQUIRES n >= MIN_N POSES (decided 2026-10-04 with Harish). The ratios in
    this file are the UNIFORM refold pH -- one code path, every submitted design, median
    over every ESMFold2 human-leg pose of that exact binder sequence, pooled across runs
    (n ranges 5 to 31). They replace a column that mixed provenances: several entries were
    never-gated placeholder zeros, and three were measured by older per-design scripts that
    disagreed with the uniform value by up to 1.8x.
        rimA01_r02  2.151 -> **0.940** (n=31)   falls out of tier 1
        bg04_r05    1.786 -> **0.889** (n=11)   falls out of tier 1
        bg02_r01    1.041 -> **0.718** (n=10)   was already tier 2
        bg04_r03    2.472 ->   2.499  (n=16)   confirmed
    The n floor exists because a single pose is not a measurement -- `rimA01/d3_rimA_14`
    read 0.2097/0.4255 on one seed and 0.0000/0.0000 on five.
    """
    n = r.get("ratio_n", 0)
    # Affinity no longer EXCLUDES (PK, 2026-10-04): it orders within tier 2 and is reported with
    # a flag. A hard gate here would launder a compromised control into a yes/no decision.
    binds = True if MIN_AFFINITY is None else max(r["hu"], r["mo"]) >= MIN_AFFINITY
    tier1 = r["ratio"] >= RATIO_BAR and n >= MIN_N and binds
    return (0 if tier1 else 1,
            -r["ratio"] if tier1 else 0,
            -r["mo"], -r["hu"])


def main():
    rows = json.load(open(SUB))
    ms, idx = load_scores()
    scored = []
    for x in rows:
        hu, _, _, _ = ms.score(idx, x["seq"], "hu", "v2")
        mo, _, _, _ = ms.score(idx, x["seq"], "mo", "v2")
        scored.append(dict(name=x["name"], sequence=x["seq"],
                           molecule_class=x.get("molecule_class", MOLECULE_CLASS),
                           ratio=float(x["ratio"]), hu=hu or 0.0, mo=mo or 0.0,
                           # ratio_n must be carried through: rank_key fails CLOSED without
                           # it, so dropping it here silently emptied tier 1 and put a
                           # non-switching binder at rank 1.
                           ratio_n=int(x.get("ratio_n", 0)),
                           # per-design assessment text must be carried through, or PK's
                           # uncertainty column silently collapses to one default string.
                           # Same class of bug as ratio_n, which emptied tier 1 earlier today.
                           assessment=x.get("assessment", "computational candidate"),
                           iptm_generator=x.get("iptm", ""), rmsd=x.get("rmsd", "")))

    for r in scored:
        n = len(r["sequence"])
        assert MIN_AA <= n <= MAX_AA, f"{r['name']}: {n} aa is outside {MIN_AA}-{MAX_AA}"
    assert len({r["sequence"] for r in scored}) == len(scored), "duplicate sequences"
    assert len({r["name"] for r in scored}) == len(scored), "duplicate names"
    for r in scored:
        assert r["molecule_class"] in VALID_CLASSES, f"{r['name']}: bad molecule_class {r['molecule_class']}"

    scored.sort(key=rank_key)
    dropped = []
    while len(scored) > LIMIT:
        dropped.append(scored.pop())            # the tail is the weakest under the hierarchy
    if dropped:
        print(f"dropped {len(dropped)} to meet the {LIMIT}-design limit:")
        for d in dropped:
            print(f"   ratio {d['ratio']:.3f}  mouse {d['mo']:.4f}  human {d['hu']:.4f}  {d['name']}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    # `refold_rmsd` WAS SHIPPED AND IS NOW DROPPED (2026-10-04). It mixed two
    # different measurements under one name:
    #   * 18 rows carried BoltzGen's own `filter_rmsd` -- its designed backbone against
    #     its OWN refold. Self-consistency. Range across the submission: 0.78-2.0 A.
    #   * 2 rows carried BoltzGen pose against an INDEPENDENT ESMFold2 refold.
    #     Cross-predictor agreement. Recomputed consistently where both poses exist:
    #       rank  1 rimA01_r15    shipped 4.99  ->  8.23 A
    #       rank  3 bg04_r03      shipped 0.96  ->  4.15 A
    #       rank  4 rimA01_r02    shipped 1.80  ->  2.66 A
    #       rank 16 bg01_r02      shipped 0.89  -> 23.26 A
    #       rank 17 mechA01_r01   shipped 1.47  -> 10.38 A
    #       rank  2 short_031     1.88 (already cross-predictor, most consistent design)
    # The shipped column therefore understated pose disagreement for most rows and
    # made the two honestly-measured designs look WORSE than the rest. The four
    # Mosaic designs have no generator pose at all, so cross-predictor RMSD is
    # undefined for them and a uniform column is not available before the deadline.
    # The column is optional metadata; shipping one name for two quantities is worse
    # than shipping neither, so it is dropped and the finding goes in the methods doc.
    cols = ["name", "sequence", "molecule_class",
            "ph_ratio_6p5_over_7p4", "ipsae_min_human", "ipsae_min_mouse",
            "iptm_generator", "ph_poses_n", "affinity_above_null", "assessment"]
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in scored:
            w.writerow([r["name"], r["sequence"], r["molecule_class"],
                        f"{r['ratio']:.3f}", f"{r['hu']:.4f}", f"{r['mo']:.4f}",
                        r["iptm_generator"], r.get("ratio_n", 0),
                        # PK: separate columns for binding evidence, pH hypothesis and uncertainty
                        "yes" if max(r["hu"], r["mo"]) >= AFFINITY_FLAG else "no",
                        r.get("assessment", "computational candidate")])
    print(f"\nwrote {OUT}: {len(scored)} designs, {len(cols)} columns")
    print(f"{'#':>3} {'ratio':>7} {'mouse':>7} {'human':>7}  name")
    for i, r in enumerate(scored, 1):
        print(f"{i:>3} {r['ratio']:>7.3f} {r['mo']:>7.4f} {r['hu']:>7.4f}  {r['name'][:44]}")


def selftest():
    # NOTE: `a` used to carry mo=hu=0.1 and still counted as tier 1. Under the affinity
    # precondition it no longer does, which is the point of the precondition -- a 4.5x ratio on
    # a design scoring 0.1 is not a ratio of affinities. The fixture is raised above the bar so
    # it tests tier ordering rather than the gate it was silently exempt from.
    a = dict(ratio=4.5, mo=0.3, hu=0.3, ratio_n=5); b = dict(ratio=0.0, mo=0.9, hu=0.9, ratio_n=5)
    assert rank_key(a) < rank_key(b), "a real switch must outrank a non-switching binder"
    c = dict(ratio=0.0, mo=0.5, hu=0.9, ratio_n=5); d = dict(ratio=0.0, mo=0.9, hu=0.1, ratio_n=5)
    assert rank_key(d) < rank_key(c), "mouse must outrank human affinity in tier 2"
    # Every tier-1 fixture below must also satisfy the AFFINITY precondition, or it is testing
    # the wrong thing. BIND = just over the bar; the three gates are then varied one at a time.
    BIND = (MIN_AFFINITY if MIN_AFFINITY is not None else AFFINITY_FLAG) + 0.01
    # the n floor: a big ratio on too few poses must NOT reach tier 1
    thin = dict(ratio=9.9, mo=BIND, hu=BIND, ratio_n=1)
    solid = dict(ratio=1.3, mo=BIND, hu=BIND, ratio_n=5)
    assert rank_key(solid) < rank_key(thin), "n<MIN_N must be excluded from tier 1"
    # a missing n must fail closed, not silently qualify
    assert rank_key(solid) < rank_key(dict(ratio=9.9, mo=BIND, hu=BIND)), "absent n must fail closed"
    # the affinity precondition: a huge ratio on a design that does not bind is NOT tier 1
    nonbinder = dict(ratio=9.9, mo=0.0, hu=0.0, ratio_n=9)
    if MIN_AFFINITY is None:
        # gate retired: a 0.0000/0.0000 design with a confirmed switch MAY hold tier 1, because
        # we no longer have a control that licenses calling it a nonbinder. It is reported with
        # an uncertainty flag instead of being excluded.
        assert rank_key(nonbinder) < rank_key(solid), "with the gate retired, ratio orders tier 1"
    else:
        assert rank_key(solid) < rank_key(nonbinder), "a 0.0000/0.0000 design must not reach tier 1"
    assert MOLECULE_CLASS in VALID_CLASSES
    print("selftest OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit(0)
    selftest(); main()
