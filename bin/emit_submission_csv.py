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

THE CUT -- settled 2026-10-04 with Harish, on PK's advice ("I would not fill the
allocation simply to reach twenty"). LIMIT is 10, not the allowed 20.

Ranks 11-20 of the 20-design build did not stand on a measurement. Eight of them read
BELOW 1.0x -- no switch at all -- sitting on the 0.702x steric floor that every design
touching a histidine without a nearby carboxylate returns. That floor is a constant of
the method, so those rows reported the method back to itself. The weakest, bg04_r03,
was 1.94x pooled with 0.0000/0.0000 on both species and three expression-QC flags
including an unpaired cysteine.

All 10 that remain stand on something measured: the four S15D designs, their matched
wild-type control, rimA01_r15, d2c S88D, mpnn11, and the two VHH-format designs that
section 4.5 of the methods document shows this instrument cannot score (reported as
such, not as scores).

SIDE EFFECT, recorded because it was an open question: `rank_key` sorts BOTH tiers by
(-mo, -hu) while the recorded decision says tier 2 should sort by HUMAN, since mouse is
measured at pH 7.4 only and a design that switches off at 7.4 has no defined mouse
ratio. That contradiction only ever changed which null-value rows padded the tail. At
LIMIT = 10 no tier-2 design ships, so it is moot and is left unchanged rather than
edited blind at submission time.

Usage:
    python3 bin/emit_submission_csv.py            # writes submissions/01-egfr.csv
    python3 bin/emit_submission_csv.py --selftest
"""
import csv, json, os, sys, importlib.util

SUB   = "analysis/01-egfr/submission_final.json"
OUT   = "submissions/01-egfr.csv"
LIMIT = 11        # 10 cut designs + A22D added 2026-10-04. NOT the allowed 20 -- see THE CUT.
RATIO_BAR = 1.20
MIN_AA, MAX_AA = 10, 250
MOLECULE_CLASS = "protein"     # DEFAULT only -- per-design `molecule_class` overrides it.
# Was hardcoded for every row. That is a compliance error the moment a non-protein format
# enters the file: `rimA02_d3_rimA_14_vhh` is a VHH and must ship as `nanobody`, both because
# the label must be true and because Adaptyv score ANTIBODY novelty by a different rule
# (CDRH3 < 70% AND global >= 70% = Level 3). Under the general-protein rule that design reads
# 77.5% identity and looks rejected; under the correct rule it is Level 3 and eligible.
VALID_CLASSES = ("protein", "nanobody", "scfv", "fab_kappa", "fab_lambda")
# Formats whose affinity this instrument cannot read -- see METHODS 4.5 and the G532 inversion.
ANTIBODY_CLASSES = {"nanobody", "scfv", "fab_kappa", "fab_lambda"}


def load_scores():
    spec = importlib.util.spec_from_file_location("ms", "bin/make_submission.py")
    ms = importlib.util.module_from_spec(spec); spec.loader.exec_module(ms)
    return ms, ms.build_index()


# AFFINITY NOW COMES FROM master_rank.json, NOT make_submission.score() -- fixed 2026-10-04.
#
# WHY. make_submission.score() calls pick_run(), which selects ONE run ("most seeds wins,
# ties broken lexicographically") and takes the median over that run's seeds. For a molecule
# that appears in several run directories -- and the same molecule appears here under up to
# three run names -- the reported affinity is therefore a function of WHICH RUN NAME sorted
# first, which is the keystone trap of this project stated verbatim in HANDOFF: "That file is
# the single source for both columns; do not re-derive either one from a run name."
#
# The emitter was re-deriving from a run name. Measured consequence: 5 of the 10 shipped
# affinity cells matched neither the pooled median nor the pooled max, because they were one
# run's median. All 5 were molecules with poses in more than one run (ratio_n 6, 6, 6, 11, 11);
# the five single-run molecules agreed exactly. METHODS 11 meanwhile told the reader the
# columns were medians, so the graded file disagreed with its own methods document.
#
# master_rank.py pools every pose of a binder SEQUENCE across all runs and is the file both
# HANDOFF and METHODS name as canonical. Reading it here makes the CSV, master_rank.json and
# METHODS agree, and removes the last run-name-keyed join on the submission path.
MASTER = "analysis/01-egfr/master_rank.json"
MULTISITE = "analysis/01-egfr/multisite_pooled.json"

# RANKING BASIS CHANGED 2026-10-04 19:50 EDT, with Harish, from the target-only pH ratio to the
# ALL-TITRATABLE-SITE product over BOTH partners.
#
# WHY. The competition's primary objective is KD(7.4)/KD(6.5). Until tonight we estimated it with
# bin/ph_gate_refolds.py, which measures only the TARGET's histidines. It never measured our own
# binders' titratable groups, and six of the then-ten submitted designs carry two or three
# histidines of their own. Those get buried at the interface, lose 1.5-2.5 pKa units, and by
# thermodynamic linkage that OPPOSES acid-tightening. Composing honestly over every titratable
# site on both partners (bin/ph_gate_multisite.py, 65 poses, n=5-11 per design):
#
#   rimA02_d3_rimA_14_vhh   5.186 -> 4.838     no binder histidine
#   rimA01_r15              4.582 -> 4.256     no binder histidine
#   bc_s360518_mpnn9_A22D   5.630 -> 3.738     one, dragging 0.776
#   d2c_mpnn13_S88D         4.572 -> 3.526     two, worst 0.983
#   h370_020_vhh            2.289 -> 2.101     no binder histidine
#   bc_s831683_mpnn6_S15D   5.397 -> 1.835     three, worst 0.661
#   bc_s831683_mpnn19_S15D  5.435 -> 1.774     three, worst 0.660
#   bc_s831683_mpnn9_S15D   5.428 -> 1.057     three, worst 0.344
#   bc_s831683_mpnn8_S15D   5.461 -> 1.023     three, worst 0.333
#   bc_d3acid_..._mpnn11    4.010 -> 0.737     three, worst 0.359
#   bc_s831683_mpnn9_WT     3.522 -> 0.593     three, worst 0.338
#
# Every binder histidine moves DOWN, 0.33 to 0.98, none up -- PROPKA noise would scatter both
# ways. Guard 4 of the combined gate fires on nearly all of them (nearest counter-charge 6.9-8.8
# A), so these are desolvation shifts with no electrostatic partner: the same mechanism as the
# 0.702x steric floor of METHODS 6 and the same physics that killed mechanism A.
#
# This is NOT a different objective. It is a less wrong estimate of the same one. The old number
# is retained as `ph_ratio_target_only_SUPERSEDED` so the change is auditable rather than silent.
#
# WHAT IT COSTS. Two VHH-format rows rise to ranks 1 and 5 on a pH estimate while their affinity
# is UNASSESSABLE by this instrument -- METHODS 4.5 shows it scores a measured 294 nM antibody
# below its own non-switching comparator. PK warned that "a nonbinding pose must not rise to the
# top through apparent selectivity." We are not demoting them for it, because demoting a design
# on an affinity reading we have shown to be inverted would be treating 0.219 as a measurement,
# which is the error METHODS 4.2 documents. Instead his other instruction is followed literally:
# "keep a diverse, eligible panel with separate columns for human binding evidence, mouse
# compatibility, pH hypothesis and uncertainty." The uncertainty is a column, not a demotion.


def multisite():
    """seq -> all-site pH product, pooled as the median over poses. The ranking key."""
    import json as _j
    return {k: v["allsite"] for k, v in _j.load(open(MULTISITE)).items()}


def pooled_affinity():
    """seq -> (hu_med, mo_med) pooled over every pose of that sequence, from master_rank.json."""
    out = {}
    for d in json.load(open(MASTER)):
        out[d["seq"]] = (d.get("hu_med"), d.get("mo_med"))
    return out


# `iptm_generator` DROPPED from the CSV, 2026-10-04. It was a leftover of the retired
# ipTM-descending order. Measured state before removal: blank for 8 of the 10 shipped rows,
# `0.23313` for one, and **`0.0` for `d2c_mpnn13_S88D_serasp`** — which to an outside reader
# says "the generator scored this design zero" when it means "no value was recorded". That is
# this project's signature failure (a 0.0 meaning ABSENT read as MEASURED) sitting in the one
# file Adaptyv grades, on a row carrying a real causal result. A column that is empty 80% of
# the time and misleading the rest is worse than no column.
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
    # rank on the all-site product when we have it; fall back to target-only and say so
    # A design with NO multi-site measurement cannot be ranked on the multi-site basis, and
    # must not be compared against one that can. Falling back to the target-only ratio looked
    # harmless and was not: it let bg04_r03 (1.937 target-only, no multi-site value) outrank
    # bc_s831683_mpnn6_S15D (5.397 target-only but 1.835 all-site), i.e. it compared two numbers
    # computed on different bases -- the error class this whole file exists to prevent. So a
    # missing multi-site value FAILS CLOSED out of tier 1 rather than borrowing the old number.
    key_ratio = r.get("allsite")
    if key_ratio is None:
        return (2, 0, -r["mo"], -r["hu"])      # unrankable on this basis: below both tiers
    tier1 = key_ratio >= RATIO_BAR and n >= MIN_N and binds

    # ASSESSABLE DESIGNS RANK AHEAD OF UNASSESSABLE ONES WITHIN TIER 1 (Harish, 2026-10-04).
    #
    # PK's ranking instruction was "apply eligibility and credible-interface checks FIRST, then
    # use the challenge priorities", with the guardrail "a nonbinding pose must not rise to the
    # top through apparent selectivity." On a pure pH ordering, rimA02_d3_rimA_14_vhh leads the
    # submission on a human ipSAE of 0.219 -- and we cannot say whether that is a weak interface
    # or an unreadable one, because METHODS 4.5 shows this instrument scores a measured 294 nM
    # antibody (G532, 0.0135) BELOW its own non-switching comparator (G532Ctrl, 0.2503) while
    # folding the Fv at 0.85. An antibody-format affinity reading here is uninterpretable.
    #
    # We do not demote them on the pH axis -- their ratios stand and are reported unchanged --
    # and we do not score them at 0.0000, which is the error METHODS 4.2 documents. We order
    # them after the designs whose affinity we CAN assess, so the row a reader reaches first is
    # one where both axes mean something.
    #
    # THE COST, stated rather than hidden: rimA02 carries the second-highest honest pH ratio in
    # the submission (4.838x all-site) and now sits below bc_s831683_mpnn19_S15D at 1.774x. If
    # Adaptyv rank strictly by the primary objective, this ordering costs us. It is a judgement
    # that credible-interface-first is the more defensible frame, not a claim that rimA02 is worse.
    unassessable = r.get("molecule_class") in ANTIBODY_CLASSES
    return (0 if tier1 else 1,
            (1 if unassessable else 0) if tier1 else 0,
            -key_ratio if tier1 else 0,
            -r["mo"], -r["hu"])


def main():
    rows = json.load(open(SUB))
    pooled = pooled_affinity()
    ms = multisite()
    missing = [x["name"] for x in rows if x["seq"] not in pooled]
    if missing:
        # Fail LOUD. The silent-zero path is this project's signature failure: `hu or 0.0`
        # turned an ABSENT measurement into a reported 0.0000, and a measured non-binder
        # scoring 0.0000 is indistinguishable from a design that was never scored.
        raise SystemExit(
            f"{len(missing)} submitted sequence(s) are absent from {MASTER}: "
            f"{', '.join(missing[:5])}. Re-run bin/master_rank.py. Refusing to emit a CSV "
            f"with affinity columns that would silently read 0.0000.")
    scored = []
    for x in rows:
        hu, mo = pooled[x["seq"]]
        scored.append(dict(name=x["name"], sequence=x["seq"],
                           molecule_class=x.get("molecule_class", MOLECULE_CLASS),
                           ratio=float(x["ratio"]), allsite=ms.get(x["seq"]),
                           hu=hu or 0.0, mo=mo or 0.0,
                           # ratio_n must be carried through: rank_key fails CLOSED without
                           # it, so dropping it here silently emptied tier 1 and put a
                           # non-switching binder at rank 1.
                           ratio_n=int(x.get("ratio_n", 0)),
                           # per-design assessment text must be carried through, or PK's
                           # uncertainty column silently collapses to one default string.
                           # Same class of bug as ratio_n, which emptied tier 1 earlier today.
                           assessment=x.get("assessment", "computational candidate"),
                           rmsd=x.get("rmsd", "")))

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
            "ph_poses_n", "ph_ratio_target_only_SUPERSEDED", "affinity_assessable",
            "affinity_above_null", "assessment"]
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in scored:
            # The headline pH column is the ALL-SITE product (both partners). The old
            # target-only number rides alongside as _SUPERSEDED so the re-rank is auditable.
            head = r.get("allsite")
            if head is None: head = r["ratio"]
            w.writerow([r["name"], r["sequence"], r["molecule_class"],
                        f"{head:.3f}", f"{r['hu']:.4f}", f"{r['mo']:.4f}",
                        r.get("ratio_n", 0),
                        f"{r['ratio']:.3f}",
                        # PK: "a failed run must not silently become a valid score of zero" --
                        # a VHH affinity reading here is not low, it is UNINTERPRETABLE. METHODS
                        # 4.5: this instrument scores a measured 294 nM antibody below its own
                        # non-switching comparator. So the column says so rather than implying
                        # the number means something.
                        "no -- antibody format, see METHODS 4.5"
                            if r["molecule_class"] in ("nanobody", "scfv", "fab_kappa", "fab_lambda")
                            else "yes",
                        # PK: separate columns for binding evidence, pH hypothesis and uncertainty
                        "yes" if max(r["hu"], r["mo"]) >= AFFINITY_FLAG else "no",
                        r.get("assessment", "computational candidate")])
    print(f"\nwrote {OUT}: {len(scored)} designs, {len(cols)} columns")
    print(f"{'#':>3} {'ratio':>7} {'mouse':>7} {'human':>7}  name")
    for i, r in enumerate(scored, 1):
        print(f"{i:>3} {r['ratio']:>7.3f} {r['mo']:>7.4f} {r['hu']:>7.4f}  {r['name'][:44]}")


def selftest():
    """Every fixture now carries `allsite`, because that is the ranking key as of 2026-10-04.

    A fixture WITHOUT it is not a tier-1 candidate at all -- see rank_key. The old fixtures
    carried only `ratio` and so silently fell to tier 2 the moment the basis changed, which is
    how this selftest caught the switch rather than the CSV catching it.
    """
    A = lambda **k: dict(ratio_n=5, **k)
    # tier ordering on the NEW basis
    a = A(ratio=4.5, allsite=4.5, mo=0.3, hu=0.3)
    b = A(ratio=0.0, allsite=0.0, mo=0.9, hu=0.9)
    assert rank_key(a) < rank_key(b), "a real switch must outrank a non-switching binder"
    c = A(ratio=0.0, allsite=0.0, mo=0.5, hu=0.9)
    d = A(ratio=0.0, allsite=0.0, mo=0.9, hu=0.1)
    assert rank_key(d) < rank_key(c), "mouse must outrank human affinity in tier 2"

    # THE BASIS ITSELF. A design whose target-only ratio is large but whose all-site product is
    # small must rank BELOW one whose all-site product is large. This is the whole re-rank, and
    # it is the live case: bc_s831683_mpnn8_S15D reads 5.461 target-only and 1.023 all-site,
    # while rimA02_d3_rimA_14_vhh reads 5.186 target-only and 4.838 all-site.
    s15d  = A(ratio=5.461, allsite=1.023, mo=0.764, hu=0.776)
    vhh   = A(ratio=5.186, allsite=4.838, mo=0.447, hu=0.219)
    # ASSESSABLE-FIRST, within tier 1 only. The live case: rimA02 (nanobody, 4.838 all-site)
    # must rank BELOW bc_s831683_mpnn19_S15D (protein, 1.774 all-site) despite a 2.7x higher
    # ratio, because the nanobody's affinity is unreadable on this instrument. That is the
    # deliberate cost recorded in rank_key.
    vhh_t1  = A(ratio=5.186, allsite=4.838, mo=0.447, hu=0.219, molecule_class="nanobody")
    prot_t1 = A(ratio=5.435, allsite=1.774, mo=0.786, hu=0.808, molecule_class="protein")
    assert rank_key(prot_t1) < rank_key(vhh_t1), \
        "an assessable tier-1 design outranks an unassessable tier-1 design"
    # but a tier-1 nanobody still outranks a tier-2 protein -- the format penalty does not
    # override the switch/no-switch split
    prot_t2 = A(ratio=5.461, allsite=1.023, mo=0.764, hu=0.776, molecule_class="protein")
    assert rank_key(vhh_t1) < rank_key(prot_t2), "tier 1 still beats tier 2 regardless of format"
    # and among assessable designs the all-site product drives the order
    hi = A(ratio=4.582, allsite=4.256, mo=0.567, hu=0.594, molecule_class="protein")
    assert rank_key(hi) < rank_key(prot_t1), "all-site product drives order among assessable designs"

    # MIXED BASIS MUST FAIL CLOSED. A candidate with no multi-site measurement cannot be
    # compared against one that has it. Before this guard, bg04_r03 (target-only 1.937, no
    # all-site value) outranked a design measured at 1.835 all-site -- two numbers on different
    # bases. Unrankable now means below both tiers.
    unmeasured = dict(ratio=9.9, mo=0.9, hu=0.9, ratio_n=9)          # no `allsite` key at all
    assert rank_key(s15d) < rank_key(unmeasured), "no multi-site value must fail closed"
    assert rank_key(A(ratio=0.1, allsite=0.1, mo=0.0, hu=0.0)) < rank_key(unmeasured), \
        "even a non-switching MEASURED design outranks an unmeasured one on this basis"

    BIND = (MIN_AFFINITY if MIN_AFFINITY is not None else AFFINITY_FLAG) + 0.01
    thin  = A(ratio=9.9, allsite=9.9, mo=BIND, hu=BIND); thin["ratio_n"] = 1
    solid = A(ratio=1.3, allsite=1.3, mo=BIND, hu=BIND)
    assert rank_key(solid) < rank_key(thin), "n<MIN_N must be excluded from tier 1"
    assert rank_key(solid) < rank_key(dict(ratio=9.9, allsite=9.9, mo=BIND, hu=BIND)), \
        "absent n must fail closed"

    nonbinder = A(ratio=9.9, allsite=9.9, mo=0.0, hu=0.0); nonbinder["ratio_n"] = 9
    if MIN_AFFINITY is None:
        assert rank_key(nonbinder) < rank_key(solid), "with the gate retired, ratio orders tier 1"
    else:
        assert rank_key(solid) < rank_key(nonbinder), "a 0.0000/0.0000 design must not reach tier 1"
    assert MOLECULE_CLASS in VALID_CLASSES
    print("selftest OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit(0)
    selftest(); main()
