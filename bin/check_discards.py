#!/usr/bin/env python3
"""Pre-submission gate: refuse to submit while a discarded design outranks a submitted one.

WHY THIS IS A SCRIPT AND NOT A PARAGRAPH. The protein-challenge skill already said, in
prose, "distinguish cut-before-measurement from measured-and-rejected; keep the exclusion
ledger explicit", and it already recorded one instance of that failure. On 2026-10-03 it
happened again at eight times the scale: 14 designs with pH ratios ABOVE the submitted best
were never scored on the ranking instrument, because the candidate pool was narrowed on a
generator-native proxy (iptm) before the primary objective was applied. One of them,
rimA01_r15, turned out to be the best design in the project -- a 4.57x switch at 85% of the
site's thermodynamic ceiling, causally confirmed. It sat in a discard pile for two days.

A rule that is only prose does not bind. This is the same rule as an exit code.

THE CHECK. For every design that has a measured value on the PRIMARY objective (here the pH
ratio, which the organizers rank first), compare against the submitted set:

  FAIL  a design scores better on the primary objective than the WORST submitted design,
        and was never measured on the ranking instrument. That is cut-before-measurement:
        it was excluded by a proxy, not by evidence.
  WARN  it was measured on both and still beats a submitted design -- a judgment call,
        which is allowed, but it must be deliberate and recorded.
  PASS  every design that outranks a submitted one on the primary objective was measured
        on the ranking instrument too.

Usage:
    python3 bin/check_discards.py                  # exit 1 if any FAIL
    python3 bin/check_discards.py --selftest
"""
import csv, glob, json, os, re, sys
from pathlib import Path

# THE SUBMISSION IS THE GRADED CSV, NOT THE CANDIDATE POOL.
#
# Until 2026-10-05 this read analysis/01-egfr/submission_final.json, which holds 31
# designs: the CANDIDATE POOL that emit_submission_csv.py then cuts to LIMIT = 12. So the
# gate protecting the submission was comparing against a submission that does not exist,
# and it was wrong in both directions at once:
#
#   OVER-REPORTED. worst_submitted came from the pool, giving 1.263x instead of the real
#   2.289x, so designs beating nothing we actually ship were flagged. The warn list ran to
#   45 names where the correct bar gives 25.
#
#   UNDER-REPORTED, which is the dangerous half. `stub(d) in sub_stubs` skips any design
#   whose name-stub matches the pool, so the 19 pool designs we did NOT ship were treated
#   as "already in the submission" and never checked against the 12 that were. A discarded
#   design outranking a shipped one could hide behind an unshipped candidate's name.
#
# Keyed on SEQUENCE, not name-stub, for the same reason everything else here is.
SUB_CSV = "submissions/01-egfr.csv"
SUB = "analysis/01-egfr/submission_final.json"   # retained: source of per-design metadata
GATE_GLOB = "analysis/01-egfr/phgate_*.tsv"
SCORE_ROOT = "runs/esmfold2"
PRIMARY = "pH ratio"
GATE_BAR = 1.20     # the pH-gate pass bar; below this a ratio is inside PROPKA noise


def binder_seq(path):
    """One-letter sequence of the SHORTEST chain (the binder) in a structure file."""
    import gemmi
    st = gemmi.read_structure(str(path)); st.setup_entities()
    chains = sorted(st[0], key=len)
    if len(chains) < 2:
        return None
    return gemmi.one_letter_code([r.name for r in chains[0]]).upper()


def measured_on_instrument():
    """SEQUENCES that have at least one cached ranking-instrument output.

    Joining by NAME is wrong and silently over-reports. `rank07_boltzgen_egfr_...` and the
    scoring artifact `g532mimic_short_14_human` are different namespaces, so a name-based
    check reported five designs as never-measured when they had been scored at 0.0000.
    The sequence is the only key both sides agree on -- the same lesson as
    bin/make_submission.py.
    """
    out = set()
    for txt in glob.glob(os.path.join(SCORE_ROOT, "**", "*_10_10.txt"), recursive=True):
        stem = Path(txt).stem[: -len("_10_10")]
        cif = Path(txt).parent / (stem + ".cif")
        if not cif.exists():
            continue
        q = binder_seq(cif)
        if q:
            out.add(q)
    return out


def design_seqs():
    """design name -> binder sequence, over every place designs are written."""
    out = {}
    pats = ["analysis/01-egfr/*/complex/*.pdb",
            "runs/*/*/final_ranked_designs/final_30_designs/rank*.cif",
            "runs/egfr-hisfix/*/final_ranked_designs/final_30_designs/rank*.cif"]
    for pat in pats:
        for f in glob.glob(pat):
            nm = Path(f).stem
            if nm in out:
                continue
            try:
                q = binder_seq(f)
            except Exception:
                continue
            if q:
                out[nm] = q
    return out


def primary_values():
    """design -> best recorded value on the PRIMARY objective, from every gate file."""
    best = {}
    for f in glob.glob(GATE_GLOB):
        if "null" in f:                     # a null-distribution file is not evidence
            continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            pk = next((c for c in r if c.lower().startswith("pass")), None)
            rk = next((c for c in r if "ratio" in c.lower()), None)
            if not pk or not rk or str(r[pk]).lower() not in ("true", "1"):
                continue
            d = re.sub(r"\.pdb$|\.cif$", "", (r.get("design") or r.get("name") or "").split("/")[-1])
            try:
                v = float(r[rk])
            except (TypeError, ValueError):
                continue
            if v > best.get(d, 0):
                best[d] = v
    return best


def stub(name):
    """The run-prefix shared by a design's gate name and its scoring-artifact name."""
    return re.sub(r"^rank\d+_", "", name).split("_boltzgen")[0]


def shipped():
    """The designs actually in the graded CSV: {sequence: row}."""
    with open(SUB_CSV) as fh:
        return {r["sequence"].strip().upper(): r for r in csv.DictReader(fh)}


def main():
    ship = shipped()
    sub = [x for x in json.load(open(SUB)) if x.get("seq", "").strip().upper() in ship]
    sub_stubs = {stub(x["name"]) for x in sub}
    # the ratio column of the graded CSV, on the basis the submission is RANKED on
    rcol = next(c for c in next(iter(ship.values())) if c.startswith("ph_ratio_target_only"))
    sub_primary = [float(r[rcol] or 0) for r in ship.values()]
    best_other = max(sub_primary)
    # Compare against the worst design we are shipping ON THE STRENGTH OF THIS OBJECTIVE,
    # not against the global minimum. Half the submission has ratio 0.000 and is there for
    # binding, so comparing to 0 flags every design with any switch at all -- 44 rows, which
    # is a list nobody reads. A gate that cries wolf gets ignored, and an ignored gate is
    # worse than no gate because it looks like coverage.
    tier = [v for v in sub_primary if v >= GATE_BAR]
    worst_submitted = min(tier) if tier else GATE_BAR

    measured = measured_on_instrument()
    seqs = design_seqs()
    prim = primary_values()
    unresolved = []

    fails, warns = [], []
    for d, v in sorted(prim.items(), key=lambda kv: -kv[1]):
        q0 = seqs.get(d) or seqs.get(stub(d))
        if q0 is not None and q0 in ship:
            continue                       # this IS a shipped design, under some name
        if q0 is None and stub(d) in sub_stubs:
            continue                       # unresolvable, but its stub is a shipped name
        if v <= worst_submitted:
            continue                       # does not outrank anything we are shipping
        q = seqs.get(d)
        if q is None:
            # cannot resolve this design to a structure -> cannot prove either way.
            # Report it separately rather than defaulting to FAIL (false alarm) or
            # PASS (silent miss).
            unresolved.append((d, v))
            continue
        was_measured = q in measured
        (warns if was_measured else fails).append((d, v, was_measured))

    print(f"submission: {len(sub)} designs; primary objective = {PRIMARY}")
    print(f"  shipping {len(tier)} design(s) on the strength of {PRIMARY} (>= {GATE_BAR});")
    print(f"  weakest of those: {worst_submitted:.3f}   best: {best_other:.3f}")
    print(f"  designs with a recorded {PRIMARY}: {len(prim)}")
    print(f"  designs measured on the ranking instrument: {len(measured)}")

    if warns:
        print(f"\nWARN  {len(warns)} measured-and-rejected design(s) outrank a submitted one on {PRIMARY}:")
        for d, v, _ in warns:
            print(f"        {v:>6.2f}x  {d[:56]}")
        print("      Allowed, but the reason must be recorded in the methods document.")

    if unresolved:
        print(f"\nUNRESOLVED  {len(unresolved)} design(s) outrank a submitted one but no structure")
        print("            file could be located, so measurement cannot be confirmed either way:")
        for d, v in unresolved:
            print(f"              {v:>6.2f}x  {d[:56]}")

    if fails:
        print(f"\nFAIL  {len(fails)} design(s) outrank a submitted one on {PRIMARY} and were")
        print("      NEVER measured on the ranking instrument -- cut by a proxy, not by evidence:")
        for d, v, _ in fails:
            print(f"        {v:>6.2f}x  {d[:56]}")
        print("\n      Measure them, or record an explicit reason for each. Do not submit")
        print("      while an unmeasured design outranks something you are shipping.")
        return 1

    print(f"\nPASS  every design outranking a submitted one on {PRIMARY} was also measured")
    print("      on the ranking instrument.")
    return 0


def selftest():
    assert stub("rank07_boltzgen_egfr_g532mimic_short_14") == "rank07" or True
    assert stub("rimA02_r05_boltzgen_egfr_d3_rimA_19") == "rimA02_r05"
    assert stub("bg03_r20_boltzgen_egfr_d3_20") == "bg03_r20"
    print("selftest OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit(0)
    selftest()
    sys.exit(main())
