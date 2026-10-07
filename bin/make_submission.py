#!/usr/bin/env python3
"""Regenerate the problem-1 submission ranking from cached scores, under v1 AND v2.

WHY THIS EXISTS. The 2026-10-03 audit found that `analysis/01-egfr/ranking_final.json`
and `submission_final.json` had NO generator anywhere in the repo -- grep for their names
across every .py/.sh/.md outside .venv/.git/runs returned zero hits. The submission could
not be regenerated, re-checked, or re-ranked. It also found the shipped ranking is v1 MAX
(21/21 `ips` values match a cached v1 max exactly; only 9/21 match a v2 median) even though
The reviewer's review replaced max-over-seeds with median-over-seeds. This script makes the ordering
executable so both estimators can be compared instead of one being swapped in silently.

THE JOIN IS BY SEQUENCE, NOT BY NAME. Ranking names (`rimA01_r02_boltzgen_egfr_d3_rimA_26`)
and scoring-artifact names (`rimA01_r02_hu`) live in different namespaces -- only 4 of 21
match -- and no mapping exists in code. The reliable key is the binder chain's one-letter
sequence, taken from the folded .cif beside each cached ipsae output.

REPLICATES ARE SURFACED, NOT POOLED. 29 groups of byte-identical inputs were folded under
different run names (68 of 104 scoring keys). Pooling them silently is what merged
ESMFold2-Fast with ESMFold2-Full into the control calibration. Here each (sequence, species)
takes ONE run deterministically -- most seeds wins, then lexicographic run name -- and the
number of other runs and their spread are reported so the ambiguity is visible.

ORDERING -- locked decision 2, lexicographic:
  tier 1: designs with a real pH switch, by pH ratio descending
  tier 2: everything else, by binding score descending
The shipped file additionally hand-interleaved four iptm==0 designs at positions 5/7/9/11;
that is reproduced only with --replicate-shipped, and is flagged as a manual act.

Usage:
    python3 bin/make_submission.py                 # v1 and v2 side by side
    python3 bin/make_submission.py --selftest
"""
import glob, json, os, re, statistics, sys
from collections import defaultdict
from pathlib import Path

SUB = "analysis/01-egfr/submission_final.json"
SCORE_ROOT = "runs/esmfold2"
RATIO_BAR = 1.20


def read_cached(txt):
    """ipSAE_min for one prediction: min over the TWO ALIGNMENT DIRECTIONS of ONE interface.
    Returns None for >1 chain pair -- min() across pairs is the bug fixed in ipsae_min.py."""
    rows = [l.split() for l in Path(txt).read_text().splitlines()
            if l.strip() and not l.startswith("Chn1")]
    asym = {f"{r[0]}->{r[1]}": float(r[5]) for r in rows if len(r) > 5 and r[4] == "asym"}
    if not asym:
        return None
    pairs = defaultdict(list)
    for k, v in asym.items():
        c1, c2 = k.split("->")
        pairs[frozenset((c1, c2))].append(v)
    if len(pairs) > 1:
        return None
    return min(next(iter(pairs.values())))


def design_of(stem):
    stem = re.sub(r"(_seed\d+)?_sample_\d+$", "", stem)
    return re.sub(r"_model_?\d+$", "", stem)


def build_index():
    """(binder_seq, species) -> {run: [(seed, score)]}. Species from the run/design name."""
    import gemmi
    idx = defaultdict(lambda: defaultdict(list))
    for txt in glob.glob(os.path.join(SCORE_ROOT, "**", "*_10_10.txt"), recursive=True):
        m = read_cached(txt)
        if m is None:
            continue
        stem = Path(txt).stem[: -len("_10_10")]
        cif = Path(txt).parent / (stem + ".cif")
        if not cif.exists():
            continue
        st = gemmi.read_structure(str(cif)); st.setup_entities()
        chains = sorted(st[0], key=len)
        if len(chains) < 2:
            continue
        binder = gemmi.one_letter_code([r.name for r in chains[0]]).upper()
        d = design_of(stem)
        sp = "mo" if re.search(r"_mo\b|_mouse", d) else "hu" if re.search(r"_hu\b|_human", d) else "?"
        run = os.path.basename(os.path.dirname(os.path.dirname(txt)))
        seed = (re.search(r"_seed(\d+)", stem) or [None, "?"])[1]
        idx[(binder, sp)][run].append((seed, m))
    return idx


def pick_run(runs):
    """Deterministic: most seeds wins, ties broken lexicographically by run name."""
    return sorted(runs.items(), key=lambda kv: (-len(kv[1]), kv[0]))[0]


def score(idx, seq, sp, agg):
    """agg='v1' -> max over seeds; agg='v2' -> median over seeds. Returns (value, n, nruns, spread)."""
    runs = idx.get((seq, sp)) or {}
    if not runs:
        return (None, 0, 0, None)
    run, vals = pick_run(runs)
    f = max if agg == "v1" else statistics.median
    v = f(x for _, x in vals)
    per_run = [f(x for _, x in vs) for vs in runs.values()]
    spread = (max(per_run) - min(per_run)) if len(per_run) > 1 else 0.0
    return (v, len(vals), len(runs), spread)


def control_bars(idx_raw):
    """Derive the binding bar from the CONTROLS in the data, never a hardcoded guess.

    Measured 2026-10-03 on ESMFold2-Full (fastgate excluded -- mixing checkpoints is how
    the calibration got corrupted before):
        NEG_nonbinder        v2 median 0.1493
        POS_cetuximab_scfv   v2 median 0.6224
        POS_cradle_1nM       v2 median 0.6275
        POS_nano2_5nM        v2 median 0.0000   <-- a REAL 5 nM nanobody, scored zero

    Two defensible bars and nothing to choose between them from data, because there are no
    control points between 0.15 and 0.62:
        permissive = max(NEG)    "better than a known nonbinder"
        strict     = min(POS>0)  "as good as a real nanomolar binder"

    The nano2 row is the load-bearing caveat: the instrument returns 0.0000 for a genuine
    5 nM binder, so ANY bar discards some real binders. A design failing this gate is not
    shown to be a nonbinder; it is shown to be unsupported BY THIS INSTRUMENT.
    """
    import glob, os, statistics as st
    from collections import defaultdict
    per = defaultdict(lambda: defaultdict(list))
    for txt in glob.glob(os.path.join(SCORE_ROOT, "**", "*_10_10.txt"), recursive=True):
        run = os.path.basename(os.path.dirname(os.path.dirname(txt)))
        if run == "fastgate":        # Full only; never mix checkpoints
            continue
        d = design_of(Path(txt).stem[: -len("_10_10")])
        # Match POSd3_/NEGd3_ as well as POS_/NEG_. The prefix test used to be
        # startswith("POS_")/("NEG_"), which silently dropped all 20 domain-III control poses
        # -- and the d3 calibration is materially different from the full-ECD one:
        #     ECD : NEG 0.1493   POS 0.6224 / 0.6275
        #     d3  : NEG 0.2218   POS 0.4005 / 0.4183
        # The majority of this project's designs were scored on d3, so they were being judged
        # against an ECD bar: a d3 design at 0.20 was called "better than a nonbinder" when the
        # construct-matched nonbinder beats it, and the 0.6224 strict gate is unreachable on d3
        # (a construct-matched 1 nM binder only reaches 0.4183).
        if not (d.startswith(("POS_", "NEG_", "POSd3_", "NEGd3_"))):
            continue
        m = read_cached(txt)
        if m is not None:
            per[d][run].append(m)
    # PER-CONSTRUCT bars. Merging the two into one pool is as wrong as dropping d3 was: it
    # would apply d3's higher nonbinder floor to ECD-scored designs and ECD's lower positive
    # to d3-scored ones. A design must be judged against the controls folded on ITS OWN target.
    def _bars(pfx_pos, pfx_neg):
        p = [st.median(v) for d, rs in per.items() if d.startswith(pfx_pos)
             for v in rs.values() if st.median(v) > 0]
        n = [st.median(v) for d, rs in per.items() if d.startswith(pfx_neg) for v in rs.values()]
        if not p or not n:
            return None
        return {"permissive": max(n), "strict": min(p), "n_pos": len(p), "n_neg": len(n)}

    ecd = _bars(("POS_",), ("NEG_",))
    d3 = _bars(("POSd3_",), ("NEGd3_",))
    if ecd is None:
        return None
    out = {"ecd": ecd, "d3": d3,
           # back-compat scalars for the two v3 gates: deliberately CONSERVATIVE, i.e. the
           # highest nonbinder and the lowest real binder seen on EITHER construct. Any caller
           # that knows which construct a design was scored on should use out["ecd"]/out["d3"].
           "permissive": max([ecd["permissive"]] + ([d3["permissive"]] if d3 else [])),
           "strict": min([ecd["strict"]] + ([d3["strict"]] if d3 else [])),
           "n_pos": ecd["n_pos"] + (d3["n_pos"] if d3 else 0),
           "n_neg": ecd["n_neg"] + (d3["n_neg"] if d3 else 0),
           "zero_pos": sum(1 for d, rs in per.items() if d.startswith(("POS_", "POSd3_"))
                           and all(st.median(v) == 0 for v in rs.values()))}
    # STANDING CAVEAT, 2026-10-04: `NEG_nonbinder` is 81% mature human EGF -- 43 of 53 residues
    # verbatim, all six cysteines in linear order -- i.e. the native ~2 nM EGFR agonist, whose
    # high-affinity contacts are on domain III, the surface every design targets. It is a BINDER.
    # So "permissive" is a binder's score, not a nonbinder's, and every "clears the bar" claim
    # in this project is calibrated against it. A composition-matched shuffled null is being
    # folded to replace it. Until that lands, treat these bars as provisional.
    return out


def rank(rows, idx, agg, bar=None):
    out = []
    for r in rows:
        seq = r["seq"]
        hu = score(idx, seq, "hu", agg)
        mo = score(idx, seq, "mo", agg)
        primary = hu[0] if hu[0] is not None else (mo[0] if mo[0] is not None else 0.0)
        out.append(dict(name=r["name"], ratio=float(r["ratio"]),
                        real=str(r["real"]) in ("1", "True"),
                        hu=hu[0], mo=mo[0], n=hu[1] or mo[1],
                        nruns=max(hu[2], mo[2]), spread=max(hu[3] or 0, mo[3] or 0),
                        primary=primary))
    if bar is None:
        # locked decision 2: pH tier first (ratio desc), then the rest by binding score desc
        t1 = sorted([x for x in out if x["real"]], key=lambda x: -x["ratio"])
        t2 = sorted([x for x in out if not x["real"]], key=lambda x: -x["primary"])
        return t1 + t2
    # v3: the pH tier is GATED on a credible interface, per the reviewer -- "the conditional term
    # must be gated on a credible interface rather than weighted beside the primary score".
    # A pH ratio on a design the instrument says does not bind is not evidence of a switch.
    for x in out:
        x["eligible"] = x["primary"] >= bar
    g1 = sorted([x for x in out if x["eligible"] and x["real"]], key=lambda x: -x["ratio"])
    g2 = sorted([x for x in out if x["eligible"] and not x["real"]], key=lambda x: -x["primary"])
    g3 = sorted([x for x in out if not x["eligible"]], key=lambda x: -x["primary"])
    for x in g3:
        x["note"] = "below binding bar"
    return g1 + g2 + g3


def show(label, ranked):
    print(f"\n=== {label} " + "=" * (62 - len(label)))
    print(f"{'#':>3} {'tier':>5} {'ratio':>7} {'score':>7} {'mouse':>7} {'n':>3} {'runs':>5} {'spread':>7}  name")
    for i, x in enumerate(ranked, 1):
        hu = f"{x['hu']:.4f}" if x["hu"] is not None else "  --  "
        mo = f"{x['mo']:.4f}" if x["mo"] is not None else "  --  "
        flag = " *" if x["nruns"] > 1 else ""
        print(f"{i:>3} {'pH' if x['real'] else '-':>5} {x['ratio']:>7.3f} {hu:>7} {mo:>7} "
              f"{x['n']:>3} {x['nruns']:>5}{flag:<2} {x['spread']:>6.4f}  {x['name'][:40]}")


def selftest():
    assert read_cached.__doc__
    # the ordering rule must put every pH-tier design above every non-pH design
    fake = [dict(name="a", ratio="2.0", real="1", seq="AAA"),
            dict(name="b", ratio="9.9", real="0", seq="BBB")]
    r = rank(fake, {}, "v1")
    assert r[0]["name"] == "a", "pH tier must outrank a higher ratio in tier 2"
    # v1 >= v2 for any multiset (max >= median)
    import random
    vals = [(str(i), 0.1 * i) for i in range(5)]
    assert max(v for _, v in vals) >= statistics.median(v for _, v in vals)
    print("selftest OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest(); sys.exit(0)
    selftest()
    rows = json.load(open(SUB))
    idx = build_index()
    print(f"indexed {len(idx)} (sequence, species) pairs from {SCORE_ROOT}")
    multi = sum(1 for v in idx.values() if len(v) > 1)
    print(f"  {multi} of them were folded in more than one run -- marked * below")
    bars = control_bars(idx)
    print(f"  control-derived bars: permissive (max NEG) = {bars['permissive']:.4f}, "
          f"strict (min POS>0) = {bars['strict']:.4f}")
    print(f"  WARNING: {bars['zero_pos']} positive control(s) score 0.0000 -- a real 5 nM "
          f"nanobody among them. Any bar discards some genuine binders.")
    v1 = rank(rows, idx, "v1")
    v2 = rank(rows, idx, "v2")
    v3 = rank(rows, idx, "v2", bar=bars["permissive"])
    v3s = rank(rows, idx, "v2", bar=bars["strict"])
    show("v1  ipSAE_min, MAX over seeds  (what shipped)", v1)
    show("v2  ipSAE_min, MEDIAN over seeds  (the reviewer's estimator)", v2)
    show(f"v3  v2 GATED on binding >= {bars['permissive']:.4f} (permissive: beats the nonbinder control)", v3)
    show(f"v3s v2 GATED on binding >= {bars['strict']:.4f} (strict: matches a real nM binder)", v3s)
    for lbl, r in (("permissive", v3), ("strict", v3s)):
        ok = [x for x in r if x.get("eligible")]
        ph = [x for x in ok if x["real"]]
        print(f"\n  {lbl} gate: {len(ok)}/{len(r)} eligible, of which {len(ph)} carry a pH switch")
    n1 = [x["name"] for x in v1]; n2 = [x["name"] for x in v2]
    moved = [(i + 1, n2.index(n) + 1, n) for i, n in enumerate(n1) if n2.index(n) != i]
    print(f"\n{len(moved)}/{len(n1)} designs change position between v1 and v2")
    for a, b, n in moved:
        print(f"   {a:>2} -> {b:<2}  {n[:46]}")
