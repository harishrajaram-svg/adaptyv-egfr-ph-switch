#!/usr/bin/env python3
"""Calibrate this project's ranking instruments against MEASURED TNF-alpha labels.

Anthropic's released campaign (CC BY 4.0) ordered 150 de novo TNF-alpha designs and
measured them: 12 bound, 138 did not, same assay vendor as our challenge
(Acro TNA-H4211). Every design carries co-folding metrics from eight models. This is
the positive control arms-backlog.md 3a called a blocking dependency and never got.

What this answers, which nothing in this repo could answer before:
  1. On the real target, does any metric separate binders from non-binders at all?
  2. Does OUR 0.45 iptm bar admit binders and reject non-binders, or neither?
  3. Does the monomer-foldability floor (protocol default 0.70) behave as a filter?
  4. Which co-folding model is the best discriminator -- i.e. which oracle to trust?
  5. Where do OUR 35 designs fall against 150 designs that were actually assayed?

HONESTY ABOUT PRE-REGISTRATION. This is a calibration of existing instruments against
external labels, not a test of a design hypothesis, and it is NOT blind: median values
for several models were inspected before the AUC statistic was chosen. The guard against
that is completeness, not blinding -- every model x every metric is printed, so no
subset can be quietly selected. Treat single-cell AUCs as descriptive.

Statistic: design-level AUC (Mann-Whitney, tie-corrected). Seeds are replicates of one
design, so they are collapsed by median per (design, model, stoichiometry) BEFORE
ranking; ranking raw seed rows would inflate n 10x and understate the error.
Significance by permutation on the labels (n=20000), because n_pos=12.

Direction: for PAE-like metrics lower is better, so they are negated before ranking.
AUC 0.5 = chance. AUC < 0.5 = anti-correlated.

  calibrate_external.py              the calibration
  calibrate_external.py --selftest   self-tests, incl. two mutation tests
"""
import csv
import glob
import os
import random
import statistics
import sys

print = __import__('functools').partial(print, flush=True)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
REF = os.path.join(REPO, "reference", "anthropic-campaign")
RUNS = os.path.join(REPO, "runs", "mosaic-p2")

# metric -> higher_is_better
METRICS = {
    "ipsae_min": True, "ipsae_max": True, "iptm_pae": True, "plddt_binder": True,
    "sc_dockq": True, "n_interface_contacts": True,
    "pae_interface_min": False, "pae_interface_mean": False,
}
OUR_BAR_IPTM = 0.45          # D-P2-1, this project's geometry gate
PROTOCOL_PLDDT_FLOOR = 70.0  # anthropic-binder-design-protocol.md:100, 0-100 scale


def num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f  # NaN -> None, never silently 0.0


def auc(pos, neg):
    """Tie-corrected Mann-Whitney AUC. None when either class is empty."""
    if not pos or not neg:
        return None
    vals = sorted(pos + neg)
    ranks, i = {}, 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[i]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        ranks[vals[i]] = r
        i = j + 1
    s = sum(ranks[v] for v in pos)
    return (s - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def perm_p(pos, neg, n=20000, seed=0):
    """One-sided permutation p on the labels: P(AUC_shuffled >= AUC_observed)."""
    obs = auc(pos, neg)
    if obs is None:
        return None
    allv, k = pos + neg, len(pos)
    rng = random.Random(seed)
    hits = 0
    for _ in range(n):
        rng.shuffle(allv)
        if auc(allv[:k], allv[k:]) >= obs:
            hits += 1
    return (hits + 1) / (n + 1)


def load_external():
    """(labels, per_design) -- per_design[(model, stoich, metric)][uuid] = median over seeds."""
    dpath = os.path.join(REF, "tnfa_designs.tsv")
    cpath = os.path.join(REF, "tnfa_cofold.tsv")
    for p in (dpath, cpath):
        if not os.path.exists(p):
            sys.exit(f"REFUSE: {p} missing. Run bin/fetch-anthropic-campaign.sh first.")

    labels, meta = {}, {}
    for r in csv.DictReader(open(dpath), delimiter="\t"):
        v = r["binder_final"].strip().lower()
        if v not in ("true", "false"):
            sys.exit(f"REFUSE: unparseable binder_final {r['binder_final']!r}")
        labels[r["uuid"]] = (v == "true")
        meta[r["uuid"]] = r
    if len(labels) != 150 or sum(labels.values()) != 12:
        sys.exit(f"REFUSE: expected 150 designs / 12 binders, got "
                 f"{len(labels)} / {sum(labels.values())}")

    acc = {}
    for r in csv.DictReader(open(cpath), delimiter="\t"):
        if r["uuid"] not in labels:
            continue
        for m in METRICS:
            v = num(r.get(m))
            if v is None:
                continue
            acc.setdefault((r["cofolding_model"], str(r["stoichiometry"]), m), {}) \
               .setdefault(r["uuid"], []).append(v)
    per = {k: {u: statistics.median(vs) for u, vs in d.items()} for k, d in acc.items()}
    return labels, meta, per


def load_ours():
    """Our 35 Mosaic designs: iptm_repred + plddt_binder_repred (Boltz-2, 1 binder:trimer)."""
    out = []
    for f in sorted(glob.glob(os.path.join(RUNS, "*", "designs.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            out.append({
                "design": r.get("design", "?"),
                "cond": os.path.basename(os.path.dirname(f)),
                "iptm": num(r.get("iptm_repred")),
                "plddt": num(r.get("plddt_binder_repred")),
            })
    return out


def report():
    labels, meta, per = load_external()
    pos_ids = {u for u, b in labels.items() if b}
    print(f"External labels: {len(labels)} TNF-alpha designs, "
          f"{len(pos_ids)} measured binders ({100*len(pos_ids)/len(labels):.1f}%)")
    print("Source: Anthropic/claude-protein-binder-design, CC BY 4.0. "
          "Assay vendor Acro TNA-H4211 -- the same construct as our challenge.\n")

    # --- 1. generator and length, which is what actually moved the outcome ---
    print("=" * 78)
    print("1. WHAT PRODUCED THE 12 BINDERS")
    print("=" * 78)
    gen = {}
    for u, r in meta.items():
        g = gen.setdefault(r["generator"], [0, 0])
        g[1] += 1
        g[0] += labels[u]
    print(f"  {'generator':30}{'hits':>6}{'n':>6}{'rate':>8}")
    for g, (h, n) in sorted(gen.items(), key=lambda x: -x[1][0] / max(x[1][1], 1)):
        print(f"  {g[:29]:30}{h:>6}{n:>6}{100*h/n:>7.1f}%")
    bl = sorted(int(meta[u]["binder_length"]) for u in pos_ids)
    al = [int(r["binder_length"]) for r in meta.values() if r["binder_length"]]
    print(f"\n  binder lengths : {bl}")
    print(f"  all 150 tested : min {min(al)}  max {max(al)}")
    sd = {meta[u]["sequence_design_method"] for u in pos_ids}
    print(f"  sequence design on all 12 binders: {sorted(sd)}")
    mouse = [u for u in pos_ids if meta[u]["mouse_binding_final"].strip().lower() == "binder"]
    print(f"\n  ALSO BIND MOUSE (our objective 2): {len(mouse)} of {len(pos_ids)}")
    for u in sorted(mouse, key=lambda x: float(meta[x]["kd_nM_final"] or 9e9)):
        r = meta[u]
        print(f"    {r['generator'][:12]:13} L={r['binder_length']:>3}  "
              f"human KD {float(r['kd_nM_final']):8.2f} nM   mouse KD {r['mouse_kd_nM_final']:>8} nM")

    # --- 2. every model x every metric, no selection ---
    print("\n" + "=" * 78)
    print("2. DOES ANY METRIC SEPARATE BINDERS FROM NON-BINDERS?  (design-level AUC)")
    print("=" * 78)
    print("   AUC 0.50 = chance.  p by permutation on labels, n=20000, one-sided.")
    best = []
    for stoich in sorted({k[1] for k in per}):
        print(f"\n--- stoichiometry {stoich}"
              f"{'   <-- our construct: 1 binder : 3 protomers' if stoich == '1to3' else ''}")
        print(f"  {'model':10}{'metric':20}{'AUC':>7}{'p':>9}"
              f"{'binder med':>12}{'non-bind med':>14}")
        for model in sorted({k[0] for k in per}):
            for m, higher in METRICS.items():
                d = per.get((model, stoich, m))
                if not d:
                    continue
                sign = 1.0 if higher else -1.0
                p_ = [sign * v for u, v in d.items() if u in pos_ids]
                n_ = [sign * v for u, v in d.items() if u not in pos_ids]
                a = auc(p_, n_)
                if a is None:
                    continue
                best.append((a, model, stoich, m, list(p_), list(n_)))
                print(f"  {model:10}{m:20}{a:>7.3f}{'':>9}"
                      f"{sign*statistics.median(p_):>12.3f}{sign*statistics.median(n_):>14.3f}")

    print("\n  TOP 8 DISCRIMINATORS across all models/stoichiometries/metrics")
    print("  (permutation p computed only here -- 20000 label shuffles each):")
    ranked = sorted(best, key=lambda x: -x[0])[:8]
    for a, model, st, m, p_, n_ in ranked:
        pv = perm_p(p_, n_, n=20000)
        print(f"    AUC {a:.3f}  p={pv:.5f}   {model} / {st} / {m}")
    print(f"  Bonferroni bar for the {len(best)} cells scanned, alpha 0.05: "
          f"p < {0.05/len(best):.5f}")
    print("  WORST 3 (anti-correlated = the metric points the wrong way):")
    for a, model, st, m, p_, n_ in sorted(best, key=lambda x: x[0])[:3]:
        print(f"    AUC {a:.3f}            {model} / {st} / {m}")

    # --- 3. our own bars, on real labels ---
    print("\n" + "=" * 78)
    print("3. OUR BARS, SCORED AGAINST REAL OUTCOMES")
    print("=" * 78)
    for model in sorted({k[0] for k in per}):
        d = per.get((model, "1to3", "iptm_pae"))
        if not d:
            continue
        passed = [u for u, v in d.items() if v >= OUR_BAR_IPTM]
        tp = len([u for u in passed if u in pos_ids])
        print(f"  iptm >= {OUR_BAR_IPTM} on {model:9}: admits {len(passed):3}/{len(d)} designs, "
              f"{tp:2}/{len(pos_ids)} binders kept, precision "
              f"{(100*tp/len(passed)) if passed else float('nan'):5.1f}%")
    print()
    for model in sorted({k[0] for k in per}):
        d = per.get((model, "1to3", "plddt_binder"))
        if not d:
            continue
        passed = [u for u, v in d.items() if v >= PROTOCOL_PLDDT_FLOOR]
        tp = len([u for u in passed if u in pos_ids])
        lo = min(d[u] for u in pos_ids if u in d)
        print(f"  plddt >= {PROTOCOL_PLDDT_FLOOR:.0f} on {model:9}: admits {len(passed):3}/{len(d)}, "
              f"{tp:2}/{len(pos_ids)} binders kept, lowest binder {lo:5.1f}")

    # --- 4. where we sit ---
    print("\n" + "=" * 78)
    print("4. OUR 35 DESIGNS vs 150 THAT WERE ACTUALLY ASSAYED")
    print("=" * 78)
    ours = load_ours()
    oi = [r["iptm"] for r in ours if r["iptm"] is not None]
    op = [r["plddt"] for r in ours if r["plddt"] is not None]
    ref_i = per.get(("boltz2", "1to3", "iptm_pae"), {})
    ref_p = per.get(("boltz2", "1to3", "plddt_binder"), {})
    print("  Comparison model: boltz2, 1 binder : 3 protomers -- our construct and our oracle.")
    print("  CAVEAT: their co-folding config is not byte-identical to our re-prediction "
          "(seeds, MSA,\n  templates), so read this as an order-of-magnitude placement, "
          "not a matched measurement.\n")
    if ref_i:
        rv = sorted(ref_i.values())
        print(f"  iptm         theirs n={len(rv):3}  min {rv[0]:.3f}  median "
              f"{statistics.median(rv):.3f}  max {rv[-1]:.3f}")
        print(f"               OURS   n={len(oi):3}  min {min(oi):.3f}  median "
              f"{statistics.median(oi):.3f}  max {max(oi):.3f}")
        print(f"               our best ({max(oi):.3f}) vs their WORST ({rv[0]:.3f}): "
              f"{'ABOVE' if max(oi) > rv[0] else 'BELOW the floor of the real field'}")
    if ref_p:
        rv = sorted(ref_p.values())
        # their plddt is 0-100, ours 0-1
        print(f"\n  plddt_binder theirs n={len(rv):3}  min {rv[0]:.1f}   median "
              f"{statistics.median(rv):.1f}   max {rv[-1]:.1f}   (0-100)")
        print(f"               OURS   n={len(op):3}  min {100*min(op):.1f}   median "
              f"{100*statistics.median(op):.1f}   max {100*max(op):.1f}   (rescaled)")
        print(f"               clearing the {PROTOCOL_PLDDT_FLOOR:.0f} protocol floor: "
              f"theirs {sum(1 for v in rv if v >= PROTOCOL_PLDDT_FLOOR)}/{len(rv)}, "
              f"OURS {sum(1 for v in op if 100*v >= PROTOCOL_PLDDT_FLOOR)}/{len(op)}")
    return 0


def selftest():
    # AUC: perfect, inverted, chance, ties
    assert auc([3, 4], [1, 2]) == 1.0
    assert auc([1, 2], [3, 4]) == 0.0
    assert abs(auc([1, 2], [1, 2]) - 0.5) < 1e-12, auc([1, 2], [1, 2])
    assert auc([], [1]) is None and auc([1], []) is None
    # tie handling: one tie pair straddling the split must give exactly 0.5 there
    assert abs(auc([2], [2]) - 0.5) < 1e-12
    # num(): NaN and junk never become 0.0
    assert num("nan") is None and num("") is None and num(None) is None
    assert num("0.0") == 0.0, "a real zero must survive"
    # permutation p is bounded and never zero
    p = perm_p([3, 4], [1, 2], n=200)
    assert 0 < p <= 1, p
    # MUTATION 1: a label-blind metric must not look like a discriminator.
    rng = random.Random(7)
    noise = [rng.random() for _ in range(150)]
    a = auc(noise[:12], noise[12:])
    assert 0.25 < a < 0.75, f"random metric scored AUC {a}, tie/rank logic is wrong"
    # MUTATION 2: dropping the seed-collapse step must change the answer, i.e. the
    # collapse is load-bearing rather than cosmetic. Unequal seed counts are the case
    # that bites: raw ranking lets a design with more seeds outvote one with fewer.
    # One positive design, 1 seed @0.6. One negative design, 5 seeds @0.9x4 + 0.1.
    pos_seeds, neg_seeds = [0.6], [0.9, 0.9, 0.9, 0.9, 0.1]
    collapsed = auc([statistics.median(pos_seeds)], [statistics.median(neg_seeds)])
    raw = auc(list(pos_seeds), list(neg_seeds))
    assert collapsed == 0.0, collapsed
    assert abs(raw - 0.2) < 1e-12, raw
    assert collapsed != raw, "seed collapse must change the statistic, or it is cosmetic"
    # and the medians it rests on are the per-design values, not the pooled ones
    assert statistics.median(neg_seeds) == 0.9 and statistics.median(pos_seeds) == 0.6
    # REFUSE paths exist for a missing input rather than an empty result
    assert "REFUSE" in open(os.path.abspath(__file__)).read()
    print("calibrate_external.py --selftest PASS")
    return 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else report())
