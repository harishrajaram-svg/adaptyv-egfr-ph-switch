#!/usr/bin/env python3
"""Does `his_cation_gate` discriminate where `dddG_elec` does not?

THE QUESTION THIS EXISTS TO ANSWER. s17 and s18 are converging on dddG_elec having no signal on
SINGLE-POINT histidine substitutions -- Adafre 2/3 with the single-His miss, Pertuzumab 5/9 below
a 6/9 null, bH1 5/8 at the null with ZERO of three hits found. If that holds, the ranking has to
rest on something else, and the obvious candidate is `his_cation_gate`: binary, geometric, no
energy term in it, and the only component of this instrument that has never failed a control.

\U0001F534 BUT IT HAS NEVER BEEN TESTED AGAINST A MEASURED SERIES EITHER. Its record is:
  s6c  0 PASS / 9 FAIL / 3 REFUSE across 12 natural TNF complexes -- a mechanism-AVAILABILITY
       measurement, explicitly NOT a selection yield, and at 0% it is not shown to be a filter.
  s6d  its first and only real-coordinate POSITIVE, from a synthetic rigid-body placement.
So its entire positive evidence is one constructed pose. Promoting it to primary ranking key on
that basis would repeat exactly the mistake s14 warns about -- adopting a gate whose pass rate on
real data nobody has measured.

The SIpHAB benchmark fixes that, and for free: 279 single-point variants with hit/no-hit labels,
four campaigns with deposited complexes, and the gate costs under a second per pose on CPU.

WHAT IS MEASURED. For each variant the gate is asked whether the INSTALLED histidine accepts a
hydrogen bond from an Arg or Lys across the interface -- Ahn's second criterion, heavy-atom
geometry only, the `designable` direction. That is a strictly geometric question and it makes no
reference to pH, energy or protonation state.

⚠️ WHAT A PASS WOULD AND WOULD NOT MEAN. The gate was built for OUR mechanism: a binder
histidine against a TARGET cation. These campaigns are antibody CDR scans against four unrelated
antigens, and their switches need not run through that mechanism at all -- Schroter's do not
(s6b: no cation within reach of the installed histidines). So a LOW hit rate here is weak
evidence against the gate, while a hit rate that TRACKS the labels would be strong evidence for
it. The asymmetry is the point and it is why this is worth running even though the targets differ.

Run after campaign_runner.py, which prepares and caches the cleaned structures.
"""
import collections, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))
import his_cation_gate as gate                       # noqa: E402
import campaign_runner as cr                         # noqa: E402


def main():
    bench = json.load(open(cr.BENCH))
    keys = sys.argv[1:] or list(cr.CAMPAIGNS)
    pooled = []
    for name in keys:
        if name not in cr.CAMPAIGNS:
            print(f"unknown campaign {name!r}"); continue
        pid, chains, antigen = cr.CAMPAIGNS[name]
        binder = "".join(sorted(set(chains.values())))
        path, renum, orig = cr.prepare(pid, chains, antigen)
        mapped, dropped, _ = cr.map_positions(bench[name]["rows"], chains, orig)
        print(f"\n### {name}  {pid}   {len(mapped)} mappable variants, "
              f"{sum(1 for m in mapped if m[5])} hits"
              + (f"   ({len(dropped)} dropped)" if dropped else ""))
        rows = []
        for pos, ch, num, ic, wt, hit in mapped:
            vp = cr.build(path, ch, renum[(ch, num, ic)], wt, f"gate_{pid}_{pos}")
            try:
                atoms, _chains = gate.load(vp)   # load() returns (atoms, chain summary)
                hits = gate.scan(atoms, list(binder), list(antigen))
            except Exception as e:
                print(f"   {pos:<8} gate error: {type(e).__name__} {str(e)[:60]}")
                os.remove(vp); continue
            # scan() returns a list of dicts; only the 'designable' direction counts -- a
            # histidine on the BINDER accepting from a cation on the ANTIGEN. 'target-his'
            # is the antigen's own histidine and is not ours to place (s6c).
            n_des = sum(1 for h in hits if h.get("direction") == "designable")
            rows.append((pos, hit, n_des))
            pooled.append((hit, n_des > 0))
            os.remove(vp)
        if not rows:
            continue
        fired = [r for r in rows if r[2] > 0]
        print(f"   gate fired on {len(fired)} of {len(rows)} variants")
        if fired:
            for pos, hit, n in fired:
                print(f"      {pos:<8}{'HIT' if hit else '-':<6}designable={n}")
        tp = sum(1 for _, h, n in rows if h and n > 0)
        fp = sum(1 for _, h, n in rows if not h and n > 0)
        hits = sum(1 for _, h, _ in rows if h)
        print(f"   of {hits} labelled hits the gate fired on {tp}; "
              f"of {len(rows)-hits} non-hits it fired on {fp}")
    if pooled:
        n = len(pooled); hits = sum(1 for h, _ in pooled if h)
        tp = sum(1 for h, f in pooled if h and f); fp = sum(1 for h, f in pooled if not h and f)
        print(f"\n{'='*80}\nPOOLED: {n} variants, {hits} hits. Gate fired on {tp+fp} "
              f"({tp} hits, {fp} non-hits).")
        if tp + fp == 0:
            print("Gate fired on NOTHING. At a 0% rate it is not shown to be a filter (s14) --")
            print("the same verdict s6c reached on natural complexes, now on labelled data.")
        else:
            prec = tp / (tp + fp)
            base = hits / n
            print(f"precision {prec:.1%} against a {base:.1%} base rate "
                  f"-> {'ENRICHED' if prec > base else 'NOT enriched'}")


if __name__ == "__main__":
    main()
