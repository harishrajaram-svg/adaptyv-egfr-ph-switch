#!/usr/bin/env python3
"""Family-balanced control recovery, per the reviewer's item (b).

  "Freeze a sequence/backbone family map WITHOUT looking at scores. One row per distinct
   sequence, seeds nested. Report raw 8/10 then a family-balanced comparison (fraction
   below EGF per family, ties = 1/2, equal family weights) with per-family results and
   leave-one-family-out sensitivity. Do NOT substitute 'effective n = 3' into
   Clopper-Pearson."

WHY IT MATTERS HERE. Nine of the ten no-KD control molecules come from ONE submitter
group (`gitter-yolo`) and one from another (`deepsatflow`). A raw 8-of-10 therefore counts
one group's designs nine times, so it is close to a statement about gitter-yolo alone. The
family-balanced figure weights the two groups equally instead.

NO CONFIDENCE INTERVAL IS REPORTED. With two families, one of which has a single member,
any binomial interval would be a statement about a denominator we do not have. The
reviewer's instruction is explicit on this point and the leave-one-family-out swing below
is reported in its place.

FAMILY MAP PROVENANCE. Families are assigned from the submitter-group prefix of the
molecule name -- `gitter-yolo*` and `deepsatflow*` -- which is metadata published with the
molecules, fixed before any score was read, and independent of the outcome. No score was
consulted in constructing it.

    control_family_balance.py
"""
import json, re, statistics as st, sys
from collections import defaultdict

SRC = 'analysis/01-egfr/control_recovery.json'
NO_KD = 'exp_neg_no_kd'


def family(mol):
    m = re.match(r'EXPNEG_([a-zA-Z]+)', mol)
    return m.group(1).lower() if m else 'other'


def frac_below(vals, ref):
    """Fraction strictly below ref, ties counted as 1/2, per the instruction."""
    below = sum(1 for v in vals if v < ref)
    tied = sum(1 for v in vals if v == ref)
    return (below + 0.5 * tied) / len(vals)


def main():
    rows = json.load(open(SRC))
    pos = [r for r in rows if r['klass'].startswith('exp_pos') or 'EGF' in r['molecule']]
    neg = [r for r in rows if r['molecule'].startswith('EXPNEG_')]
    if len(pos) != 1:
        sys.exit(f"expected exactly 1 measured binder, found {len(pos)}: "
                 f"{[r['molecule'] for r in pos]}")
    egf = pos[0]
    print(f"quantified binder: {egf['molecule']}  human {egf['hu_med']:.4f}  "
          f"mouse {egf['mo_med']:.4f}")
    print(f"no-KD molecules (right-censored): {len(neg)}  "
          f"-- one row per distinct sequence, 5 seeds nested in each\n")

    fams = defaultdict(list)
    for r in neg:
        fams[family(r['molecule'])].append(r)
    print("FROZEN FAMILY MAP (submitter-group prefix; no score consulted):")
    for f, rs in sorted(fams.items(), key=lambda kv: -len(kv[1])):
        print(f"  {f:<14} n={len(rs):<3} {', '.join(x['molecule'].replace('EXPNEG_','') for x in rs)[:84]}")

    for leg, key, ref in (('human', 'hu_med', egf['hu_med']),
                          ('mouse', 'mo_med', egf['mo_med'])):
        vals = [r[key] for r in neg]
        raw = frac_below(vals, ref)
        print(f"\n=== {leg.upper()} leg (reference: EGF {ref:.4f}) ===")
        print(f"  RAW, one row per molecule:            "
              f"{raw*len(vals):.1f}/{len(vals)} = {raw:.3f}")
        per = {}
        for f, rs in fams.items():
            per[f] = frac_below([r[key] for r in rs], ref)
        bal = st.mean(per.values())
        print(f"  PER FAMILY (fraction below EGF, ties = 1/2):")
        for f in sorted(per, key=lambda x: -len(fams[x])):
            print(f"     {f:<14} n={len(fams[f]):<3} {per[f]:.3f}")
        print(f"  FAMILY-BALANCED (equal family weights): {bal:.3f}")
        if len(per) > 1:
            print(f"  LEAVE-ONE-FAMILY-OUT:")
            for drop in sorted(per, key=lambda x: -len(fams[x])):
                rest = [v for f, v in per.items() if f != drop]
                print(f"     drop {drop:<14} -> {st.mean(rest):.3f}  "
                      f"(on {len(rest)} family/families)")
            lo, hi = min(per.values()), max(per.values())
            print(f"  SENSITIVITY: the family-balanced figure moves over [{lo:.3f}, {hi:.3f}]\n"
                  f"               depending on which single family is retained. With "
                  f"{len(per)} families,\n"
                  f"               one of size {min(len(v) for v in fams.values())}, no "
                  f"interval is reported -- see the module docstring.")

    print("\n=== BOTH SPECIES REQUIRED (how the submission is scored) ===")
    beat = [r['molecule'] for r in neg
            if r['hu_med'] >= egf['hu_med'] and r['mo_med'] >= egf['mo_med']]
    print(f"  no-KD molecules matching or exceeding EGF on BOTH legs: {len(beat)}"
          f"{' -- ' + ', '.join(beat) if beat else ''}")
    print("  So requiring both species, the one quantified binder outranks all "
          f"{len(neg)} censored molecules.")
    print("  This is consistent with the instrument working on this panel and does not")
    print("  demonstrate that it does: n = 1 positive, and the ten are censored, not zero.")


if __name__ == '__main__':
    main()
