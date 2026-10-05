#!/usr/bin/env python3
"""Every shipped design must have a levelled novelty record that clears the gate.

WHY THIS IS A GATE AND NOT A NOTE
Adaptyv's novelty filter is a HARD eligibility rule, not a ranking input: the submission
gate is level >= 3. A shipped design with no levelled record is not "probably fine", it is
an unknown sitting in a graded upload. This was discovered on 2026-10-05 with 4 of 18
shipped designs never levelled -- one of them rank 4, the best-corroborated design in the
submission.

THE JOIN IS THE HARD PART, AND IT IS WHY THE GAP WENT UNNOTICED.
The novelty TSVs are keyed on the PDB path of the pose that was scanned; the CSV is keyed
on the submission name. They do not match and cannot be matched by prefix -- a fuzzy
`name in stem or stem in name` match was tried first and it silently reported a shipped
VHH as Level 1 / does-not-clear by matching the WRONG row. So the map below is EXPLICIT.
A shipped design absent from it fails the gate rather than being guessed at.

LEVEL RULE (bin/novelty_gate.py implements it; restated here only to read the margins)
  HIGH structural TM >= 0.80;  moderate TM >= 0.50;  sequence bars at 70% and 30%.
  L3 = exactly one of (seq > 30%, >= moderate structural)  -> clears
  L4 = seq <= 30% AND less than moderate                   -> clears
  L2 = HIGH structural ALONE is enough to land here        -> DOES NOT CLEAR
So a design at 12-20% identity clears iff its TM stays under 0.80. The whole submission
lives in a band from 0.59 to 0.79, which means the margin to the cliff is the number that
matters, and `cf_short120_r031` at TM 0.811 is the precedent for landing on the wrong side.

    check_novelty_coverage.py            # report + exit non-zero on any unlevelled design
    check_novelty_coverage.py --margins  # also print each design's distance to the cliff
"""
import csv, glob, os, sys

CSV = 'submissions/01-egfr.csv'
HIGH_TM = 0.80

# submission name -> stem of the pose actually scanned, verified one at a time by reading
# the TSV rows. `None` = no record exists anywhere; that is a FAILURE, not an exemption.
STEM = {
    'c5_cf_short__boltzgen_egfr_cropfree_short_48': 'c5_cf_short__boltzgen_egfr_cropfree_short_48',
    'c5_cr_crop_patch__boltzgen_egfr_crop_patch_05': 'c5_cr_crop_patch__boltzgen_egfr_crop_patch_05',
    'rimA01_r15_boltzgen_egfr_d3_rimA_20': 'rimA01_r15_boltzgen_egfr_d3_rimA_20',
    'bc_s360518_mpnn9_A22D': None,
    'ss_bc_s831683_mpnn6_S15D_S62H_routeA': 'ss_bc_s831683_mpnn6_S15D_S62H_routeA',
    'd2c_mpnn13_S88D_serasp': None,
    'cons_gap_h370_only__boltzgen_egfr_h370_018': 'cons_gap_h370_only__boltzgen_egfr_h370_018',
    'bcr_d3acid3_l60_s647537_mpnn3': 'bcr_d3acid3_l60_s647537_mpnn3',
    'bcr_d3acid3_l60_s647537_mpnn11': 'bcr_d3acid3_l60_s647537_mpnn11',
    'bc_s831683_mpnn6_S15D': 'mpnn6_S15D',
    'bc_s831683_mpnn19_S15D': 'mpnn19_S15D',
    'rimA02_d3_rimA_14_vhh': 'rimA02__d3_rimA_14',
    'h370_020_vhh': 'boltzgen_egfr_h370_020',
    'rimA01_r15_L133E': None,
    'bc_s831683_mpnn9_S15D': 'mpnn9_S15D',
    'bc_s831683_mpnn9_WT': 'mpnn9_WT',
    'bc_d3acid_l65_s831683_mpnn11': None,
    'bc_s831683_mpnn8_S15D': 'mpnn8_S15D',
}

# For an unlevelled design, the nearest scanned relative and why it is or is not a safe
# inference. A single point mutation moves sequence identity by ~1/L and TM by very little,
# so a parent with wide margin to HIGH_TM carries its mutant; a parent sitting ON the cliff
# does not. These stems are read from the TSVs, never assumed.
# A measurement that EXISTS but is NOT REPRODUCIBLE is recorded here rather than silently
# credited. rank 4's own assessment string in the graded CSV states "Novelty re-measured on
# the MUTANT pose, not the wild-type backbone: TM 0.792, identity 15.2%, Level 3 -- it
# clears", and git commit 332b09e (2026-10-04 19:53) says the same and adds "clears by
# 0.008, so it is also the row most exposed to a domain-wise novelty rejection". So the run
# happened. Its output was never persisted: no novelty TSV at any commit contains TM 0.792
# for this design, and foldseek plus ~/fsdb are both gone, so it cannot be re-run locally.
# The gate therefore treats rank 4 as UNRESOLVED, not as clearing -- an unreproducible
# number is not evidence in a document where every other number is reproducible.
RECORDED_BUT_UNREPRODUCIBLE = {
    'bc_s360518_mpnn9_A22D': ('TM 0.792, identity 15.2%, Level 3, clears by 0.008',
                              'commit 332b09e and the CSV assessment string; no TSV anywhere'),
}

PROXY = {
    'bc_s360518_mpnn9_A22D': ('d3acid3_l65_s360518_mpnn9_model1', 'single A22D from this parent'),
    'rimA01_r15_L133E': ('rimA01_r15_boltzgen_egfr_d3_rimA_20', 'single L133E from this parent'),
    'bc_d3acid_l65_s831683_mpnn11': ('d3acid_l65_s831683_mpnn11_model1', 'same backbone, MPNN sibling'),
    'd2c_mpnn13_S88D_serasp': ('d2c_101_l147_s144898_m_S88D', 'same lineage and same S88D mutation'),
}


def records():
    """stem -> dict of every novelty row found for it, across all TSVs."""
    out = {}
    for f in sorted(glob.glob('analysis/01-egfr/novelty_*.tsv')):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                stem = os.path.basename(r.get('design', '')).replace('.pdb', '')
                if stem and stem != 'design':
                    r['_src'] = os.path.basename(f)
                    out.setdefault(stem, []).append(r)
    return out


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    rows = list(csv.DictReader(open(CSV)))
    rec = records()
    unmapped = [r['name'] for r in rows if r['name'] not in STEM]
    if unmapped:
        sys.exit("FAIL: no STEM entry for " + ", ".join(unmapped) +
                 "\n  Add it to bin/check_novelty_coverage.py after reading the TSV row.\n"
                 "  There is deliberately no fuzzy fallback: it matched the wrong row once\n"
                 "  and reported a clearing design as Level 1.")

    cleared, unlevelled = [], []
    show = '--margins' in sys.argv
    if show:
        print(f"{'#':>3} {'design':<46}{'lvl':>4}{'gate':>6}{'TM':>8}{'ident':>7}"
              f"  {'margin/rule':>11}")
    for i, r in enumerate(rows, 1):
        stem = STEM[r['name']]
        hit = None
        for h in rec.get(stem, []) if stem else []:
            if h.get('clears_gate') not in (None, '', 'None'):
                hit = h
                break
        if hit:
            cleared.append((i, r['name'], hit))
            if show:
                tm, fid = num(hit.get('best_qtm')), num(hit.get('best_fident'))
                # ANTIBODIES ARE LEVELLED ON A DIFFERENT RULE BRANCH (CDRH3 identity, not
                # whole-chain TM), so the distance to HIGH_TM is not their cliff and
                # printing it is actively misleading: both shipped VHHs sit ABOVE 0.80
                # and clear at level 3 and 4. Report their CDRH3 identity instead.
                ab = str(hit.get('antibody')).strip().lower() == 'true' or hit.get('cdrh3_id')
                if ab:
                    c = num(hit.get('cdrh3_id'))
                    mg = f"CDRH3 {c:.3f}" if c is not None else 'CDRH3 rule'
                else:
                    mg = f"{HIGH_TM - tm:+.4f}" if tm is not None else '—'
                print(f"{i:>3} {r['name'][:45]:<46}{hit.get('level'):>4}"
                      f"{hit.get('clears_gate'):>6}{tm if tm is not None else '—':>8}"
                      f"{fid if fid is not None else '—':>7}  {mg:>11}")
        else:
            unlevelled.append((i, r['name']))

    bad = [(i, n, h) for i, n, h in cleared if h.get('clears_gate') != 'True']
    print(f"\nlevelled and clearing : {len(cleared) - len(bad)} of {len(rows)}")
    if bad:
        print(f"LEVELLED AND FAILING  : {len(bad)}")
        for i, n, h in bad:
            print(f"   rank {i:>2}  {n}  level {h.get('level')}")
    print(f"NOT LEVELLED          : {len(unlevelled)}")

    risky = []
    for i, n in unlevelled:
        px, why = PROXY.get(n, (None, 'no nearest scanned relative identified'))
        h = next((x for x in rec.get(px, []) if num(x.get('best_qtm')) is not None), None) if px else None
        tm = num(h.get('best_qtm')) if h else None
        fid = num(h.get('best_fident')) if h else None
        if tm is None:
            print(f"   rank {i:>2}  {n}\n            UNKNOWN -- {why}")
            risky.append(n)
            continue
        margin = HIGH_TM - tm
        verdict = 'AT RISK' if margin < 0.05 else 'inference safe'
        if margin < 0.05:
            risky.append(n)
        print(f"   rank {i:>2}  {n}\n            proxy {px} TM {tm:.4f} "
              f"ident {fid:.3f} -> margin to HIGH {margin:+.4f}  {verdict}\n"
              f"            ({why})")

    thin = []
    for i, n, h in cleared:
        if str(h.get('antibody')).strip().lower() == 'true' or h.get('cdrh3_id'):
            continue
        tm = num(h.get('best_qtm'))
        if tm is not None and HIGH_TM - tm < 0.04:
            thin.append((i, n, HIGH_TM - tm))
    if thin:
        print(f"\nTHIN MARGIN (<0.04 TM from the level-2 cliff), levelled and clearing "
              f"but with little room: ")
        for i, n, m in sorted(thin, key=lambda t: t[2]):
            print(f"   rank {i:>2}  {n}  {m:+.4f}")
        print("  These already have a verdict; the note is that the verdict is not robust\n"
              "  to a re-scan against a larger database than ours.")
    for n, (val, src) in RECORDED_BUT_UNREPRODUCIBLE.items():
        if any(n == x for _, x in unlevelled):
            print(f"\n   RECORDED BUT UNREPRODUCIBLE -- {n}")
            print(f"     claims: {val}")
            print(f"     source: {src}")
            print(f"     The run happened; its output was not persisted, and foldseek and")
            print(f"     ~/fsdb are absent, so it cannot be reproduced. Treated as")
            print(f"     unresolved rather than clearing.")
    if risky:
        print(f"\nAT RISK, needs an actual levelling run: {', '.join(risky)}")
        print("  Local re-levelling is not possible: foldseek and ~/fsdb are absent.")
        print("  Resolve in the Adaptyv portal's novelty check before nominating.")
    if unlevelled or bad:
        print("\nFAIL: a hard eligibility gate is unresolved for "
              f"{len(unlevelled) + len(bad)} shipped design(s).")
        sys.exit(1)
    print("\nPASS: every shipped design is levelled and clears level >= 3.")


if __name__ == '__main__':
    main()
