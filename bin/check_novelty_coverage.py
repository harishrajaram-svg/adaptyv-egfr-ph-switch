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
    'bc_s360518_mpnn9_A22D': 'bc_s360518_mpnn9_A22D',
    'ss_bc_s831683_mpnn6_S15D_S62H_routeA': 'ss_bc_s831683_mpnn6_S15D_S62H_routeA',
    'd2c_mpnn13_S88D_serasp': 'd2c_mpnn13_S88D_serasp',
    'cons_gap_h370_only__boltzgen_egfr_h370_018': 'cons_gap_h370_only__boltzgen_egfr_h370_018',
    'bcr_d3acid3_l60_s647537_mpnn3': 'bcr_d3acid3_l60_s647537_mpnn3',
    'bcr_d3acid3_l60_s647537_mpnn11': 'bcr_d3acid3_l60_s647537_mpnn11',
    'bc_s831683_mpnn6_S15D': 'mpnn6_S15D',
    'bc_s831683_mpnn19_S15D': 'mpnn19_S15D',
    'rimA02_d3_rimA_14_vhh': 'rimA02__d3_rimA_14',
    'h370_020_vhh': 'boltzgen_egfr_h370_020',
    'rimA01_r15_L133E': 'rimA01_r15_L133E',
    'bc_s831683_mpnn9_S15D': 'mpnn9_S15D',
    'bc_s831683_mpnn9_WT': 'mpnn9_WT',
    'bc_d3acid_l65_s831683_mpnn11': 'bc_d3acid_l65_s831683_mpnn11',
    'bc_s831683_mpnn8_S15D': 'mpnn8_S15D',
}

# For an unlevelled design, the nearest scanned relative and why it is or is not a safe
# inference. A single point mutation moves sequence identity by ~1/L and TM by very little,
# so a parent with wide margin to HIGH_TM carries its mutant; a parent sitting ON the cliff
# does not. These stems are read from the TSVs, never assumed.
# RESOLVED 2026-10-05. Four designs had no levelled record, and rank 4's figure existed
# only in a commit message and the CSV prose: commit 332b09e (2026-10-04 19:53) said
# "Novelty re-measured on the MUTANT pose ... TM 0.792, 15.2% identity, Level 3 -- clears
# by 0.008", and no novelty TSV anywhere contained it. The run had happened; its output was
# never persisted, and the toolchain was gone.
#
# That is a toolchain gap, not an unreproducible result, so it was closed rather than
# declared: foldseek reinstalled from the upstream static binary, the FoldSeek PDB database
# re-downloaded (2.2 GB transfer, 6.4 GB indexed at ~/fsdb), the four binder chains
# re-extracted from their ORIGINAL unrelaxed ESMFold2 poses by SEQUENCE, and
# bin/novelty_gate.py re-run. All four clear Level 3, and rank 4 reproduces at
# **TM 0.792** exactly -- the figure from the commit message, now in
# analysis/01-egfr/novelty_gap4.tsv and committed so it stays reproducible. Its measured
# identity is 0.141 against the 15.2% quoted, which does not change the level (both are
# well under the 30% bar).
#
# The 0.008 margin to the level-2 cliff is REAL and remains the submission's sharpest
# eligibility exposure; see limitation 32.

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

    # Two designs may never share a novelty stem: one scanned chain cannot be two
    # sequences. The nres check alone cannot see this when the designs are the same
    # length -- pointing both 65 aa designs at one 65 aa stem passed it.
    used = {}
    for r in rows:
        st_ = STEM.get(r['name'])
        if st_:
            used.setdefault(st_, []).append(r['name'])
    dupes = {k: v for k, v in used.items() if len(v) > 1}
    if dupes:
        print("STEM MAP IS NOT INJECTIVE -- one scanned chain claimed by several designs:")
        for k, v in dupes.items():
            print(f"   stem {k!r} <- " + ", ".join(v))
        sys.exit("  Each design needs its own novelty record; fix the STEM map.")

    cleared, unlevelled, misjoin = [], [], []
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
            # VALIDATE THE JOIN, do not trust the hand map. Pointing all four previously
            # unmapped designs at a single stem made the gate report "18 of 18 clear" --
            # a hand-maintained map that is never checked is the ELIGIBILITY: 0 failure
            # again. The scanned chain's residue count must match the submitted sequence's
            # length, which is cheap and catches every mis-join of this shape.
            try:
                nres = int(hit.get('nres'))
            except (TypeError, ValueError):
                nres = None
            if nres is not None and nres != len(r['sequence'].strip()):
                misjoin.append((i, r['name'], stem, nres, len(r['sequence'].strip())))
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

    if misjoin:
        print(f"\nJOIN VALIDATION FAILED for {len(misjoin)} design(s) -- the novelty record "
              f"found does not describe the submitted sequence:")
        for i, n, stem, nres, aa in misjoin:
            print(f"   rank {i:>2}  {n}\n            stem {stem!r} has {nres} residues; "
                  f"the submitted sequence is {aa} aa")
        print("  Fix the STEM map in bin/check_novelty_coverage.py. The map is "
              "hand-maintained and this check is what keeps it honest.")
        sys.exit(1)

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
