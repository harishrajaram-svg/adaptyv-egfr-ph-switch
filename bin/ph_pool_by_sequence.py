#!/usr/bin/env python3
"""One pH ratio per MOLECULE, pooled over every refold pose of that exact binder sequence.

WHY THIS FILE EXISTS
--------------------
`ph_gate_refolds.py` keys its output on the RUN name, and the same molecule appears
under a different run name in every arm that ever scored it:

    gap_h370_only__boltzgen_egfr_h370_020     n=1   3.104
    cons_gap_h370_only__boltzgen_egfr_h370_020 n=5  2.289
    gc_h370_only__boltzgen_egfr_h370_020       n=5  2.211   <- same molecule, 3 rows

So a design that has 11 poses on disk could read n=1, and the tier-1 floor of 5 poses
would reject it on a bookkeeping artifact. Worse, picking the lowest-n row makes the
OPTIMISTIC single-pose number the headline: 3.104 instead of the honest ~2.25.

The submission's `ratio`/`ratio_n` were produced this way once, by hand, in a throwaway
script. This is that computation written down: the binder SEQUENCE is the identity of a
molecule, and every human-leg pose of it anywhere in runs/esmfold2 counts toward one
median. Names are provenance; sequences are the thing being submitted.

Output: analysis/01-egfr/ph_pooled_by_sequence.json
    [{seq, n, median, max, min, site, poses_by_arm{...}, names[...]}]
Median, not max, is the headline -- the reviewer's v2. Max is kept beside it so the gap between
the two is visible rather than chosen.

Usage:  ph_pool_by_sequence.py [run_dir ...]      (default: every runs/esmfold2/*/)
"""
import glob, json, os, re, statistics as st, sys
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_gate_all import one, family, orient
from ph_gate_refolds import active, design_of
import gemmi

AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def binder_seq(cif):
    """One-letter sequence of the chain that is NOT the target.

    Which chain that is cannot be assumed: 29 of 69 run directories put the target in
    chain A. `orient()` finds the target by its histidine fingerprint and the binder is
    whatever is left -- the same rule the pH gate itself uses, so the two never disagree.
    """
    try:
        st_ = gemmi.read_structure(cif)
        st_.setup_entities(); st_.remove_ligands_and_waters(); st_.setup_entities()
        tgt, fam = orient(st_[0])
        if fam is None:
            return None
        for ch in st_[0]:
            if ch.name != tgt:
                return "".join(AA3to1.get(r.name.upper(), "X") for r in ch)
    except Exception:
        return None
    return None


def job(cif):
    r = one(cif)
    if not r:
        return None
    site, ratio, d = active(r['sites'])
    if site is None:
        return None
    s = binder_seq(cif)
    if not s:
        return None
    return (s, ratio, site, cif)


if __name__ == '__main__':
    dirs = sys.argv[1:] or sorted(glob.glob('runs/esmfold2/*/'))
    files = []
    for d in dirs:
        files += [p for p in glob.glob(os.path.join(d, '**', '*.cif'), recursive=True)
                  if ('_hu_' in os.path.basename(p) or '_human_' in os.path.basename(p))]
    print(f'{len(files)} human refold poses across {len(dirs)} dirs', flush=True)
    pool = {}
    with ProcessPoolExecutor(max_workers=16) as ex:
        for res in ex.map(job, files, chunksize=4):
            if not res:
                continue
            s, ratio, site, cif = res
            e = pool.setdefault(s, dict(seq=s, ratios=[], sites=[], names=set(), arms=set()))
            e['ratios'].append(ratio); e['sites'].append(site)
            e['names'].add(design_of(cif))
            e['arms'].add(cif.split('/')[2])
    rows = []
    for e in pool.values():
        rs = e['ratios']
        rows.append(dict(seq=e['seq'], n=len(rs), median=round(st.median(rs), 3),
                         max=round(max(rs), 3), min=round(min(rs), 3),
                         site=max(set(e['sites']), key=e['sites'].count),
                         n_arms=len(e['arms']), arms=sorted(e['arms']),
                         names=sorted(e['names'])))
    rows.sort(key=lambda r: -r['median'])
    OUT = 'analysis/01-egfr/ph_pooled_by_sequence.json'
    json.dump(rows, open(OUT, 'w'), indent=1)
    print(f'{len(rows)} distinct binder sequences; wrote {OUT}\n')
    print(f"{'median':>7}{'max':>7}{'n':>4}{'arms':>5} {'site':>5}  names")
    for r in rows[:50]:
        print(f"{r['median']:7.2f}{r['max']:7.2f}{r['n']:4d}{r['n_arms']:5d} {r['site']:>5}  "
              f"{', '.join(n[:34] for n in r['names'][:2])}")
    split = [r for r in rows if len(r['names']) > 1]
    print(f"\n{len(split)} molecules appear under MORE THAN ONE run name "
          f"(the fragmentation this script exists to undo)")
    thin = [r for r in rows if r['n'] < 5 and r['median'] >= 2.0]
    print(f"{len(thin)} molecules with median >= 2.0 are still measured on < 5 poses")
