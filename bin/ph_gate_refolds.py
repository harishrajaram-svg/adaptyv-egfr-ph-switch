#!/usr/bin/env python3
"""Re-gate the pH ratio on ESMFold2 refold poses and compare to the generator pose.

WHY THIS EXISTS. ph_gate_all.py gates every design on its BoltzGen pose, which is the
only pose that exists before any GPU is spent. But measured on six designs gated both
ways, the refold moved the ratio by 0.25x to 1.80x:

    mechA_d4_06        2.64 -> 0.67   (0.25x)
    d3_rimA_26         2.12 -> 0.70   (0.33x)
    cropfree_short_59  1.88 -> 0.94   (0.50x)
    cropfree_short_031 3.55 -> 3.95   (1.11x)
    cropfree_21        1.46 -> 1.98   (1.36x)
    bg04/d3_00         2.39 -> 4.30   (1.80x)  <- became the 2nd strongest switch we have

So a generator-pose verdict of "no switch" is not a verdict. The refold is the
independent prediction and it is the one that counts. This script runs the gate on
refolds and prints both numbers side by side.

Usage:  ph_gate_refolds.py <esmfold2_out_dir> [more dirs...]
"""
import glob, json, os, re, sys, statistics as st
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_gate_all import one          # same free-leg-in-place implementation, unchanged

# DO NOT COMPARE `best_ratio`. It is the max over every target histidine, and a
# histidine the binder is nowhere near returns pKa_bound == pKa_free, i.e. a ratio of
# EXACTLY 1.000. A median of 3 of 5 sites per design sit on that floor, and 1263 of
# 1584 designs report best_ratio == 1.000 for that reason alone. So best_ratio cannot
# go below 1.0 in practice and comparing it across poses compares one floor to another.
# The site that matters is the one the binder actually MOVED.
def active(sites):
    """The histidine that most changes the OBSERVABLE, and its ratio.

    Selected by largest |log(ratio)|, NOT by largest |dpKa|. A histidine can shift three
    pKa units and change nothing: if pKa_free and pKa_bound both sit far BELOW 6.5 the
    site is deprotonated at both pH values and the ratio stays ~1.0. Several rs_crop6
    designs show exactly that -- H370 moving +2.96 and +2.98 units for a ratio of 0.98 --
    and selecting on dpKa reported those designs as flat while a genuinely switching site
    on the same design was ignored. Selecting on |log(ratio)| cannot make that error.
    """
    import math
    cand = [(abs(math.log(v['ratio'])) if v['ratio'] > 0 else 0.0, k, v) for k, v in sites.items()]
    cand.sort(reverse=True)
    if not cand or cand[0][0] < 0.01:          # nothing moved the observable
        return None, None, 0.0
    _, k, v = cand[0]
    return k, v['ratio'], v['bound'] - v['free']

def design_of(path):
    """Strip the ESMFold2 _seed<N>_sample_<N> suffix and the _hu/_mo species tag."""
    b = os.path.basename(path)[:-4]
    import re
    b = re.sub(r'_seed\d+_sample_\d+$', '', b)
    return re.sub(r'_(hu|mo|human|mouse)$', '', b)

if __name__ == '__main__':
    files = []
    for d in sys.argv[1:]:
        # Human leg only, one pose per seed. TWO naming conventions exist in this project:
        # today's runs use `_hu_` / `_mo_`, yesterday's used `_human_` / `_mouse_`. Matching
        # only `_hu_` silently skipped every older run (mimic14, confirm_top2, reval_*).
        files += [p for p in glob.glob(os.path.join(d, '**', '*.cif'), recursive=True)
                  if ('_hu_' in os.path.basename(p) or '_human_' in os.path.basename(p))]
    print(f'{len(files)} human refold poses across {len(sys.argv)-1} dirs', flush=True)
    res, skipped = {}, []
    with ProcessPoolExecutor(max_workers=16) as ex:
        for p, r in zip(files, ex.map(one, files, chunksize=4)):
            if not r:
                skipped.append(p); continue
            site, ratio, d = active(r['sites'])
            if site is None: continue
            res.setdefault(design_of(p), []).append((ratio, site, d))
    if skipped:
        print(f'\n*** {len(skipped)} poses SKIPPED (unrecognised target family or no pKa). '
              f'A silent skip here hid 32 ECD designs once; it is now reported.')
        for p in skipped[:5]: print(f'      {os.path.basename(p)}')
    # generator-pose ratios, keyed by the design name embedded in the faa/run name
    gen = {}
    for r in json.load(open('analysis/01-egfr/ph_gate_all.json')):
        site, ratio, d = active(r['sites'])
        if site is None: continue
        prev = gen.get(r['file'])
        if prev is None or ratio > prev[0]: gen[r['file']] = (ratio, site)
    rows = []
    for d, v in res.items():
        stem = d.split('__')[-1]
        g = gen.get(stem)
        rows.append((st.median([x[0] for x in v]), len(v), g[0] if g else None,
                     v[0][1], g[1] if g else '-', d))
    rows.sort(key=lambda r: -r[0])
    print(f'\n{len(rows)} designs re-gated on refolds\n')
    print(f"{'refold':>8}{'n':>3}{'gen':>7}{'move':>8} {'site(rf/gen)':>14}  design")
    for med, n, g, srf, sgn, d in rows:
        mv = f'{med/g:.2f}x' if g and g > 0 else '   -  '
        gs = f'{g:.2f}' if g else '  -  '
        flag = ''
        if g and med >= 1.20 and g < 1.20: flag = '  *** NEW SWITCH -- generator said no'
        elif g and med < 1.20 and g >= 1.20: flag = '  (generator false positive)'
        print(f"{med:8.2f}{n:3d}{gs:>7}{mv:>8} {srf+'/'+sgn:>14}  {d[:44]}{flag}")
    new = [r for r in rows if r[2] is not None and r[0] >= 1.20 and r[2] < 1.20]
    fp  = [r for r in rows if r[2] is not None and r[0] < 1.20 and r[2] >= 1.20]
    print(f'\nNEW switches the generator pose missed: {len(new)}')
    print(f'generator false positives killed by the refold: {len(fp)}')
    # MERGE, never overwrite. This file was written four times on 2026-10-04 and each run
    # replaced the last, so it held 30 designs after re-gating 1,200. The consolidated
    # record is the thing the methods document cites; a partial one is worse than none.
    OUT = 'analysis/01-egfr/ph_refold_regate.json'
    prev = {}
    if os.path.exists(OUT):
        for r in json.load(open(OUT)):
            prev[r['design']] = r
    for m, n, g, srf, sgn, d in rows:
        prev[d] = {'design': d, 'refold': m, 'n': n, 'generator': g,
                   'site_refold': srf, 'site_gen': sgn}
    json.dump(sorted(prev.values(), key=lambda r: -r['refold']), open(OUT, 'w'), indent=1)
    print(f'consolidated record: {len(prev)} designs in {OUT}')
