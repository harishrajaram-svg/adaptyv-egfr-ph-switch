#!/usr/bin/env python3
"""The project's one ranking table, keyed on the thing we actually submit: the sequence.

WHY
---
pH and affinity have been joined by NAME all project. Names carry the arm that scored
them, so the same molecule is several rows -- `ut_rimA02__d3_rimA_14` and
`rest_rimA01__boltzgen_egfr_d3_rimA_14` are one 129-residue VHH, and reading either row
alone reports n=5 when 6 poses exist and a median 0.36x off the pooled value. Earlier
name joins in this project silently dropped 337 designs from affinity scoring and
mislabelled 83 generator values.

So: pH comes from `ph_pooled_by_sequence.json` (median over every human-leg refold pose
of that exact binder) and affinity is pooled the same way here -- every ESMFold2 pose of
that binder against each species, max and median over poses. The join key is the binder
sequence. Names are kept only as provenance.

ipSAE_min is the minimum over the two ALIGNMENT DIRECTIONS of one interface (the two
`asym` rows), never a minimum across chain pairs.

Output: analysis/01-egfr/master_rank.json  + a printed table.
Ranking follows the stated competition hierarchy: pH selectivity, then mouse
cross-reactivity, then human affinity -- with a floor of MIN_N poses before a pH ratio
is allowed to put a design in tier 1, because one pose is not a measurement.
"""
import glob, json, os, re, statistics as st, collections

MIN_N = 5
RATIO_BAR = 1.20
OUT = 'analysis/01-egfr/master_rank.json'
TARGET_LENS = None   # binder is simply the SHORTER record; target constructs are 170/609/621


def faa_binder(path):
    """Binder sequence from a two-record complex faa: the shorter record."""
    recs, cur = [], []
    for ln in open(path):
        if ln.startswith('>'):
            if cur: recs.append(''.join(cur)); cur = []
        else:
            cur.append(ln.strip())
    if cur: recs.append(''.join(cur))
    return min(recs, key=len) if len(recs) >= 2 else (recs[0] if recs else None)


def ipsae_min(txt):
    vals = [float(q[5]) for q in (ln.split() for ln in open(txt))
            if len(q) > 5 and q[4] == 'asym']
    if len(vals) >= 2: return min(vals)
    return vals[0] if vals else None


def main():
    # index every scored pose directory by its basename, once
    posedirs = {}
    for d in glob.glob('runs/esmfold2/*/**/', recursive=True):
        posedirs.setdefault(os.path.basename(d.rstrip('/')), []).append(d)

    aff = collections.defaultdict(lambda: collections.defaultdict(list))
    names = collections.defaultdict(set)
    for faa in glob.glob('analysis/01-egfr/score_*/*.faa'):
        stem = os.path.basename(faa)[:-4]
        m = re.match(r'(.+)_(hu|mo|human|mouse)$', stem)
        if not m: continue
        base, sp = m.group(1), ('hu' if m.group(2) in ('hu', 'human') else 'mo')
        seq = faa_binder(faa)
        if not seq: continue
        names[seq].add(base)
        for d in posedirs.get(stem, []):
            for t in glob.glob(os.path.join(d, '*_10_10.txt')):
                v = ipsae_min(t)
                if v is not None: aff[seq][sp].append(v)

    ph = {r['seq']: r for r in json.load(open('analysis/01-egfr/ph_pooled_by_sequence.json'))}
    rows = []
    for seq in set(ph) | set(aff):
        p = ph.get(seq)
        a = aff.get(seq, {})
        f = lambda sp, g: (g(a[sp]) if a.get(sp) else None)
        rows.append(dict(
            seq=seq, aa=len(seq),
            ratio=p['median'] if p else None, ratio_max=p['max'] if p else None,
            ratio_n=p['n'] if p else 0, site=p['site'] if p else None,
            hu=f('hu', max), hu_med=f('hu', st.median), hu_n=len(a.get('hu', [])),
            mo=f('mo', max), mo_med=f('mo', st.median), mo_n=len(a.get('mo', [])),
            names=sorted((p['names'] if p else []) + sorted(names.get(seq, [])))[:6]))

    def key(r):
        tier1 = (r['ratio'] or 0) >= RATIO_BAR and r['ratio_n'] >= MIN_N
        return (0 if tier1 else 1, -(r['ratio'] or 0) if tier1 else 0,
                -(r['mo'] or 0), -(r['hu'] or 0))
    rows.sort(key=key)
    json.dump(rows, open(OUT, 'w'), indent=1)

    print(f"{len(rows)} distinct binder sequences; wrote {OUT}\n")
    hdr = f"{'#':>3} {'pH':>6}{'n':>4} {'mo_max':>7}{'mo_med':>7} {'hu_max':>7}{'hu_med':>7} {'aa':>4} {'site':>5}  name"
    print(hdr); print('-' * len(hdr))
    shown = 0
    for i, r in enumerate(rows, 1):
        if r['ratio_n'] < MIN_N or (r['ratio'] or 0) < 2.0: continue
        shown += 1
        if shown > 30: break
        g = lambda v: f"{v:7.4f}" if v is not None else "      -"
        print(f"{shown:>3} {r['ratio']:6.2f}{r['ratio_n']:4d} {g(r['mo'])}{g(r['mo_med'])} "
              f"{g(r['hu'])}{g(r['hu_med'])} {r['aa']:4d} {str(r['site']):>5}  {r['names'][0][:44] if r['names'] else '?'}")
    print(f"\n(showing designs with ratio >= 2.0 confirmed on >= {MIN_N} poses)")
    unscored = [r for r in rows if r['ratio_n'] >= MIN_N and (r['ratio'] or 0) >= 2.0
                and not (r['hu_n'] or r['mo_n'])]
    print(f"{len(unscored)} of those have NO affinity measurement at all")
    for r in unscored[:10]:
        print(f"    {r['ratio']:5.2f} n={r['ratio_n']:<3} {r['names'][0][:54] if r['names'] else '?'}")


main()
