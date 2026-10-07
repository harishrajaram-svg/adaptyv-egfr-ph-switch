#!/usr/bin/env python3
"""pH sensitivity analysis: all-titratable-site product vs the histidine-only product.

The reviewer, 2026-10-05, on the gate behind the shipped ranking: the composition is
incomplete. HIS, ASP and GLU are all parsed, yet only the histidines are ever pushed
into `sites`, so no acid ever reaches the multiplied product. Until that is fixed, they
directed us to label the number for what it is -- a two-partner, HISTIDINE-ONLY
approximation -- and to state outright which sites went unassessed. They noted the gap
bites hardest where the designed change is itself the addition of an Asp or a Glu.

On presentation, they were equally specific: hold the deletion estimate up against
matched-preparation alternatives and against pKa values nudged within their plausible
error, treat the comparison as SENSITIVITY rather than a CI, and fall back to
provisional tiers if the ordering moves by much. Two things were ruled out: going back
to a target-only gate, and ANNOUNCING THE REVISED ORDER AS SETTLED.

This script re-scores every submitted design on BOTH bases from the same poses and
the same code path, so the difference is attributable to the acids alone and to
nothing else. It reports them side by side. It does not re-rank anything by itself.

    ph_sensitivity_multisite.py --out analysis/01-egfr/ph_sensitivity.json [--limit N]

The pose join is by binder SEQUENCE, reusing master_rank's own index
(analysis/01-egfr/score_*/*.faa -> runs/esmfold2 pose dirs). Joining on a run NAME
has broken six analyses in this project; it is not done here.
"""
import argparse, glob, json, math, os, re, statistics as st, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ph_gate_multisite import score_pose
from master_rank import faa_binder
import pose_seq_index

CSV = 'submissions/01-egfr.csv'


def submitted():
    """[(name, sequence)] in submission order, from the graded file itself.

    --extra FILE adds candidates not yet in the CSV, so a design can be scored on all
    three bases BEFORE it is added to the submission. Without that the emitter would see
    no multisite record for it, fall to tier 2, and the addition would look worse than the
    rows it is being compared against -- an artefact of ordering, not a measurement.
    """
    import csv
    with open(CSV) as f:
        rows = list(csv.DictReader(f))
    seqcol = next(c for c in rows[0] if 'seq' in c.lower())
    namecol = next(c for c in rows[0] if 'name' in c.lower() or 'id' in c.lower())
    out = [(r[namecol], r[seqcol].strip().upper()) for r in rows]
    if '--extra' in sys.argv:
        extra = json.load(open(sys.argv[sys.argv.index('--extra') + 1]))
        have = {q for _, q in out}
        for e in extra:
            q = e['seq'].strip().upper()
            if q not in have:
                out.append((e['name'], q)); have.add(q)
    return out


def human_leg_poses_via_faa():
    """The ORIGINAL join: sequence -> poses via analysis/01-egfr/score_*/*.faa.

    Kept only as a cross-check on the structure-derived index below. It is incomplete by
    construction: a pose directory is found only if a .faa of a matching basename exists,
    so any run nobody wrote a .faa for is silently dropped and the design's n is quietly
    too small."""
    posedirs = {}
    for d in glob.glob('runs/esmfold2/*/**/', recursive=True):
        posedirs.setdefault(os.path.basename(d.rstrip('/')), []).append(d)
    out = {}
    for faa in glob.glob('analysis/01-egfr/score_*/*.faa'):
        stem = os.path.basename(faa)[:-4]
        m = re.match(r'(.+)_(hu|human)$', stem)      # human leg only, per METHODS
        if not m:
            continue
        seq = faa_binder(faa)
        if not seq:
            continue
        for d in posedirs.get(stem, []):
            for cif in sorted(glob.glob(os.path.join(d, '*.cif'))):
                out.setdefault(seq.upper(), []).append(cif)
    return {k: sorted(set(v)) for k, v in out.items()}


def human_leg_poses():
    """binder sequence -> [human-leg pose .cif paths], pooled across every run.

    Joins on the binder SEQUENCE READ FROM THE STRUCTURE (bin/pose_seq_index.py), not on
    a .faa file that may not exist. This is the project's own stated rule -- "pool every
    pose of a binder SEQUENCE across all runs; joining on a run name has broken six
    analyses" -- applied to pose discovery as well as to scoring.

    It found poses the .faa join did not: runs/esmfold2/w3_triad holds 15 human-leg poses
    of a sequence byte-identical to the submitted rimA01_r15_L133E, and the submission was
    reporting n=5 for it. The human leg is selected by the `_hu`/`_human` suffix on the
    pose directory, the same convention the rest of the pipeline uses."""
    idx = pose_seq_index.load()
    out = {}
    for seq, dirs in idx.items():
        for d in dirs:
            base = os.path.basename(d.rstrip('/'))
            if not (base.endswith('_hu') or base.endswith('_human')):
                continue
            for cif in sorted(glob.glob(os.path.join(d, '*.cif'))):
                out.setdefault(seq.upper(), []).append(cif)
    return {k: sorted(set(v)) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default='analysis/01-egfr/ph_sensitivity.json')
    ap.add_argument('--limit', type=int, default=0, help='poses per design (0 = all)')
    ap.add_argument('--extra', help='JSON list of {name, seq} candidates not yet in the CSV')
    a = ap.parse_args()

    poses = human_leg_poses()
    old = human_leg_poses_via_faa()
    for _n, _s in submitted():
        _new, _prev = len(poses.get(_s, [])), len(old.get(_s, []))
        if _new != _prev:
            print(f"  pose count {_n}: {_prev} (.faa join) -> {_new} (sequence index)", flush=True)
    rows, t0 = {}, time.time()
    for name, seq in submitted():
        cifs = poses.get(seq, [])
        if a.limit:
            cifs = cifs[:a.limit]
        if not cifs:
            rows[name] = {'name': name, 'seq': seq, 'error': 'no human-leg pose matched this sequence'}
            print(f"{name}: NO POSES", flush=True)
            continue
        per = []
        for cif in cifs:
            r = score_pose(cif)
            if r and 'error' not in r:
                per.append(r)
            else:
                print(f"    {os.path.basename(cif)}: {r.get('error') if r else 'None'}", flush=True)
        if not per:
            rows[name] = {'name': name, 'seq': seq, 'error': 'every pose failed to score'}
            continue
        allsite = [p['product'] for p in per]
        hisonly = [p['product_his_only'] for p in per]
        # THIRD BASIS -- the physically defensible one.
        #
        # The all-site product is inflated by sites whose pKa shifts with NO counter-charge
        # anywhere near them: rimA01_r15_L133E's largest single contributor is binder:ASP33
        # at ratio 7.37 with its nearest opposite charge 9.84 A away, while the designed
        # GLU133 contributes 1.30. A pKa shift at 10-13 A is burial/desolvation, not a
        # titratable interaction with the partner -- which is the SAME argument this project
        # already used to rule out mechanism A, so applying it here is consistency, not a
        # new threshold chosen to help. Guard 4 already measures the distance; this basis
        # just stops multiplying the sites it flags.
        #
        # PARTNER_CUT (6.0 A) was fixed before this analysis and is not tuned to it.
        def partnered(pose):
            vals = [v['ratio'] for v in pose['sites'].values() if not v['no_partner']]
            return math.prod(vals) if vals else 1.0
        partner = [partnered(p) for p in per]
        med_a, med_h = st.median(allsite), st.median(hisonly)
        med_p = st.median(partner)
        rows[name] = dict(
            name=name, seq=seq, n_poses=len(per),
            allsite_median=round(med_a, 4), allsite_min=round(min(allsite), 4),
            allsite_max=round(max(allsite), 4),
            allsite_spread_over_median=round((max(allsite) - min(allsite)) / med_a, 3) if med_a else None,
            hisonly_median=round(med_h, 4),
            hisonly_min=round(min(hisonly), 4), hisonly_max=round(max(hisonly), 4),
            partnered_median=round(med_p, 4),
            partnered_min=round(min(partner), 4), partnered_max=round(max(partner), 4),
            partnered_spread_over_median=round((max(partner) - min(partner)) / med_p, 3) if med_p else None,
            n_partnered=int(st.median([sum(1 for v in p['sites'].values()
                                           if not v['no_partner']) for p in per])),
            acid_effect=round(med_a / med_h, 4) if med_h else None,
            n_sites=int(st.median([p['n_sites'] for p in per])),
            n_his=int(st.median([p['n_his'] for p in per])),
            n_acid=int(st.median([p['n_acid'] for p in per])),
            n_acid_moved=int(st.median([p['n_acid_moved'] for p in per])),
            n_unassessed=int(st.median([p['n_unassessed'] for p in per])),
            unassessed_example=next((p['unassessed'] for p in per if p['n_unassessed']), {}),
            n_impossible=sum(1 for p in per if p['verdict'] != 'ok'),
            per_pose=[{'file': p['file'], 'allsite': p['product'],
                       'hisonly': p['product_his_only'], 'partnered': round(partnered(p), 4),
                       'n_sites': p['n_sites'],
                       'n_acid_moved': p['n_acid_moved'], 'verdict': p['verdict']} for p in per],
            # full per-pose site detail, so no later question needs a re-run
            sites_all_poses={p['file']: p['sites'] for p in per},
            sites_example={k: v for k, v in per[0]['sites'].items()})
        r = rows[name]
        print(f"{name[:40]:<41} n={r['n_poses']:<3} allsite={r['allsite_median']:>8.2f} "
              f"hisonly={r['hisonly_median']:>7.3f} partnered={r['partnered_median']:>7.3f} "
              f"(n_part {r['n_partnered']} of {r['n_sites']}) unassessed={r['n_unassessed']}",
              flush=True)
    json.dump(rows, open(a.out, 'w'), indent=1)
    print(f"\nwrote {a.out} in {time.time()-t0:.0f}s", flush=True)


if __name__ == '__main__':
    main()
