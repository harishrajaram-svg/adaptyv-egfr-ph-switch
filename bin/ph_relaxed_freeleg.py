#!/usr/bin/env python3
"""Relaxed free leg: the third rung of the reviewer's free-state ladder.

The reviewer's position, 2026-10-05, in our own words: predicting the unbound chain on
its own yields one more approximation and carries no guarantee of being the right one.
Their instruction was to line the deletion number up against apo and relaxed numbers
built under a matched preparation, and to read the resulting spread as SENSITIVITY --
never as a measured confidence interval.

Three ways to estimate the FREE-state pKa, with the BOUND leg held identical in all three:

  deletion  the bound complex with the partner's atoms deleted in place. Conformation is
            the bound one; only burial changes. (shipped basis)
  relaxed   the same chain, extracted from the same bound pose, pdbfixer-repaired, then
            energy-minimised under amber14/ff14SB + GBn2 with a 10 kcal/mol/A^2 harmonic
            restraint on N/CA/C/O. Side chains relax; the backbone stays near the bound
            pose. (runs/relaxed/w6)
  apo       the same chain folded ALONE by the same predictor, same 5 seeds.
            (runs/esmfold2/w5_apo, see ph_apo_freeleg.py)

The ladder is ordered by how far the free state is allowed to move from the bound pose:
deletion (not at all) -> relaxed (side chains only) -> apo (everything). It is a
sensitivity analysis over that one axis, not a confidence interval.

Unlike the apo arm, the relaxed arm also covers the TARGET leg for three designs
(31 poses), so the target-side free pKa can be moved too -- the bounded target-leg
comparison the reviewer asked for.

    ph_relaxed_freeleg.py [--out analysis/01-egfr/ph_relaxed_freeleg.json]
"""
import csv, glob, json, os, re, statistics as st, subprocess, sys, tempfile

RELAX = 'runs/relaxed/w6'
SENS = 'analysis/01-egfr/ph_sensitivity.json'
APO = 'analysis/01-egfr/ph_apo_freeleg.json'
CSV = 'submissions/01-egfr.csv'
OUT = 'analysis/01-egfr/ph_relaxed_freeleg.json'
PH_LO, PH_HI = 6.5, 7.4


def link(free, bound):
    def K(ph):
        return (1 + 10 ** (bound - ph)) / (1 + 10 ** (free - ph))
    return K(PH_LO) / K(PH_HI)


def propka_his(pdb_path):
    """pKa of every HIS in an already-written single-chain PDB. {resnum: pKa}.

    Fails closed: a propka crash or an unparseable .pka returns {} and the caller
    falls back to the deletion pKa rather than silently inventing one.
    """
    tag = os.path.basename(pdb_path)[:-4]
    with tempfile.TemporaryDirectory() as td:
        local = os.path.join(td, tag + '.pdb')
        with open(pdb_path) as fi, open(local, 'w') as fo:
            fo.write(fi.read())
        r = subprocess.run([sys.executable, '-m', 'propka', tag + '.pdb'],
                           capture_output=True, cwd=td)
        f = os.path.join(td, tag + '.pka')
        if r.returncode != 0 or not os.path.exists(f):
            return {}
        out = {}
        for ln in open(f):
            q = ln.split()
            if len(q) > 3 and q[0] == 'HIS':
                try:
                    out[int(q[1])] = float(q[3])
                except ValueError:
                    pass
        return out


def main():
    out_path = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else OUT
    sens = json.load(open(SENS))
    ship = {r['sequence'].strip().upper() for r in csv.DictReader(open(CSV))}
    by_name = {v['name']: v for v in sens.values()
               if v.get('seq', '').strip().upper() in ship}
    if not by_name:
        sys.exit("no submitted designs in ph_sensitivity.json")

    # Index the relaxed structures: design -> pose -> {binder,target} -> path.
    # Design names themselves contain '__' (c5_cf_short__boltzgen_egfr_...), so the
    # design/pose boundary CANNOT be found by splitting on '__'. Split by testing the
    # known design names as prefixes, longest first, and refuse to guess otherwise.
    names = sorted(by_name, key=len, reverse=True)
    idx, orphan = {}, []
    for p in sorted(glob.glob(os.path.join(RELAX, '*.pdb'))):
        b = os.path.basename(p)
        m = re.match(r'(.+)__(binder|target)\.pdb$', b)
        if not m:
            orphan.append(b); continue
        stem, leg = m.group(1), m.group(2)
        for n in names:
            if stem.startswith(n + '__'):
                idx.setdefault(n, {}).setdefault(stem[len(n) + 2:], {})[leg] = p
                break
        else:
            orphan.append(b)
    if orphan:
        print(f"WARNING: {len(orphan)} relaxed file(s) match no shipped design, "
              f"e.g. {orphan[0]}")

    cache = {}

    def pkas(path):
        if path not in cache:
            cache[path] = propka_his(path)
        return cache[path]

    rows = []
    hdr = (f"{'design':<42}{'deletion':>10}{'relaxed':>9}{'fold':>7}"
           f"{'bHis':>6}{'poses':>7}{'tgt':>5}")
    print(hdr); print('-' * len(hdr))
    for name, v in by_name.items():
        poses = idx.get(name, {})
        del_m, rel_m, relb_m = [], [], []
        nbind = ntgt = 0
        have = set()
        for pose, sites in (v.get('sites_all_poses') or {}).items():
            paths = poses.get(pose, {})
            if not paths:
                continue
            have.add(pose)
            rb = pkas(paths['binder']) if 'binder' in paths else {}
            rt = pkas(paths['target']) if 'target' in paths else {}
            if not rb and not rt:
                # The structure relaxed fine; it simply carries no histidine on the
                # relaxed chain(s). NOT missing data -- the switch is target-borne and a
                # binder-only relaxation cannot move it. Recorded, excluded from stats.
                continue
            pd = pr = prb = 1.0
            for s in sites.values():
                if s['resname'] != 'HIS':
                    continue
                pd *= link(s['free'], s['bound'])
                # binder-only relaxation
                fb = rb.get(s['resnum']) if s['partner'] == 'binder' else None
                prb *= link(fb if fb is not None else s['free'], s['bound'])
                # both legs relaxed where available
                src = rb if s['partner'] == 'binder' else rt
                f2 = src.get(s['resnum'])
                pr *= link(f2 if f2 is not None else s['free'], s['bound'])
            nbind = sum(1 for s in sites.values() if s['resname'] == 'HIS'
                        and s['partner'] == 'binder' and s['resnum'] in rb)
            ntgt = sum(1 for s in sites.values() if s['resname'] == 'HIS'
                       and s['partner'] == 'target' and s['resnum'] in rt)
            del_m.append(pd); rel_m.append(pr); relb_m.append(prb)
        if not del_m:
            if have:
                nhis = sum(1 for st_ in (v.get('sites_all_poses') or {}).values()
                           for s in st_.values() if s['resname'] == 'HIS'
                           and s['partner'] == 'binder')
                why = ('relaxed, no movable site: 0 binder histidines and no relaxed '
                       'target -- this design switches on the target'
                       if not nhis else
                       'relaxed, but no binder histidine survived the relaxed chain')
                rows.append(dict(name=name, status='no_movable_site', reason=why,
                                 n_poses_relaxed=len(have), deletion=None))
                print(f"{name[:41]:<42}  {why[:44]}")
            else:
                rows.append(dict(name=name, status='no_relaxed_structure'))
                print(f"{name[:41]:<42}  no relaxed structure on disk")
            continue
        d = st.median(del_m); r_ = st.median(rel_m); rb_ = st.median(relb_m)
        rows.append(dict(name=name, deletion=round(d, 4),
                         relaxed_binder_only=round(rb_, 4), relaxed=round(r_, 4),
                         fold=round(r_ / d, 4) if d else None,
                         n_binder_his_relaxed=nbind, n_target_his_relaxed=ntgt,
                         n_poses=len(del_m)))
        print(f"{name[:41]:<42}{d:>10.3f}{r_:>9.3f}{r_/d:>7.2f}"
              f"{nbind:>6}{len(del_m):>7}{ntgt:>5}")

    ok = [x for x in rows if x.get('fold')]
    nomove = [x for x in rows if x.get('status') == 'no_movable_site']
    nostruct = [x for x in rows if x.get('status') == 'no_relaxed_structure']
    summary = dict(n_compared=len(ok), n_no_movable_site=len(nomove),
                   n_no_relaxed_structure=len(nostruct))
    print(f"\nCOVERAGE of the {len(rows)} shipped designs: {len(ok)} compared, "
          f"{len(nomove)} relaxed with no movable site, "
          f"{len(nostruct)} with no relaxed structure.")
    if ok:
        folds = [x['fold'] for x in ok]
        moved = [x for x in ok if abs(x['fold'] - 1.0) >= 0.10]
        print(f"\n{len(ok)} designs compared, {sum(x['n_poses'] for x in ok)} poses.")
        print(f"relaxed/deletion fold: median {st.median(folds):.3f}, "
              f"range {min(folds):.3f}-{max(folds):.3f}")
        print(f"{len(moved)} design(s) move by >=10%"
              + (": " + ", ".join(f"{x['name'][:28]} {x['fold']:.2f}x" for x in moved)
                 if moved else ""))
        import itertools
        do = sorted(ok, key=lambda x: -x['deletion'])
        ro = sorted(ok, key=lambda x: -x['relaxed'])
        di = {x['name']: i for i, x in enumerate(do)}
        ri = {x['name']: i for i, x in enumerate(ro)}
        c = dc = 0
        for x, y in itertools.combinations(ok, 2):
            if (di[x['name']] - di[y['name']]) * (ri[x['name']] - ri[y['name']]) > 0:
                c += 1
            else:
                dc += 1
        tau = (c - dc) / (c + dc)
        shift = max(abs(di[x['name']] - ri[x['name']]) for x in ok)
        print(f"Kendall tau deletion vs relaxed ordering: {tau:+.3f}; "
              f"largest single-design rank shift: {shift}")
        summary.update(n=len(ok), fold_median=round(st.median(folds), 4),
                       fold_min=round(min(folds), 4), fold_max=round(max(folds), 4),
                       n_moved_10pct=len(moved), kendall_tau=round(tau, 4),
                       max_rank_shift=shift)

        # three-rung ladder, where the apo arm also has a number
        if os.path.exists(APO):
            apo = {x['name']: x for x in json.load(open(APO))['rows'] if 'apo' in x}
            lad = [(x['name'], x['deletion'], x['relaxed'], apo[x['name']]['apo'])
                   for x in ok if x['name'] in apo]
            if lad:
                print(f"\nTHREE-RUNG LADDER (bound leg identical in all three)")
                print(f"{'design':<42}{'deletion':>10}{'relaxed':>9}{'apo':>9}{'span':>8}")
                for n, dd, rr, aa in sorted(lad, key=lambda t: -t[1]):
                    sp = max(dd, rr, aa) / min(dd, rr, aa) if min(dd, rr, aa) else None
                    print(f"{n[:41]:<42}{dd:>10.3f}{rr:>9.3f}{aa:>9.3f}"
                          + (f"{sp:>8.2f}" if sp else f"{'-':>8}"))
                spans = [max(t[1:]) / min(t[1:]) for t in lad if min(t[1:])]
                print(f"\nspan across the three free legs: median {st.median(spans):.3f}x, "
                      f"max {max(spans):.3f}x over {len(spans)} designs")
                summary['ladder_n'] = len(spans)
                summary['ladder_span_median'] = round(st.median(spans), 4)
                summary['ladder_span_max'] = round(max(spans), 4)

    json.dump(dict(rows=rows, summary=summary), open(out_path, 'w'), indent=1)
    print(f"\nwrote {out_path}")


if __name__ == '__main__':
    main()
