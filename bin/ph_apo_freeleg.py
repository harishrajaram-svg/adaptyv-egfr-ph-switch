#!/usr/bin/env python3
"""Apo free leg vs partner-deletion free leg: the other half of the reviewer's sensitivity ask.

What the reviewer told us on 2026-10-05, in substance:

  - Taking the partner out of the file is by itself enough to drop its atoms from the
    burial sum, which is why a histidine in the interface can read as solvent-exposed
    even when the binder has not shifted at all.
  - WHAT THE METHOD CANNOT ESCAPE IS THAT IT KEEPS THE BOUND GEOMETRY. Everything the
    genuinely unbound chain does -- side chains settling into new rotamers, water
    working its way in, the redistribution of conformational populations, protonation
    states coupling to one another -- is simply absent from the model. That shortfall is
    symmetric: it hits the target leg and the binder leg alike.
  - FOLDING THE CHAIN BY ITSELF AND CALLING THAT THE APO TRUTH DOES NOT REPAIR THIS. It
    is a second approximation, with no guarantee that it is the right one.
  - The prescribed handling: set the deletion number beside apo and relaxed numbers
    produced under the same preparation, and read the spread as SENSITIVITY -- never as
    a measured confidence interval.

So this is a COMPARISON OF TWO APPROXIMATIONS, not a correction of one by the other.

  deletion free leg  the bound complex with the partner's atoms removed in place. The
                     conformation is the bound one; only burial changes.
  apo free leg       the same chain folded ALONE by the same predictor, same 5 seeds,
                     same pipeline (runs/esmfold2/w5_apo). The conformation is the
                     predictor's unbound guess.

Neither is the free protein. They differ in which error they make: deletion holds the
backbone and side chains in a conformation the free protein does not adopt; the apo fold
gives a plausible unbound conformation but one that has no particular relationship to the
bound pose, so pairing its pKa with the complex's bound pKa mixes two structures.

The bound leg is identical in both, so every difference below is attributable to the free
leg alone.

    ph_apo_freeleg.py [--out analysis/01-egfr/ph_apo_freeleg.json]
"""
import csv, glob, json, os, re, statistics as st, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

APO_ROOT = 'runs/esmfold2/w5_apo'
SENS = 'analysis/01-egfr/ph_sensitivity.json'
CSV = 'submissions/01-egfr.csv'
OUT = 'analysis/01-egfr/ph_apo_freeleg.json'
PH_LO, PH_HI = 6.5, 7.4


def link(free, bound):
    def K(ph):
        return (1 + 10 ** (bound - ph)) / (1 + 10 ** (free - ph))
    return K(PH_LO) / K(PH_HI)


def propka_his(pdb_dir, tag, chain='A'):
    """pKa of every HIS on `chain` of an already-written PDB. {resnum: pKa}."""
    import subprocess
    subprocess.run([sys.executable, '-m', 'propka', tag + '.pdb'],
                   capture_output=True, cwd=pdb_dir)
    f = os.path.join(pdb_dir, tag + '.pka')
    out = {}
    if os.path.exists(f):
        for ln in open(f):
            q = ln.split()
            if len(q) > 3 and q[0] == 'HIS' and q[2] == chain:
                try:
                    out[int(q[1])] = float(q[3])
                except ValueError:
                    pass
    return out


def apo_pkas(tag):
    """median HIS pKa per residue over the apo monomer's seeds. {resnum: pKa}."""
    import gemmi, tempfile
    cifs = sorted(glob.glob(os.path.join(APO_ROOT, '**', tag, '*.cif'), recursive=True))
    if not cifs:
        cifs = sorted(glob.glob(os.path.join(APO_ROOT, '**', tag + '*', '*.cif'), recursive=True))
    if not cifs:
        return None, 0
    per = {}
    with tempfile.TemporaryDirectory() as td:
        for i, c in enumerate(cifs):
            st_ = gemmi.read_structure(c)
            st_.setup_entities(); st_.remove_ligands_and_waters()
            # the monomer must be chain A for the propka parse
            st_[0][0].name = 'A'
            st_.write_pdb(os.path.join(td, f'apo{i}.pdb'))
            for num, v in propka_his(td, f'apo{i}').items():
                per.setdefault(num, []).append(v)
    return {k: st.median(v) for k, v in per.items()}, len(cifs)


def main():
    out_path = OUT
    if '--out' in sys.argv:
        out_path = sys.argv[sys.argv.index('--out') + 1]
    sens = json.load(open(SENS))
    ship = {r['sequence'].strip().upper(): r['name'] for r in csv.DictReader(open(CSV))}
    by_name = {v['name']: v for v in sens.values()
               if v.get('seq', '').strip().upper() in ship}
    if not by_name:
        sys.exit("no submitted designs in ph_sensitivity.json")

    rows = []
    print(f"{'design':<44}{'deletion':>10}{'apo':>9}{'fold':>7}{'nHis':>6}{'seeds':>6}")
    for name, v in by_name.items():
        tag = 'apo_' + re.sub(r'[^A-Za-z0-9]+', '_', name)[:48]
        apo, nseed = apo_pkas(tag)
        if apo is None:
            rows.append(dict(name=name, error=f'no apo structure for {tag}'))
            print(f"{name[:43]:<44}  no apo structure yet ({tag})")
            continue
        # per pose: recompute the BINDER histidine product with the apo free leg.
        # The target leg is left on the deletion estimate -- the target apo is a single
        # shared structure and is reported separately below.
        del_meds, apo_meds, nhis = [], [], 0
        for pose, sites in (v.get('sites_all_poses') or {}).items():
            pd_, pa_ = 1.0, 1.0
            for s in sites.values():
                if s['resname'] != 'HIS':
                    continue
                pd_ *= link(s['free'], s['bound'])
                if s['partner'] == 'binder' and s['resnum'] in apo:
                    pa_ *= link(apo[s['resnum']], s['bound'])
                else:
                    pa_ *= link(s['free'], s['bound'])   # target sites: unchanged here
                nhis += 1
            del_meds.append(pd_); apo_meds.append(pa_)
        if not del_meds:
            rows.append(dict(name=name, error='no histidine sites')); continue
        d, a = st.median(del_meds), st.median(apo_meds)
        rows.append(dict(name=name, deletion=round(d, 4), apo=round(a, 4),
                         fold=round(a / d, 4) if d else None,
                         n_binder_his_remapped=sum(1 for s in
                             next(iter(v['sites_all_poses'].values())).values()
                             if s['resname'] == 'HIS' and s['partner'] == 'binder'
                             and s['resnum'] in apo),
                         apo_seeds=nseed))
        print(f"{name[:43]:<44}{d:>10.3f}{a:>9.3f}{a/d:>7.2f}"
              f"{rows[-1]['n_binder_his_remapped']:>6}{nseed:>6}")

    ok = [r for r in rows if 'fold' in r and r['fold']]
    if ok:
        folds = [r['fold'] for r in ok]
        moved = [r for r in ok if abs(r['fold'] - 1.0) >= 0.10]
        print(f"\n{len(ok)} designs compared. apo/deletion fold: median "
              f"{st.median(folds):.3f}, range {min(folds):.3f}-{max(folds):.3f}")
        print(f"{len(moved)} design(s) move by >=10%")
        # does the ORDER change?
        do = [r['name'] for r in sorted(ok, key=lambda r: -r['deletion'])]
        ao = [r['name'] for r in sorted(ok, key=lambda r: -r['apo'])]
        import itertools
        c = dc = 0
        for x, y in itertools.combinations(ok, 2):
            s1 = (do.index(x['name']) - do.index(y['name']))
            s2 = (ao.index(x['name']) - ao.index(y['name']))
            if s1 * s2 > 0: c += 1
            else: dc += 1
        print(f"Kendall tau between the two orderings: {(c-dc)/(c+dc):+.3f}")
    json.dump(dict(rows=rows), open(out_path, 'w'), indent=1)
    print(f"wrote {out_path}")


if __name__ == '__main__':
    main()
