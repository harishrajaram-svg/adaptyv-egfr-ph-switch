#!/usr/bin/env python3
"""Index every scored pose directory by its BINDER SEQUENCE, read from the structure.

WHY. Every sequence->pose join in this project went through
`analysis/01-egfr/score_*/*.faa`: a pose directory was found only if a .faa file of a
matching basename existed in a score_* directory. That index is incomplete by
construction -- a new Modal run has no .faa until someone writes one -- and it is held
together by 8 untracked symlinks (score_w1_*, score_w2_*) that exist purely so a glob
would match. PK flagged the symlinks as missing from the public tree; the deeper problem
is that a missing index entry is silent: the pose is simply not pooled, and the design's
n is quietly too small.

That bit tonight. runs/esmfold2/w3_triad holds 15 poses of a binder sequence BYTE-
IDENTICAL to the submitted rimA01_r15_L133E, and no .faa existed for them, so the
submission kept reporting n=5 for a sequence with 20 poses on disk.

This index reads the sequence out of the structure itself, so any run is picked up with
no bookkeeping. The binder is the SHORTER chain (the EGFR constructs are 170/621 aa; no
submitted binder exceeds 150). Cached, because it is 4.5k structure reads.

    pose_seq_index.py                 # build/refresh the cache
    pose_seq_index.py --stats         # what it found
"""
import glob, json, os, sys

CACHE = 'analysis/01-egfr/pose_seq_index.json'

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}


def chains_of(cif):
    """{chain: one-letter sequence} from the atom_site loop, in file order.

    A light text parse rather than gemmi: this runs over thousands of structures, and the
    only fields needed are chain id, residue number and residue name."""
    hdr, out, seen, inloop, done = [], {}, set(), False, False
    try:
        fh = open(cif)
    except OSError:
        return {}
    with fh:
        for ln in fh:
            if ln.startswith('_atom_site.'):
                hdr.append(ln.strip()); inloop = True; continue
            if not inloop:
                continue
            if done:
                break
            if ln.startswith('#') or (ln.strip() and ln.startswith('_')):
                done = True; continue
            if not ln.strip():
                continue
            f = ln.split()
            if len(f) < len(hdr):
                continue
            try:
                ci = hdr.index('_atom_site.auth_asym_id')
            except ValueError:
                ci = hdr.index('_atom_site.label_asym_id')
            try:
                ri = hdr.index('_atom_site.auth_seq_id')
            except ValueError:
                ri = hdr.index('_atom_site.label_seq_id')
            try:
                ni = hdr.index('_atom_site.label_comp_id')
            except ValueError:
                ni = hdr.index('_atom_site.auth_comp_id')
            key = (f[ci], f[ri])
            if key in seen:
                continue
            seen.add(key)
            out.setdefault(f[ci], []).append(AA3.get(f[ni].upper(), 'X'))
    return {k: ''.join(v) for k, v in out.items()}


def build():
    dirs = sorted({os.path.dirname(c) for c in glob.glob('runs/**/*.cif', recursive=True)})
    idx, skipped = {}, 0
    for i, d in enumerate(dirs):
        cifs = sorted(glob.glob(os.path.join(d, '*.cif')))
        ch = chains_of(cifs[0]) if cifs else {}
        if len(ch) < 2:
            skipped += 1; continue
        binder = min(ch.values(), key=len)
        idx.setdefault(binder.upper(), []).append(d)
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{len(dirs)} dirs", flush=True)
    json.dump({'by_sequence': {k: sorted(v) for k, v in idx.items()}},
              open(CACHE, 'w'))
    print(f"wrote {CACHE}: {len(idx)} distinct binder sequences over "
          f"{len(dirs)-skipped} pose dirs ({skipped} skipped: <2 chains)")
    return idx


def load():
    if not os.path.exists(CACHE):
        return build()
    return json.load(open(CACHE))['by_sequence']


if __name__ == '__main__':
    idx = build()
    if '--stats' in sys.argv:
        import csv
        sub = list(csv.DictReader(open('submissions/01-egfr.csv')))
        print(f"\n{'design':<42}{'pose dirs':>10}{'cifs':>7}")
        for r in sub:
            ds = idx.get(r['sequence'].strip().upper(), [])
            n = sum(len(glob.glob(os.path.join(d, '*.cif'))) for d in ds)
            print(f"{r['name'][:41]:<42}{len(ds):>10}{n:>7}")
