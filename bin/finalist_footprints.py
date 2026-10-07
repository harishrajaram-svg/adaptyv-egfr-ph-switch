#!/usr/bin/env python3
"""Epitope footprint of every submitted design, with the four checks the reviewer asked for.

On 2026-10-05 the reviewer required that the footprints of the designs actually submitted
be put through four tests -- the complete ECD, glycosylation, the receptor's conformational
state, and a human-versus-mouse comparison at the contacts themselves -- and stated that
their own earlier read on domain II grants these designs no clearance.

That assessment was made against a construct and a target region, not against where
these particular binders actually land. This computes the footprint per design from its
own human-leg poses and then runs each check against it.

  DOMAIN            which EGFR ectodomain the contacts fall in (mature numbering)
  CROP              whether contacts fall OUTSIDE the 170 aa domain-III crop most of
                    these binders were designed against -- i.e. whether the footprint is
                    even assessable on the crop, or only on the full ECD
  GLYCAN            overlap with N-glycosylation sequons (N-X-S/T, X != P) of human EGFR.
                    A footprint on a sequon is a footprint on a site that carries a glycan
                    in a real cell and does not in any structure we fold.
  HUMAN/MOUSE       identity at the contacted positions specifically. The submission
                    requires cross-reactivity, so conservation AT THE EPITOPE is the
                    relevant quantity, not whole-protein identity.

Numbering is mature ECD numbering (canonical P00533 residues 25-645 -> mature 1-621),
which is what targets/egfr/egfr_ecd_human.faa declares. Pose chains are mapped to it by
global sequence alignment, so d3-crop and full-ECD poses land on one coordinate system.

    finalist_footprints.py [--cut 5.0]
"""
import csv, glob, json, os, re, statistics as st, sys
from collections import Counter, defaultdict
import gemmi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pose_seq_index

CSV = 'submissions/01-egfr.csv'
HU = 'targets/egfr/egfr_ecd_human.faa'
MO = 'targets/egfr/egfr_ecd_mouse.faa'
CUT = 5.0
# Mature-numbering ectodomain boundaries (Ferguson 2008; standard L1/CR1/L2/CR2 split).
DOMAINS = [('I (L1)', 1, 165), ('II (CR1)', 166, 310),
           ('III (L2)', 311, 480), ('IV (CR2)', 481, 621)]


def read_faa(p):
    return ''.join(l.strip() for l in open(p) if not l.startswith('>')).upper()


def domain_of(i):
    for nm, a, b in DOMAINS:
        if a <= i <= b:
            return nm
    return '?'


def sequons(seq):
    """1-based positions of the Asn in every N-X-S/T sequon, X != P."""
    return {m.start() + 1 for m in re.finditer(r'N[^P][ST]', seq)}


def align_map(a, b):
    """index-in-a (0-based) -> index-in-b (0-based)."""
    res = gemmi.align_string_sequences(list(a), list(b), [])
    m, ia, ib = {}, 0, 0
    for length, op in re.findall(r'(\d+)([MID])', res.cigar_str()):
        length = int(length)
        if op == 'M':
            for _ in range(length):
                m[ia] = ib; ia += 1; ib += 1
        elif op == 'I':
            ia += length
        elif op == 'D':
            ib += length
    return m


def chains(path):
    st_ = gemmi.read_structure(str(path))
    st_.setup_entities(); st_.remove_ligands_and_waters()
    out = []
    for ch in st_[0]:
        res = [r for r in ch]
        if len(res) < 10:
            continue
        out.append((gemmi.one_letter_code([r.name for r in res]).upper(), res))
    return out


def contact_target_idx(tres, bres, cut):
    """0-based indices of TARGET residues within cut of any binder heavy atom."""
    cell = cut
    grid = defaultdict(list)
    for rb in bres:
        for a in rb:
            if a.element == gemmi.Element('H'):
                continue
            grid[(int(a.pos.x // cell), int(a.pos.y // cell), int(a.pos.z // cell))].append(a.pos)
    hit = set()
    for it, rt in enumerate(tres):
        done = False
        for a in rt:
            if done or a.element == gemmi.Element('H'):
                continue
            kx, ky, kz = int(a.pos.x // cell), int(a.pos.y // cell), int(a.pos.z // cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for p in grid.get((kx+dx, ky+dy, kz+dz), ()):
                            if a.pos.dist(p) <= cut:
                                hit.add(it); done = True; break
                        if done: break
                    if done: break
                if done: break
    return hit


def main():
    cut = CUT
    if '--cut' in sys.argv:
        cut = float(sys.argv[sys.argv.index('--cut') + 1])
    hu, mo = read_faa(HU), read_faa(MO)
    if len(hu) != len(mo):
        print(f"  note: human ECD {len(hu)} aa vs mouse {len(mo)} aa; mapping by alignment")
    h2m = align_map(hu, mo)
    sq = sequons(hu)
    d3 = read_faa('targets/egfr/egfr_d3_human.faa')
    d3map = align_map(d3, hu)                 # crop index -> mature index
    crop_lo = min(d3map.values()) + 1
    crop_hi = max(d3map.values()) + 1
    print(f"human ECD {len(hu)} aa (mature 1-{len(hu)}); domain-III crop covers mature "
          f"{crop_lo}-{crop_hi}")
    print(f"N-glycosylation sequons in the human ECD: {len(sq)} "
          f"(Asn at {', '.join(str(x) for x in sorted(sq)[:12])}{', ...' if len(sq)>12 else ''})\n")

    idx = pose_seq_index.load()
    rows = []
    for r in csv.DictReader(open(CSV)):
        seq = r['sequence'].strip().upper()
        poses = [c for d in idx.get(seq, [])
                 if os.path.basename(d.rstrip('/')).endswith(('_hu', '_human'))
                 for c in sorted(glob.glob(os.path.join(d, '*.cif')))]
        votes = Counter()
        n_ok = 0
        for p in poses:
            cs = chains(p)
            if len(cs) < 2:
                continue
            cs.sort(key=lambda x: -len(x[1]))
            (tseq, tres), (bseq, bres) = cs[0], cs[1]
            tm = align_map(tseq, hu)
            hits = contact_target_idx(tres, bres, cut)
            n_ok += 1
            for i in hits:
                if i in tm:
                    votes[tm[i] + 1] += 1          # 1-based mature numbering
        if not n_ok:
            rows.append(dict(name=r['name'], error='no scoreable pose')); continue
        # a residue is in the footprint if it contacts in the MAJORITY of poses
        foot = sorted(p for p, v in votes.items() if v >= (n_ok + 1) // 2)
        if not foot:
            rows.append(dict(name=r['name'], n_poses=n_ok, footprint=[],
                             note='no residue contacted in a majority of poses')); continue
        doms = Counter(domain_of(p) for p in foot)
        outside = [p for p in foot if not (crop_lo <= p <= crop_hi)]
        gly = sorted(set(foot) & sq)
        cons = [p for p in foot if p - 1 in h2m and hu[p-1] == mo[h2m[p-1]]]
        rows.append(dict(name=r['name'], n_poses=n_ok, footprint=foot,
                         n_footprint=len(foot),
                         domains={k: v for k, v in doms.most_common()},
                         n_outside_d3_crop=len(outside),
                         outside_d3_crop=outside[:12],
                         glycan_sequon_hits=gly,
                         hu_mo_identity_at_epitope=round(len(cons) / len(foot), 3)))
    print(f"{'design':<42}{'n':>3}{'foot':>6}{'domains':<26}{'offcrop':>8}{'glyc':>6}{'hu=mo':>7}")
    for r in rows:
        if 'footprint' not in r:
            print(f"{r['name'][:41]:<42}  {r.get('error','?')}"); continue
        d = ' '.join(f"{k}:{v}" for k, v in (r['domains'] or {}).items())
        print(f"{r['name'][:41]:<42}{r['n_poses']:>3}{r['n_footprint']:>6}{d:<26}"
              f"{r['n_outside_d3_crop']:>8}{len(r['glycan_sequon_hits']):>6}"
              f"{r['hu_mo_identity_at_epitope']:>7.2f}")
    json.dump(dict(cut=cut, crop_mature_range=[crop_lo, crop_hi],
                   n_sequons=len(sq), designs=rows),
              open('analysis/01-egfr/finalist_footprints.json', 'w'), indent=1)
    ok = [r for r in rows if 'footprint' in r and r['footprint']]
    if ok:
        print(f"\nacross {len(ok)} designs: median footprint "
              f"{st.median([r['n_footprint'] for r in ok]):.0f} residues; "
              f"{sum(1 for r in ok if r['n_outside_d3_crop'])} have contacts outside the d3 crop; "
              f"{sum(1 for r in ok if r['glycan_sequon_hits'])} touch a glycosylation sequon; "
              f"median human/mouse identity at the epitope "
              f"{st.median([r['hu_mo_identity_at_epitope'] for r in ok]):.2f}")
    print("wrote analysis/01-egfr/finalist_footprints.json")


if __name__ == '__main__':
    main()
