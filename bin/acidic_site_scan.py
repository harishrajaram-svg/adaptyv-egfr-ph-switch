#!/usr/bin/env python3
"""Rank target acidic clusters as Mechanism A sites -- with glycan clearance measured
against SEQUONS IN THE SEQUENCE, not glycans in the coordinates.

WHY THIS FILE EXISTS (2026-10-03)
---------------------------------
The 2026-10-01 scan called all six candidate clusters "conserved / exposed / glycan-free",
and it left no script behind, so its logic could not be audited. It was wrong: it could
only see the 4 N-glycans actually MODELLED in 6ARU out of 11 N-X-S/T sequons in the
sequence. Cluster 5 looked pristine at 31.2 A from the nearest modelled sugar while
sitting 8.8 A from an unmodelled sequon -- underneath a glycan nobody had resolved.
An N-glycan reaches 10-20 A+ from its attachment Asn. Crystal structures systematically
under-model carbohydrate, so a coordinate-based clearance check is biased clean.

THE RULE THIS ENFORCES: clearance is distance to the SEQUON ASN, from the sequence.

Numbering, verified residue-by-residue against egfr_ecd_6aru_renum.pdb:
    renum = canonical - 27 = mature - 3
    H433 canonical = renum 406;  E543/E545 = renum 516/518;
    E319/E320/D321 = renum 292/293/294;  tether Y246/D563 mature = renum 243/560.

Usage:
    acidic_site_scan.py --pdb targets/egfr/egfr_ecd_6aru_renum.pdb \
        --human targets/egfr/egfr_ecd_human.faa --mouse targets/egfr/egfr_ecd_mouse.faa
    acidic_site_scan.py --self-test
"""
import argparse, math, re, sys
from collections import defaultdict
from pathlib import Path

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
CARBOX = {('ASP','OD1'),('ASP','OD2'),('GLU','OE1'),('GLU','OE2')}
CLUSTER_D   = 11.0   # A, carboxylate-O to carboxylate-O; single linkage. Set a priori to
                     # match the span one binder face can engage (cf. E543/E545 at 10.9 A).
BURIAL_D    = 10.0   # A, neighbour-residue count radius
GLYCAN_REACH= 20.0   # A, an N-glycan's outer reach from its Asn (10-20 A+; use the ceiling)
TETHER      = (243, 560)   # renum of Y246 / D563 mature
OFF_CANON   = 27     # canonical = renum + OFF_CANON
OFF_MATURE  = 3      # mature    = renum + OFF_MATURE

DOMAINS = [("I",1,165),("II",166,310),("III",311,480),("IV",481,620)]  # mature numbering


def sequons(seq, first_idx=1, atypical=True):
    """N-glycosylation sequons. Returns [(1-based index of the Asn, motif)].

    Canonical N-X-S/T with X != P, PLUS the atypical N-X-C motif when `atypical`.
    PK flagged the omission on 2026-10-03: EGFR carries a documented atypical
    N-X-C site at mature N32 / canonical N56 (motif NNC) that an N-X-S/T-only
    scan cannot see. Scanning for it finds FOUR N-X-C sequons on the ECD
    (canonical 56, 234, 468, 497), lifting the sequon count from 11 to 15.
    Overlapping matches are both real, so scan positions rather than using
    re.finditer, which consumes the first match's characters.

    A sequon is POTENTIAL occupancy, not established occupancy -- it is the
    conservative direction for a clearance check, which is the point.
    """
    third = 'STC' if atypical else 'ST'
    out = []
    for i in range(len(seq) - 2):
        n, x, t = seq[i], seq[i+1], seq[i+2]
        if n == 'N' and x != 'P' and t in third:
            out.append((first_idx + i, n + x + t))
    return out


def load(pdb):
    atoms = defaultdict(list)   # renum -> [(atomname, resname, (x,y,z))]
    for l in Path(pdb).read_text().splitlines():
        if l.startswith('ATOM') and l[21] == 'A' and l[76:78].strip() != 'H':
            atoms[int(l[22:26])].append(
                (l[12:16].strip(), l[17:20].strip(),
                 (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return atoms


def domain_of(renum):
    m = renum + OFF_MATURE
    for name, lo, hi in DOMAINS:
        if lo <= m <= hi:
            return name
    return "?"


def site_domain(group):
    """Majority domain. Using the lowest-numbered member mislabels any site that
    straddles a domain boundary -- E543/E545 (domain IV) came out as 'III'."""
    c = defaultdict(int)
    for r in group:
        c[domain_of(r)] += 1
    return max(sorted(c), key=lambda k: c[k])


def sites(acidic, atoms):
    """Candidate sites = the NEIGHBOURHOOD of each acidic residue, not a transitive closure.

    Single-linkage was the first thing tried and it is wrong here: with 65 acidic residues
    on one surface it chains them into 12-member "clusters" spanning two domains, and the
    glycan clearance of such a group is the clearance of its worst member, not of the two
    or three residues a binder would actually engage. A binder face sees a local
    neighbourhood, so that is what gets scored. Sites are deduplicated by residue set.
    """
    ox = {r: [a[2] for a in atoms[r] if (a[1], a[0]) in CARBOX] for r in acidic}
    seen, out = set(), []
    for r in acidic:
        grp = tuple(sorted(s for s in acidic
                           if any(math.dist(p, q) <= CLUSTER_D
                                  for p in ox[r] for q in ox[s])))
        if len(grp) >= 2 and grp not in seen:
            seen.add(grp); out.append(list(grp))
    return out, ox


def scan(pdb, human, mouse):
    atoms = load(pdb)
    resname = {r: atoms[r][0][1] for r in atoms}
    hs = ''.join(l.strip() for l in Path(human).read_text().splitlines() if not l.startswith('>'))
    ms = ''.join(l.strip() for l in Path(mouse).read_text().splitlines() if not l.startswith('>'))
    if len(hs) != len(ms):
        sys.exit(f"human ({len(hs)}) and mouse ({len(ms)}) ECD lengths differ -- needs alignment")

    # sequons from the SEQUENCE (mature 1-based), mapped to renum
    seqs = [(i - OFF_MATURE, m) for i, m in sequons(hs, 1)]
    modelled = [(r, m) for r, m in seqs if r in atoms]
    print(f"N-glycosylation sequons in the human ECD sequence (N-X-S/T and the atypical "
          f"N-X-C): {len(seqs)}  (Asn present in the structure: {len(modelled)})")
    print("  " + ", ".join(f"N{r+OFF_CANON}" for r, _ in seqs))

    acidic = [r for r in sorted(atoms) if resname[r] in ('ASP', 'GLU')]
    groups, ox = sites(acidic, atoms)
    print(f"\n{len(acidic)} acidic residues -> {len(groups)} clusters of >=2 "
          f"(single linkage at {CLUSTER_D} A)\n")

    rows = []
    for g in groups:
        pts = [p for r in g for p in ox[r]]
        # burial: residues with any heavy atom within BURIAL_D of any cluster carboxylate O
        nb = sum(1 for r2 in atoms if r2 not in g and
                 any(math.dist(a[2], p) <= BURIAL_D for a in atoms[r2] for p in pts))
        # glycan clearance: to the sequon ASN's own atoms (sequence-derived)
        gl = min(((min(math.dist(a[2], p) for a in atoms[r] for p in pts), r)
                  for r, _ in modelled if r not in g), default=(float('inf'), None))
        # conservation human vs mouse at the cluster positions
        cons = sum(1 for r in g if hs[r + OFF_MATURE - 1] == ms[r + OFF_MATURE - 1])
        teth = min(min(math.dist(a[2], p) for a in atoms[t] for p in pts)
                   for t in TETHER if t in atoms)
        rows.append(dict(renum=g, canon=[r + OFF_CANON for r in g], n=len(g), burial=nb,
                         glyc=gl[0], glyc_res=gl[1], cons=cons, teth=teth,
                         dom=site_domain(g)))

    # H433 as a calibration row: the site mechanism B used, same burial metric
    h = 406
    if h in atoms:
        hb = sum(1 for r2 in atoms if r2 != h and
                 any(math.dist(a[2], b[2]) <= BURIAL_D for a in atoms[r2] for b in atoms[h]))
        print(f"calibration -- H433 (renum 406, mechanism B's site): burial {hb}\n")

    rows.sort(key=lambda r: (-min(r['glyc'], GLYCAN_REACH * 2), -r['n'], r['burial']))
    print(f"{'canonical acidic residues':<34}{'dom':>4}{'n':>3}{'cons':>6}"
          f"{'burial':>8}{'glycan':>8}{'sequon':>9}{'tether':>8}  verdict")
    for r in rows:
        v = ("CLEAR" if r['glyc'] > GLYCAN_REACH else
             "UNDER-GLYCAN" if r['glyc'] < 10 else "MARGINAL")
        if r['cons'] < r['n']:
            v += ",NOT-CONSERVED"
        if r['teth'] < 15:
            v += ",AT-TETHER"
        print(f"{','.join(str(c) for c in r['canon'])[:33]:<34}{r['dom']:>4}{r['n']:>3}"
              f"{str(r['cons'])+'/'+str(r['n']):>6}{r['burial']:>8}{r['glyc']:>8.1f}"
              f"{('N'+str(r['glyc_res']+OFF_CANON)) if r['glyc_res'] else '-':>9}"
              f"{r['teth']:>8.1f}  {v}")
    return rows


def neutral_sites(pdb, human, mouse, top=12):
    """The inverse scan: exposed, conserved, glycan-clear patches with NO carboxylate near.

    Why this exists: the carboxylate-proximity requirement (a target Asp/Glu within 4.0 A of
    the imidazole) was called NECESSARY off 6/6 switches having one and 0/96 non-switches
    switching. That is correlational, and every Mechanism A design we ever made was aimed AT
    an acidic site -- so the comparison has no negative arm. Pinning a histidine against a
    carboxylate-free patch is the direct test. If switches appear here too, the requirement
    is an artifact of where we chose to design.
    """
    atoms = load(pdb)
    resname = {r: atoms[r][0][1] for r in atoms}
    hs = ''.join(l.strip() for l in Path(human).read_text().splitlines() if not l.startswith('>'))
    ms = ''.join(l.strip() for l in Path(mouse).read_text().splitlines() if not l.startswith('>'))
    acid_ox = [a[2] for r in atoms if resname[r] in ('ASP', 'GLU')
               for a in atoms[r] if (a[1], a[0]) in CARBOX]
    seqs = [(i - OFF_MATURE, m) for i, m in sequons(hs, 1)]
    modelled = [r for r, _ in seqs if r in atoms]

    cand = []
    for r in sorted(atoms):
        if resname[r] in ('ASP', 'GLU', 'GLY', 'PRO'):
            continue
        pts = [a[2] for a in atoms[r]]
        near_acid = min((math.dist(p, q) for p in pts for q in acid_ox), default=99.0)
        if near_acid <= CLUSTER_D:
            continue                                   # a carboxylate is in reach -- not neutral
        burial = sum(1 for r2 in atoms if r2 != r and
                     any(math.dist(a[2], p) <= BURIAL_D for a in atoms[r2] for p in pts))
        glyc = min((min(math.dist(a[2], p) for a in atoms[n] for p in pts)
                    for n in modelled if n != r), default=99.0)
        cons = hs[r + OFF_MATURE - 1] == ms[r + OFF_MATURE - 1]
        if glyc <= GLYCAN_REACH or not cons or burial > 45:
            continue
        cand.append((near_acid, r, resname[r], burial, glyc))

    cand.sort(key=lambda t: (-t[0], t[3]))
    print(f"carboxylate-free, conserved, exposed, glycan-clear residues "
          f"(no Asp/Glu within {CLUSTER_D} A, burial <= 45, sequon > {GLYCAN_REACH} A)\n")
    print(f"{'canonical':>10}{'aa':>5}{'dom':>5}{'burial':>8}{'nearest acid':>14}{'sequon':>9}")
    for na, r, rn, b, g in cand[:top]:
        print(f"{r+OFF_CANON:>10}{rn:>5}{domain_of(r):>5}{b:>8}{na:>14.1f}{g:>9.1f}")
    return cand


def self_test():
    assert sequons("AANATA", 1) == [(3, 'NAT')], sequons("AANATA", 1)
    assert sequons("AANPSA", 1) == [], "X=P must be rejected"
    assert sequons("AANASA", 1) == [(3, 'NAS')], sequons("AANASA", 1)
    # N-A-C used to be asserted as a non-sequon here. That assertion encoded the
    # very omission PK flagged: N-X-C is a real, if atypical, glycosylation motif.
    assert sequons("AANACA", 1) == [(3, 'NAC')], sequons("AANACA", 1)
    assert sequons("AANACA", 1, atypical=False) == [], "S/T-only mode must still reject it"
    assert sequons("AANAVA", 1) == [], "+2 must be S, T or C"
    # overlapping sequons are both real and both must be reported; a single
    # re.finditer pass consumes the first match's characters and misses the second.
    assert sequons("ANNSTX", 1) == [(2, 'NNS'), (3, 'NST')], sequons("ANNSTX", 1)
    assert sequons("AANAT", 5) == [(7, 'NAT')], "first_idx offset"
    # the atypical N-X-C motif PK flagged -- EGFR mature N32 is NNC
    assert sequons("AANNCA") == [(3, 'NNC')], sequons("AANNCA")
    assert sequons("AANNCA", atypical=False) == [], "N-X-C must be opt-out-able"
    assert sequons("AANPCA") == [], "X=P must be rejected for N-X-C too"
    assert domain_of(292) == "II" and domain_of(516) == "IV" and domain_of(406) == "III"
    print("self-test OK: sequon rule (N-X-S/T, X!=P, overlapping) + domain mapping")


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        self_test(); sys.exit()
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True); ap.add_argument("--human", required=True)
    ap.add_argument("--mouse", required=True)
    ap.add_argument("--neutral", action="store_true",
                    help="inverse scan: carboxylate-FREE control patches")
    a = ap.parse_args()
    (neutral_sites if a.neutral else scan)(a.pdb, a.human, a.mouse)
