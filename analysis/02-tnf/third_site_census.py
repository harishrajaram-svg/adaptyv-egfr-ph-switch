#!/usr/bin/env python3
"""Problem 2: census of target-side cations that could anchor a BINDER histidine.

Why this exists
---------------
The switch mechanism is binder-side: a histidine placed against a target cation has its pKa
suppressed on binding (the His+ / Arg+ or Lys+ repulsion destabilises the protonated form), so
the complex tolerates only the neutral His. Acidify to 6.0, the His protonates, binding is lost.

The working design target named R107/R108/K166. R107 is Gln in mouse, so it cannot carry a
load-bearing contact (objective 2 is assayed at pH 7.4, where that contact pays for affinity).
That leaves two anchors -> 51.6x, not the 370x that three independent sites would give.

This script answers: is there a third anchor, and is it independent of the first two?

Run it:  python3 analysis/02-tnf/third_site_census.py
Depends: per_partner.json (from per_partner.py), structures/1TNF.pdb (fetched if absent).

Self-tests abort on failure -- every expensive failure in problem 1 was a clean, plausible,
wrong number (playbook Part III).
  T1  the three named residues are classified as expected: R107 divergent, R108/K166 conserved
  T2  the consensus core recomputes to the same 21 residues the day-1 analysis published
  T3  every candidate's residue type in the PDB matches the UniProt sequence at that position
"""
import json, os, sys, urllib.request, itertools

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from species_and_histidines import HUMAN, UNIPROT_START, species_differences

import gemmi

CATIONIC = {"ARG": "CZ", "LYS": "NZ", "HIS": "CE1"}   # charged-group proxy atom
COUPLING_A = 12.0   # below this, two titratable sites are not safely independent (playbook s3)
BURIAL_R = 10.0     # heavy atoms within this radius = burial proxy

# The ten receptor copies the day-1 analysis used, as (entry, partner-chain) pairs.
#
# \U0001F534 FOUR OF THE TEN ARE A MUTEIN, found 2026-10-06 by the widened identity guard
# (analysis/02-tnf/target_identity.py). 3ALQ's own PDB entity record carries
# pdbx_mutation = "K11M, K65S, K90P, K98R, K112N, K128P" = K87M, K141S, K166P, K174R, K188N,
# K204P in our numbering. THREE of those are cationic sites this census ranks, and K166 -- the
# primary anchor -- IS PROLINE in 3ALQ.
#
# What that costs, and what it does not. The CONTACT POSITIONS are still real: 166 is in the
# receptor-contact set of 8ZUI and 7KPB as well, both wild type and both clean through the
# guard, so the epitope does not rest on the mutein. But a conservation count over these ten
# copies measures POSITION, not residue IDENTITY, at the 3ALQ rows. Anything phrased as
# "K166 is present in 10 of 10 receptor copies" must be restated as SIX of ten wild-type copies
# plus four copies at the same position in a mutein.
#
# This census's T3 asserts residue identity "only at cationic positions" -- against a single
# 1TNF reference, not per entry, which is why it passed while four of its ten rows carried a
# proline there. Scoped-guard failure, playbook s24, same shape as the 1TNF error itself.
WILD_TYPE_ENTRIES = {"8ZUI", "7KPB"}
MUTEIN_ENTRIES = {"3ALQ": "K87M K141S K166P K174R K188N K204P"}
RECEPTOR_COPIES = [("3ALQ", c) for c in ("T", "R", "V", "U")] + \
                  [("8ZUI", c) for c in ("E", "K", "D", "J")] + \
                  [("7KPB", c) for c in ("E", "F")]
N_WILD_TYPE_COPIES = sum(1 for e, _ in RECEPTOR_COPIES if e in WILD_TYPE_ENTRIES)


def fetch_1tnf(path):
    if not os.path.exists(path):
        url = "https://files.rcsb.org/download/1TNF.pdb"
        sys.stderr.write(f"fetching {url}\n")
        urllib.request.urlretrieve(url, path)
    return path


def load_epitope(pp):
    """-> (consensus core set, union set, {uniprot: n_copies_containing_it})"""
    per_copy = [set(pp[e][c]["total"]) for e, c in RECEPTOR_COPIES]
    counts = {}
    for s in per_copy:
        for u in s:
            counts[u] = counts.get(u, 0) + 1
    return set.intersection(*per_copy), set.union(*per_copy), counts


def protomer_roles(pp):
    """For each receptor copy, which uniprot positions sit on the major vs minor protomer.
    Returns {uniprot: {'major': n, 'minor': n}} counted over the ten copies."""
    role = {}
    for e, c in RECEPTOR_COPIES:
        per = pp[e][c]["per_protomer"]
        ranked = sorted(per.items(), key=lambda kv: -len(kv[1]))
        for rank, (_chain, residues) in enumerate(ranked[:2]):
            tag = "major" if rank == 0 else "minor"
            for u in residues:
                role.setdefault(u, {"major": 0, "minor": 0})[tag] += 1
    return role


def tnf_chain_offset(ch):
    """Solve uniprot = seqid + k by matching to HUMAN; returns (k, identity)."""
    obs = []
    for r in ch:
        info = gemmi.find_tabulated_residue(r.name)
        if info and info.is_amino_acid():
            obs.append((r.seqid.num, info.one_letter_code.upper()))
    if len(obs) < 40:
        return None, 0.0
    best = (None, 0.0)
    for k in range(40, 120):
        m = t = 0
        for num, c in obs:
            u = num + k - UNIPROT_START
            if 0 <= u < len(HUMAN):
                t += 1
                m += (c == HUMAN[u])
        if t >= 40 and m / t > best[1]:
            best = (k, m / t)
    return best


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    pp = json.load(open(os.path.join(here, "per_partner.json")))
    diffs, _, _ = species_differences()
    core, union, counts = load_epitope(pp)
    roles = protomer_roles(pp)
    aa = lambda u: HUMAN[u - UNIPROT_START]

    # ---- self-tests -------------------------------------------------------
    assert 107 in diffs and diffs[107] == "R107Q", f"T1 FAIL: R107 not divergent, got {diffs.get(107)}"
    assert 108 not in diffs, "T1 FAIL: R108 should be conserved"
    assert 166 not in diffs, "T1 FAIL: K166 should be conserved"
    assert len(core) == 21, f"T2 FAIL: consensus core is {len(core)}, day-1 published 21"

    pdb = fetch_1tnf(os.path.join(here, "structures", "1TNF.pdb"))
    st = gemmi.read_structure(pdb)
    st.setup_entities()
    st.remove_ligands_and_waters()
    model = st[0]

    # map (chain, uniprot) -> charged-group position, for cationic residues only
    sites, offsets = {}, {}
    for ch in model:
        k, acc = tnf_chain_offset(ch)
        if not k or acc < 0.90:
            continue
        offsets[ch.name] = k
        for r in ch:
            u = r.seqid.num + k
            if r.name in CATIONIC and 0 <= u - UNIPROT_START < len(HUMAN):
                assert aa(u) == gemmi.find_tabulated_residue(r.name).one_letter_code.upper(), \
                    f"T3 FAIL: {ch.name}:{r.name}{r.seqid.num} -> U{u} expected {aa(u)}"
                at = r.find_atom(CATIONIC[r.name], "*")
                if at:
                    sites[(ch.name, u)] = at.pos
    assert offsets, "T3 FAIL: no TNF chain matched P01375"
    print(f"self-tests T1 T2 T3 passed (core={len(core)}, chains={sorted(offsets)}, "
          f"offset +{sorted(set(offsets.values()))[0]}, cationic sites={len(sites)})\n")

    # burial proxy: heavy atoms within BURIAL_R of the charged group
    allatoms = [a.pos for c in model for r in c for a in r if a.element != gemmi.Element("H")]
    def burial(p):
        return sum(1 for q in allatoms if p.dist(q) <= BURIAL_R)

    # ---- the census -------------------------------------------------------
    cand = sorted({u for (_c, u) in sites} & union)
    print(f"CATIONIC RESIDUES (Arg/Lys/His) IN THE RECEPTOR-CONTACT UNION: {len(cand)}")
    print(f"  union={len(union)} residues, consensus core={len(core)}\n")
    print(f"  {'site':<7}{'core?':<7}{'copies':<8}{'mouse':<10}{'protomer':<12}{'burial':<8}verdict")
    rows = []
    for u in cand:
        conserved = u not in diffs
        incore = u in core
        r = roles.get(u, {"major": 0, "minor": 0})
        side = f"{r['major']}maj/{r['minor']}min"
        b = int(sum(burial(p) for (c, uu), p in sites.items() if uu == u) /
                max(1, sum(1 for (c, uu) in sites if uu == u)))
        if not conserved:
            verdict = f"EXCLUDED ({diffs[u]})"
        elif incore:
            verdict = "CANDIDATE"
        else:
            verdict = "peripheral"
        rows.append((u, incore, counts.get(u, 0), conserved, side, b, verdict))
        print(f"  {aa(u)}{u:<6}{'YES' if incore else '-':<7}{counts.get(u,0):<8}"
              f"{('same' if conserved else diffs[u]):<10}{side:<12}{b:<8}{verdict}")

    # ---- independence from the two surviving anchors ----------------------
    anchors = [108, 166]
    # Any conserved cation is eligible, core or peripheral -- the core has only the two anchors,
    # so a third site necessarily comes from the periphery.
    cands = [u for (u, _ic, _n, cons, _s, _b, _v) in rows if cons and u not in anchors]
    print(f"\nINDEPENDENCE CHECK -- min distance to the anchors R108 / K166, "
          f"across all protomer pairs (coupling risk below {COUPLING_A:.0f} A)")
    print(f"  {'candidate':<12}{'min d(R108)':<14}{'min d(K166)':<14}verdict")
    for u in cands:
        best = {}
        for a in anchors:
            ds = [p.dist(q) for (c1, u1), p in sites.items() if u1 == u
                  for (c2, u2), q in sites.items() if u2 == a]
            best[a] = min(ds) if ds else float("nan")
        ok = all(best[a] >= COUPLING_A for a in anchors)
        print(f"  {aa(u)}{u:<11}{best[108]:<14.1f}{best[166]:<14.1f}"
              f"{'independent' if ok else 'COUPLED -- not a third site'}")

    # ---- anchor geometry: are R108 and K166 even on the same protomer? ----
    print("\nANCHOR GEOMETRY (the two survivors)")
    for a, b in itertools.combinations(anchors, 2):
        ds = sorted((p.dist(q), c1, c2) for (c1, u1), p in sites.items() if u1 == a
                    for (c2, u2), q in sites.items() if u2 == b)
        same = [d for d, c1, c2 in ds if c1 == c2]
        cross = [d for d, c1, c2 in ds if c1 != c2]
        print(f"  {aa(a)}{a} <-> {aa(b)}{b}: closest same-protomer "
              f"{min(same) if same else float('nan'):.1f} A, "
              f"closest cross-protomer {min(cross) if cross else float('nan'):.1f} A")
    r = roles.get(108, {}), roles.get(166, {})
    print(f"  protomer role over the {len(RECEPTOR_COPIES)} receptor copies: "
          f"R108 {r[0]}, K166 {r[1]}")
    print(f"  \u26a0\ufe0f  {N_WILD_TYPE_COPIES} of those {len(RECEPTOR_COPIES)} are WILD TYPE "
          f"({sorted(WILD_TYPE_ENTRIES)}). The other "
          f"{len(RECEPTOR_COPIES) - N_WILD_TYPE_COPIES} are 3ALQ, a six-lysine mutein carrying "
          f"{MUTEIN_ENTRIES['3ALQ']} -- so K166 is PROLINE in them. A count over all ten "
          f"measures contact POSITION, not residue identity. Report conservation as "
          f"{N_WILD_TYPE_COPIES} of {len(RECEPTOR_COPIES)} wild-type copies, never as "
          f"{len(RECEPTOR_COPIES)} of {len(RECEPTOR_COPIES)}.")


if __name__ == "__main__":
    main()
