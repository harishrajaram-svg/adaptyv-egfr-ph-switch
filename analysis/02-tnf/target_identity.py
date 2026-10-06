#!/usr/bin/env python3
"""
Assert that every target structure matches the CANONICAL sequence at EVERY position.

WHY THIS EXISTS, AND WHY THE EXISTING GUARDS MISSED IT
    On 2026-10-06 Ken Osumi reported in #anthropic_adaptyv_competition that 1TNF carries LEU at
    PDB residue 143 in all three chains where canonical human TNF-alpha (UniProt P01375) has
    ASP -- Asp219 of the full-length protein, uniprot 219 in our numbering. One mismatch in 157
    residues, and it sits in the CONSENSUS CORE of the receptor-contact epitope, present in
    10 of 10 receptor copies.

    He attached the consequence rather than just the observation: binders he designed against
    1TNF failed re-prediction against canonical TNF-alpha, and after redesigning against the
    canonical sequence they bound by BLI at Adaptyv Bio.

    Three of our own guards were in a position to catch this and all three declined:
      per_partner.py  best_offset() accepts a chain at >= 0.90 identity. 156/157 is 0.994.
      third_site_census.py  T3 asserts residue identity ONLY at CATIONIC positions.
      anchor_reach.py       T3 asserts identity ONLY at the five positions it reports on.
    Every one of them verified identity exactly where we were already looking. The residue we
    were not looking at is the one that was wrong. That is playbook s24 in a new costume: a
    guard scoped to the subset you suspect cannot tell you about the subset you do not.

    So this file checks the WHOLE sequence, every chain, every target file, and refuses rather
    than warns. It is cheap and it is the only guard here that is not scoped to a hypothesis.

    One mismatch is all it takes, so there is no tolerance parameter on purpose.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

from species_and_histidines import HUMAN, UNIPROT_START   # noqa: E402

OFFSET = 76          # uniprot = pdb_seqid + OFFSET, verified residue-by-residue

# Target files that generation and scoring actually consume, and what each must equal.
TARGETS = {
    "targets/tnf/tnf_canonical_trimer.pdb": ("human", 3),
}

# Files kept for provenance that are KNOWN to differ. Listing them here is a statement that
# the difference is understood, not that it is acceptable to design against them.
KNOWN_BAD = {
    "targets/tnf/tnf_1tnf_trimer.pdb":
        "raw 1TNF: LEU at uniprot 219 where canonical has ASP, all three chains. "
        "Superseded by tnf_canonical_trimer.pdb. Do not point a generator at this.",
}


def chains_of(path):
    import gemmi
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    out = {}
    for ch in st[0]:
        seq = []
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info and info.is_amino_acid():
                seq.append((r.seqid.num, info.one_letter_code.upper()))
        if seq:
            out[ch.name] = seq
    return out


def mismatches(seq, ref=HUMAN, start=UNIPROT_START, offset=OFFSET):
    bad, checked = [], 0
    for num, c in seq:
        u = num + offset
        i = u - start
        if 0 <= i < len(ref):
            checked += 1
            if ref[i] != c:
                bad.append(f"{ref[i]}{u}{c}")
    return bad, checked


def main():
    fail = []
    print(f"{'file':<40}{'chain':<7}{'checked':>8}  verdict")
    print("-" * 78)
    for rel, (species, want_chains) in TARGETS.items():
        path = os.path.join(ROOT, rel)
        if not os.path.exists(path):
            print(f"{rel:<40}{'-':<7}{'-':>8}  MISSING")
            fail.append(f"{rel}: file not found")
            continue
        chs = chains_of(path)
        if len(chs) != want_chains:
            fail.append(f"{rel}: {len(chs)} chains, expected {want_chains}")
        for name, seq in sorted(chs.items()):
            bad, n = mismatches(seq)
            verdict = "OK" if not bad else f"MISMATCH {bad}"
            print(f"{os.path.basename(rel):<40}{name:<7}{n:>8}  {verdict}")
            if bad:
                fail.append(f"{rel} chain {name}: {bad}")

    print()
    for rel, why in KNOWN_BAD.items():
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            continue
        chs = chains_of(p)
        got = {m for seq in chs.values() for m in mismatches(seq)[0]}
        print(f"known-divergent, retained for provenance: {os.path.basename(rel)}")
        print(f"    {why}")
        print(f"    observed now: {sorted(got) or 'none -- the note is STALE, update it'}")
        if not got:
            fail.append(f"{rel} is listed as known-divergent but now matches canonical; "
                        f"the KNOWN_BAD note is stale and must be corrected")

    # Nothing under targets/tnf/ may point a generator at a known-divergent file.
    print()
    tdir = os.path.join(ROOT, "targets", "tnf")
    bad_refs = []
    for fn in sorted(os.listdir(tdir)):
        if not fn.endswith(".yaml"):
            continue
        for ln_no, ln in enumerate(open(os.path.join(tdir, fn)), 1):
            s = ln.strip()
            if s.startswith("#"):
                continue                      # a comment may name the bad file; a path may not
            for bad in KNOWN_BAD:
                if os.path.basename(bad) in s:
                    bad_refs.append(f"{fn}:{ln_no} -> {os.path.basename(bad)}")
    if bad_refs:
        print("GENERATOR SPECS POINTING AT A KNOWN-DIVERGENT TARGET:")
        for b in bad_refs:
            print(f"    {b}")
        fail.extend(bad_refs)
    else:
        print(f"no generator spec in targets/tnf/ points at a known-divergent target "
              f"({sum(1 for f in os.listdir(tdir) if f.endswith('.yaml'))} specs checked)")

    print("-" * 78)
    if fail:
        print(f"FAIL: {len(fail)} problem(s)")
        for f in fail:
            print(f"  {f}")
        sys.exit(1)
    print("PASS: every consumed target matches canonical at every position")


if __name__ == "__main__":
    main()
