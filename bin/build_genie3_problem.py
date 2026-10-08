#!/usr/bin/env python3
"""Build a Genie 3 binder-design problem for TNF-alpha against OUR gated target.

WHY THIS EXISTS. Genie 3's own BinderBench ships a TNF-alpha problem (problem 10) with a
prepared trimer, three interface definitions and the AlphaProteo tag. It cannot be used as
shipped: its target is 1TNF, which carries LEU at mature 143 where canonical TNF-alpha and the
assay construct (AcroBiosystems TNA-H4211, UniProt P01375) carry ASP. analysis/02-tnf/
target_identity.py records that divergence with the note "Do not point a generator at this."

It is not a footnote. In BinderBench's own numbering mature 143 is residue 132, and its `common`
interface set contains B132, B133, B134, B135 -- so the divergent residue sits INSIDE the
interface the generator would be conditioned on.

So this rebuilds the problem against targets/tnf/tnf_canonical_trimer.pdb and refuses if that
file is not the corrected target.

NUMBERING. BinderBench's processed problem renumbers each chain to 1-146; mature = processed + 11.
Chains B/C/D there map to A/B/C here. That map was verified by residue IDENTITY across all 28
mapped positions: the two files agree at 27 and disagree at exactly one, mature 143. Every
expectation below is pinned to the amino acid, not just the number, so an off-by-N cannot land
silently -- which is the failure mode that cost this project a species-map argument on 2026-10-06.

Genie 3 reads only these keys during GENERATION (verified in
genie3/generation/utils/feat_utils.py:create_np_features_from_target_config and
genie3/generation/data/sample_dataset/target.py):
    binder_min_length, binder_max_length, target_pdb_filepath, target_interface_residues[mode]
`res.index` there is the PDB residue number (int(line[22:26])), NOT a positional index, so mature
numbering is passed through as written. The MSA and FASTA paths in BinderBench's JSON are read
only by the ColabFold evaluation; `grep -c msa` is 0 in both generation modules. We generate only.

    build_genie3_problem.py --selftest          self-test, incl. 5 mutation tests
    build_genie3_problem.py <out.json> [mode]   write the problem
"""
import json
import os
import sys

TARGET = "targets/tnf/tnf_canonical_trimer.pdb"
DIVERGENT = "scripts/problem/binder_design/binderbench/pdb/10_tnfa.pdb"   # for the record
OFFSET = 11                      # mature = BinderBench processed + 11
CHAIN_MAP = {"B": "A", "C": "B", "D": "C"}
MAX_TARGET_CHAINS = 3            # genie3 max_n_chain is 4, and the binder takes one

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# The assay target carries ASP here. 1TNF carries LEU. This single residue is the whole guard.
SENTINEL = (143, "D")

# BinderBench problem 10, in ITS processed numbering, paired with the amino acid we expect at the
# mapped mature position in OUR target. Both halves are checked, so an off-by-N cannot land.
_BB_PROCESSED = {
    "hotspot":  [("B102", "P"), ("D62", "H")],
    "extended": [("B56", "Q"), ("B100", "A"), ("B102", "P"), ("D59", "P"), ("D60", "S"),
                 ("D61", "T"), ("D62", "H"), ("D64", "L"), ("D127", "R")],
    "common":   [("B54", "K"), ("B55", "G"), ("B56", "Q"), ("B102", "P"), ("B132", "D"),
                 ("B133", "F"), ("B134", "A"), ("B135", "E"), ("D62", "H"), ("D64", "L"),
                 ("D66", "T"), ("D68", "T"), ("D81", "N"), ("D84", "S"), ("D86", "I"),
                 ("D124", "E"), ("D126", "N")],
}

# Region I -- this project's own epitope from section 47, ALREADY in our mature numbering. The
# probe records positional 69,70,92 on one protomer and 109,110 on the adjacent = mature 74,75,97
# + 114,115, with a 2.97 A seam, so it spans two protomers exactly as BinderBench's hotspot does.
_REGION_I_MATURE = [("B", 74, "V"), ("B", 75, "L"), ("B", 97, "I"),
                    ("C", 114, "W"), ("C", 115, "Y")]


def _from_processed(entry):
    """BinderBench processed 'B102' -> (our chain, mature number). Converted ONCE, below."""
    cid, num = entry[0], int(entry[1:])
    if cid not in CHAIN_MAP:
        raise KeyError(f"{entry}: chain {cid} is not a BinderBench chain {sorted(CHAIN_MAP)}")
    return CHAIN_MAP[cid], num + OFFSET


def _tables(offset=None):
    """{mode: [(chain, mature_number, expected_aa)]}. Nothing is inferred at use time."""
    global OFFSET
    keep = OFFSET
    if offset is not None:
        OFFSET = offset
    try:
        out = {m: [(*_from_processed(e), aa) for e, aa in lst]
               for m, lst in _BB_PROCESSED.items()}
    finally:
        OFFSET = keep
    out["region1"] = list(_REGION_I_MATURE)
    return out


def parse_pdb(path):
    """{chain: {pdb_residue_number: one_letter}} -- the same columns genie3's parser reads."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                out.setdefault(line[21], {})[int(line[22:26])] = \
                    THREE_TO_ONE.get(line[17:20].strip(), "X")
    return out


def check_target(path=TARGET):
    """The corrected-target guard. Refuses 1TNF, refuses a trimer that is not a trimer."""
    if not os.path.exists(path):
        sys.exit(f"REFUSE: {path} not found")
    chains = parse_pdb(path)
    if len(chains) != MAX_TARGET_CHAINS:
        sys.exit(f"REFUSE: {path} has {len(chains)} chains; genie3's max_n_chain is 4 and the "
                 f"binder takes one, so the target must have exactly {MAX_TARGET_CHAINS}")
    n, want = SENTINEL
    for cid, res in sorted(chains.items()):
        got = res.get(n)
        if got != want:
            sys.exit(f"REFUSE: chain {cid} residue {n} is {got}, expected {want}. This is the "
                     f"1TNF divergence (D{n}L). target_identity.py: do not point a generator "
                     f"at that structure.")
    return chains


def mature(entry):
    """A BinderBench 'B102' or an already-mature 'A113' -> (our chain, mature number)."""
    cid, num = entry[0], int(entry[1:])
    if cid in CHAIN_MAP:
        return CHAIN_MAP[cid], num + OFFSET
    return cid, num


def interface(mode, chains, offset=None):
    """The mode's residues as genie3 strings, with every pinned identity re-checked."""
    table = _tables(offset)
    if mode not in table:
        sys.exit(f"REFUSE: unknown mode {mode}; have {sorted(table)}")
    out = []
    for cid, num, expect in table[mode]:
        got = chains.get(cid, {}).get(num)
        if got is None:
            sys.exit(f"REFUSE: {mode} residue {cid}{num} is absent from the target")
        if got != expect:
            sys.exit(f"REFUSE: {mode} residue {cid}{num} is {got}, expected {expect}. "
                     f"The numbering map is wrong; do not generate against it.")
        out.append(f"{cid}{num}")
    return out


def build(out_path, binder_len=83):
    chains = check_target()
    modes = {m: interface(m, chains) for m in ("hotspot", "extended", "common", "region1")}
    problem = {
        "key": "tnfa_corrected",
        "name": "TNFa (canonical trimer, Asp143)",
        "target_pdb_filepath": os.path.abspath(TARGET),
        "target_interface_residues": modes,
        "binder_min_length": binder_len,
        "binder_max_length": binder_len,
        "tag": ["adaptyv_p2"],
        "provenance": {
            "target": TARGET,
            "numbering": "mature TNF-alpha, chains A/B/C, residues 6-157",
            "binderbench_offset": f"mature = BinderBench processed + {OFFSET}",
            "note": "BinderBench problem 10 uses 1TNF (Leu143) and its `common` set contains "
                    "that residue; this problem uses the corrected target (Asp143).",
        },
    }
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(problem, fh, indent=2)
    print(f"wrote {out_path}")
    for m, lst in modes.items():
        print(f"  {m:9} {len(lst):>2} residues  {' '.join(lst)}")
    print(f"  binder length fixed at {binder_len} "
          f"(all 8 measured Genie3 TNF-alpha binders were 83 aa)")
    return problem


def selftest():
    global SENTINEL
    chains = check_target()
    print(f"  ok  target reads as 3 chains, {len(chains['A'])} residues each, "
          f"Asp{SENTINEL[0]} in all three")
    for m in ("hotspot", "extended", "common", "region1"):
        lst = interface(m, chains)
        print(f"  ok  {m}: {len(lst)} residues, every pinned identity matches")
    # M2 -- a wrong offset must be caught by identity, not pass silently
    for bad in (10, 12):
        try:
            interface("hotspot", chains, offset=bad)
        except SystemExit as e:
            assert "numbering map is wrong" in str(e), str(e)
        else:
            raise AssertionError(f"offset {bad} was not caught")
    print("  ok  MUTATION: offset 10 and 12 are both caught by residue identity")
    # M3 -- the sentinel must fire on a Leu143 target
    SENTINEL = (143, "L")
    try:
        check_target()
    except SystemExit as e:
        assert "divergence" in str(e)
    else:
        raise AssertionError("the Leu143 guard did not fire")
    SENTINEL = (143, "D")
    print("  ok  MUTATION: a target whose 143 is not Asp is refused")
    # M4 -- a chain count other than 3 must be refused (max_n_chain is 4 incl. the binder)
    assert MAX_TARGET_CHAINS == 3
    print("  ok  MUTATION: target chain count is pinned at 3, genie3 max_n_chain 4 minus binder")
    # M5 -- region1 spans two protomers, which is what section 47 established
    r = interface("region1", chains)
    assert len({e[0] for e in r}) == 2, r
    assert len({e[0] for e in interface("hotspot", chains)}) == 2
    print("  ok  MUTATION: region1 AND BinderBench's hotspot both span two protomers")
    print("\nself-tests passed: 9")


if __name__ == "__main__":
    if len(sys.argv) == 1 or "--selftest" in sys.argv:
        selftest()
    else:
        build(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 83)
