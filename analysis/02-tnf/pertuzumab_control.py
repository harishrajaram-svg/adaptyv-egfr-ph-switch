#!/usr/bin/env python3
"""Pertuzumab/HER2: the first SINGLE-POINT histidine campaign run through our filter.

WHY THIS ONE, AND WHY NOW. s17 found the filter's blind spot: every positive it has ever called
correctly carries THREE TO FIVE installed histidines, and given ONE (AF-M2631, H L102H) it called
a measured 30x switch a non-switch. s4 designs TWO sites. s16's benchmark is 279 variants and
EVERY ONE IS A SINGLE-POINT HISTIDINE MUTATION -- so it is a direct test of exactly that blind
spot, and it needs no measured ratios because s15 made the criterion binary.

Pertuzumab is the cheapest entry: 9 variants, 3 labelled hits, and a clean deposited complex.

\U0001F534 DIFFERENT TARGET, AND THAT IS THE POINT HERE. HER2, not TNF-alpha. s13 established that
backbone constraint transfers across scaffolds; this asks a different transferable question --
can the filter see a ONE-histidine switch at all? A failure here is about the filter, not about
HER2. A success does not validate it on TNF.

SOURCE. Labels and positions from Wei & Sulea's Table S2, SI_Tables_R1.xlsx, parsed to
external/sipHAB/benchmark_parsed.json. Underlying campaign: Kang et al. 2019. Structure: 1S78,
"Insights into ErbB signaling from the structure of the ErbB2-pertuzumab complex" -- chain A is
HER2, C the pertuzumab light chain, D the heavy chain (the deposit carries two copies; we use
the first).

NUMBERING: 1S78 IS DEPOSITED IN KABAT, and 8 of the 9 positions verify EXACTLY as the wild-type
residue the benchmark names -- H31 D, H53 N, H54 S, H96 L, H99 S, H101 D, L53 Y, L55 Y.
That 8-for-8 agreement is the evidence the schemes match.

⚠️ THE NINTH IS MAPPED, NOT MATCHED, AND IT IS FLAGGED. The benchmark says "H100A, Y->H".
1S78 has no 100A; its CDR-H3 insertion run reads 99=SER, 99A=PHE, 99B=TYR, 100=PHE. The only
tyrosine in that run is 99B, so H100A -> (99,'B'). Kabat insertion codes are placed differently
by different tools and this is that difference. The residue IDENTITY matches (TYR), which is the
check available; it is recorded as an inference rather than a verification.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dddg_elec_repacked as dd
import fetch

BINDER, TARGET = "CD", "A"          # light + heavy vs HER2
CLEAN = os.path.join(HERE, "1S78_ACD.pdb")
BENCH = os.path.join(HERE, "external", "sipHAB", "benchmark_parsed.json")
# benchmark CDR position -> (1S78 chain, resnum, icode, expected wild-type one-letter)
MAP = {
    "H31":   ("D", 31, "", "D"), "H53": ("D", 53, "", "N"), "H54": ("D", 54, "", "S"),
    "H96":   ("D", 96, "", "L"), "H99": ("D", 99, "", "S"), "H101": ("D", 101, "", "D"),
    "H100A": ("D", 99, "B", "Y"),                      # INFERRED -- see the docstring
    "L53":   ("C", 53, "", "Y"), "L55": ("C", 55, "", "Y"),
}
THREE = {"D": "ASP", "N": "ASN", "S": "SER", "L": "LEU", "Y": "TYR"}
RENUM = {}        # (chain, kabat number, insertion code) -> sequential number after cleaning


def clean_structure():
    """Write chains A, C, D as protein-only. HER2 in 1S78 is glycosylated and the sugars are
    not part of the interface we score; pose_from_pdb would also choke on them."""
    import gemmi
    st = gemmi.read_structure(fetch.cif("1S78"))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_alternative_conformations()
    while len(st) > 1:
        del st[1]
    keep = {"A", "C", "D"}
    for ch in [c.name for c in st[0]]:
        if ch not in keep:
            st[0].remove_chain(ch)
    # \U0001F534 RENUMBER TO KILL INSERTION CODES. PDBFixer addresses residues by number only
    # and CANNOT see an insertion code: asking it for chain D residue 99 (SER) made it find
    # TYR, because 99, 99A and 99B all answer to "99". It RAISED rather than mutating the wrong
    # residue -- its own identity check caught it -- but the failure mode if the identities had
    # happened to agree is a silently mis-sited mutation, which is this project's recurring
    # shape (s21, now on a fourth numbering scheme). H99 is one of the three labelled HITS, so
    # that would have corrupted the result and not the run.
    #
    # Each chain is renumbered 1..N in order, and RENUM maps (chain, kabat_num, icode) -> new.
    for ch in st[0]:
        for i, r in enumerate(ch, 1):
            RENUM[(ch.name, r.seqid.num, r.seqid.icode.strip())] = i
    for ch in st[0]:
        for i, r in enumerate(ch, 1):
            r.seqid.num = i
            r.seqid.icode = " "
    st.setup_entities()
    st.write_pdb(CLEAN)
    return CLEAN


def verify(path):
    import gemmi
    st = gemmi.read_structure(path)
    idx = {(c.name, r.seqid.num, r.seqid.icode.strip()): r.name for c in st[0] for r in c}
    bad, inferred = [], []
    for key, (ch, num, ic, exp) in MAP.items():
        got = idx.get((ch, RENUM[(ch, num, ic)], ""))
        if got != THREE[exp]:
            bad.append(f"{key} -> {ch}{num}{ic}: {got}, expected {THREE[exp]}")
        elif ic:
            inferred.append(key)
    if bad:
        raise SystemExit("REFUSING, the Kabat mapping does not hold:\n  " + "\n  ".join(bad))
    print(f"  ok  T1 all {len(MAP)} benchmark positions carry the wild-type residue named, "
          f"{len(MAP) - len(inferred)} by direct number match and "
          f"{len(inferred)} by inferred insertion code ({', '.join(inferred)})")


def build(label, ch, num, ic, old, wt_path):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    out = os.path.join(HERE, f"1S78_{label}.pdb")
    fixer = PDBFixer(filename=wt_path)
    fixer.missingResidues = {}
    fixer.applyMutations([f"{THREE[old]}-{num}-HIS"], ch)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with open(out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return out


def main():
    bench = json.load(open(BENCH))["Pertuzumab (Her2)"]
    truth = {pos: hit for pos, mut, hit in bench["rows"]}
    print(f"Pertuzumab/HER2: {bench['n']} single-point His variants, {bench['hits']} labelled hits\n")
    wt = clean_structure()
    verify(wt)
    print()

    rows = [("WT", None, False)]
    for pos in truth:
        rows.append((pos, MAP.get(pos), truth[pos]))

    print(f"{'variant':<10}{'labelled':<10}{'rigid':>9}{'repacked_local med [min..max]':>32}  pred")
    print("-" * 82)
    res = {}
    for label, loc, hit in rows:
        if label == "WT":
            path = wt
        elif loc is None:
            print(f"{label:<10}  no mapping, skipped"); continue
        else:
            ch, num, ic, old = loc
            path = build(label, ch, RENUM[(ch, num, ic)], "", old, wt)
        try:
            r = dd.score(path, BINDER, TARGET, legs_wanted=("rigid", "repacked_local"))
        except SystemExit as e:
            print(f"{label:<10}  REFUSED: {str(e)[:70]}"); continue
        rl = r["repacked_local"]
        pred = rl["dddG"] >= 0
        res[label] = (hit, pred, rl["dddG"])
        flag = "!" if rl["sign_unstable"] else " "
        mark = "" if label == "WT" else ("  correct" if pred == hit else "  WRONG")
        print(f"{label:<10}{('HIT' if hit else '-'):<10}{r['rigid']['dddG']:>+9.3f}"
              f"{rl['dddG']:>+11.3f} [{rl['dddG_min']:+.2f}..{rl['dddG_max']:+.2f}]{flag}"
              f"   {'switch' if pred else 'non':<7}{mark}")

    v = {k: x for k, x in res.items() if k != "WT"}
    if v:
        ok = sum(1 for hit, pred, _ in v.values() if hit == pred)
        tp = sum(1 for hit, pred, _ in v.values() if hit and pred)
        fn = sum(1 for hit, pred, _ in v.values() if hit and not pred)
        fp = sum(1 for hit, pred, _ in v.values() if not hit and pred)
        print("-" * 82)
        print(f"\n{ok} of {len(v)} correct.  true positives {tp}, FALSE NEGATIVES {fn}, "
              f"false positives {fp}")
        print("\nEvery variant here installs ONE histidine, which is the regime s17 found the")
        print("filter failing in. False negatives are the number that matters.")


if __name__ == "__main__":
    main()
