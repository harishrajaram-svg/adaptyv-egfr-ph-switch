#!/usr/bin/env python3
"""A SECOND measured pH-switch series on TNF-alpha, from a second lab, scorable on a structure
we already hold.

Watkins & Watkins, J Immunol 209(4):829-839 (2022), PMID 35896334, PMC10580234,
doi 10.4049/jimmunol.2101180. Adafre Biosciences. Found 2026-10-06 in the reference list of the
SIpHAB compilation (s16) -- it is NOT in our 21-entry RCSB survey, NOT in the prior-art ledger,
and nothing in this project had named it.

WHY IT IS WORTH MORE THAN THE FOUR HETEROLOGOUS CAMPAIGNS s16 FOUND. Those are HER2, CTLA-4,
IL-6R and VEGF -- other targets, other epitopes. This is OUR TARGET, OUR DIRECTION (bind at 7.4,
release at 6.0), and it is built on ADALIMUMAB, so the structure is 3WD5: the same file s6b's
Schroter control already uses. Nothing has to be downloaded or modelled from scratch.

THE VARIANTS, mutations quoted from the paper:
  AF-M2637   L chain Q89H, R90H, N92H      measured SWITCH -- EC50 15x higher after a pH 6.0 wash
  AF-M2631   H chain L102H                 measured SWITCH -- EC50 30x higher. ONE histidine.
  AF-M2630   monovalent, NO His installed  measured NON-SWITCH, "no significant pH-dependent
                                           release". Its Fv is sequence-identical to adalimumab,
                                           so scoring it IS scoring wild type -- its control value
                                           lies in the FORMAT (monovalent), which an interface
                                           energy cannot see. Recorded, not scored.

\U0001F534 THE RATIOS ARE A THIRD QUANTITY. Schroter reports OFF-RATE (kd) ratios. The Wyman bound
constrains KD ratios. These are EC50 shifts after a wash in an ELISA-type format. They are not
interchangeable and this project has already paid once for conflating two of them (s-Prior-art
item 1). Treat the Adafre numbers as ORDINAL ONLY: both variants switch, wild type does not.

NUMBERING, CHECKED NOT ASSUMED. The paper gives positions without naming a scheme. 3WD5 chain L
carries GLN at 89, ARG at 90, ASN at 92, and chain H carries LEU at 102 -- all four are exactly
the wild-type residues the paper names, at exactly those numbers. That agreement across four
independent positions is itself the evidence that the schemes match; it is asserted in T1.

AF-M2631 IS THE MORE INTERESTING OF THE TWO. It installs a SINGLE histidine and reports the
LARGER shift (30x vs 15x). Every other positive this project has is multi-histidine: Schroter's
PSVs carry 3-5. A single-His positive tests whether our filter needs a network or can see one
site -- which is exactly what s4 designs, two sites rather than five.
"""
import os, sys, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import dddg_elec_repacked as dd

WT = os.path.join(HERE, "structures", "3WD5.pdb")
BINDER, TARGET = "HL", "A"
# (chain, residue number, expected wild-type residue)
VARIANTS = {
    "AF-M2637": ([("L", 89, "GLN"), ("L", 90, "ARG"), ("L", 92, "ASN")],
                 "SWITCH", "EC50 15x higher after pH 6.0 wash"),
    "AF-M2631": ([("H", 102, "LEU")],
                 "SWITCH", "EC50 30x higher; ONE histidine"),
}


def check_wt():
    """T1: every position the paper names must hold the residue the paper names."""
    import gemmi
    st = gemmi.read_structure(WT); st.setup_entities(); st.remove_ligands_and_waters()
    got = {}
    for ch in st[0]:
        for r in ch:
            got[(ch.name, r.seqid.num)] = r.name
    bad = []
    for name, (subs, _, _) in VARIANTS.items():
        for c, n, old in subs:
            if got.get((c, n)) != old:
                bad.append(f"{name}: {c}{n} is {got.get((c, n))}, paper says {old}")
    if bad:
        raise SystemExit("REFUSING, numbering schemes do not agree:\n  " + "\n  ".join(bad))
    n = sum(len(v[0]) for v in VARIANTS.values())
    print(f"  ok  T1 all {n} positions across {len(VARIANTS)} variants carry the wild-type "
          f"residue the paper names, in 3WD5's own numbering")


def build(name, subs):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    # structures/, not HERE: prior_art_ledger.py's REBUILT block reads from structures/ and the
    # first version wrote beside the script, so the two Adafre variants were added to the ledger
    # and silently not found (47 distinct sequences, unchanged).
    out = os.path.join(HERE, "structures", f"3WD5_{name.replace('-', '')}.pdb")
    fixer = PDBFixer(filename=WT)
    fixer.missingResidues = {}
    by_chain = {}
    for c, n, old in subs:
        by_chain.setdefault(c, []).append(f"{old}-{n}-HIS")
    for c, muts in by_chain.items():
        fixer.applyMutations(muts, c)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with open(out, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return out


def main():
    print(__doc__.strip().splitlines()[0] + "\n")
    check_wt()
    print()
    rows = {}
    # wild type through the SAME PDBFixer path, so the comparison is one pipeline two sequences
    wt_path = build("WT", [])
    rows["adalimumab_WT"] = ("NON-SWITCH", "EC50 232 pM, minimal pH sensitivity", wt_path)
    for name, (subs, truth, note) in VARIANTS.items():
        rows[name] = (truth, note, build(name, subs))

    print(f"{'variant':<16}{'measured':<12}{'rigid':>10}{'repacked_local med [min..max]':>34}  note")
    print("-" * 104)
    res = {}
    for name, (truth, note, path) in rows.items():
        r = dd.score(path, BINDER, TARGET, legs_wanted=("rigid", "repacked_local"))
        res[name] = (truth, r)
        rl = r["repacked_local"]
        flag = "!" if rl["sign_unstable"] else " "
        print(f"{name:<16}{truth:<12}{r['rigid']['dddG']:>+10.3f}"
              f"{rl['dddG']:>+12.3f} [{rl['dddG_min']:+.2f}..{rl['dddG_max']:+.2f}]{flag}"
              f"   n_his={r['n_his']}  {note}")

    print("-" * 104)
    print("\nAgainst the `>= 0` criterion (s15), on the repacked leg:")
    ok = 0
    for name, (truth, r) in res.items():
        pred = "SWITCH" if r["repacked_local"]["dddG"] >= 0 else "NON-SWITCH"
        hit = pred == truth
        ok += hit
        print(f"  {name:<16}measured {truth:<12}predicted {pred:<12}"
              f"{'correct' if hit else 'WRONG'}")
    print(f"\n  {ok} of {len(res)} correct on this series.")
    print("\n⚠️  EC50 shifts, not KD ratios and not off-rate ratios. Ordinal only: both engineered")
    print("   variants switch and wild type does not. Do not combine these magnitudes with")
    print("   Schroter's kd ratios.")


if __name__ == "__main__":
    main()
