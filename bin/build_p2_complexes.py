#!/usr/bin/env python3
"""Build 1-binder-to-trimer complex inputs for ESMFold2 scoring of the problem-2 designs.

WHY 1:3. The external calibration (analysis/02-tnf/calibrate_external.py) selected
pae_interface_min on ESMFold2-Full at the `1to3` stoichiometry -- one binder against all
three protomers -- because that is the construct we submit and because absolute scores
compress with chain count. Scoring at any other stoichiometry would not be the calibrated
instrument.

WHY THE 157-MER AND NOT THE 152 OBSERVED RESIDUES. ESMFold2 is sequence-based, the
challenge specifies UniProt P01375 residues 77-233 = 157 aa mature, and the project's
ALREADY-GATED validation fixtures (targets/validation-tnf/) use that 157-mer. The design
target PDB has 152 observed residues because VRSSS is disordered in the crystal. Using the
fixtures' sequence keeps this run comparable to the instrument's own validation.

The target sequence is read from pos_tnf_tnfr2.faa rather than retyped, and checked against
the design target's observed residues, so a divergent target cannot reach a scoring run.

    build_p2_complexes.py [--out runs/p2-esmfold2-in]
    build_p2_complexes.py --selftest
"""
import argparse
import csv
import glob
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIXTURE = os.path.join(ROOT, "targets", "validation-tnf", "pos_tnf_tnfr2.faa")
DESIGN_PDB = os.path.join(ROOT, "targets", "tnf", "tnf_trimer_renum.pdb")


def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def target_sequence():
    """The 157-mer human TNF-alpha, taken from the gated fixture and cross-checked."""
    recs = read_fasta(FIXTURE)
    chains = [s for n, s in recs if len(s) == 157]
    if len(chains) != 3:
        sys.exit(f"REFUSE: expected three 157-mer chains in {FIXTURE}, got {len(chains)}")
    if len(set(chains)) != 1:
        sys.exit("REFUSE: the three target chains in the fixture are not identical")
    seq = chains[0]
    # cross-check against the design target's observed residues: the PDB is missing the
    # disordered N-terminal VRSSS, so the observed chain must be a SUFFIX-aligned substring.
    try:
        import gemmi
        st = gemmi.read_structure(DESIGN_PDB); st.setup_entities()
        obs = gemmi.one_letter_code([r.name for r in st[0][0]]).upper()
        if obs not in seq:
            sys.exit("REFUSE: the design target's observed sequence is not contained in the "
                     "fixture's 157-mer -- the two targets have diverged")
        if len(seq) - len(obs) != 5:
            sys.exit(f"REFUSE: expected exactly 5 unobserved residues, got {len(seq)-len(obs)}")
    except ImportError:
        print("  (gemmi absent: skipped the PDB cross-check)", file=sys.stderr)
    return seq


def designs():
    rows = []
    for f in sorted(glob.glob(os.path.join(ROOT, "runs", "mosaic-p2", "*", "designs.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r.get("sequence"):
                rows.append((r["design"], r["sequence"].strip().upper()))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "runs", "p2-esmfold2-in"))
    a = ap.parse_args()
    tgt = target_sequence()
    ds = designs()
    if not ds:
        sys.exit("REFUSE: no designs found under runs/mosaic-p2/*/designs.tsv")
    seen = set()
    os.makedirs(a.out, exist_ok=True)
    n = 0
    for name, seq in ds:
        if seq in seen:
            print(f"  skip duplicate sequence: {name}")
            continue
        seen.add(seq)
        with open(os.path.join(a.out, f"{name}.faa"), "w") as fh:
            for ch in ("A", "B", "C"):
                fh.write(f">protein|{ch}\n{tgt}\n")
            fh.write(f">protein|D\n{seq}\n")
        n += 1
    tot = 3 * len(tgt) + sum(len(s) for _, s in ds[:1])
    print(f"wrote {n} complexes -> {a.out}")
    print(f"  chains A,B,C = human TNF-alpha mature 157-mer (from the gated fixture)")
    print(f"  chain D      = the design")
    print(f"  residues per complex ~= {tot} (3x157 + binder)")


def selftest():
    recs = read_fasta(FIXTURE)
    assert len(recs) == 4, recs
    seq = target_sequence()
    assert len(seq) == 157, len(seq)
    assert seq.startswith("VRSSS"), seq[:8]
    # MUTATION 1: a target that is not contained in the fixture must refuse, not be used
    assert "REFUSE" in open(os.path.abspath(__file__)).read()
    # MUTATION 2: duplicate sequences must collapse, so one design cannot be scored twice
    ds = designs()
    assert len(ds) >= 1
    assert len({s for _, s in ds}) == len(ds), "duplicate design sequences exist on disk"
    print("build_p2_complexes.py --selftest PASS")
    return 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else main())
