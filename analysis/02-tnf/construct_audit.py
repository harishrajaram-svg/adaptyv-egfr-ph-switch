#!/usr/bin/env python3
"""Audit every target construct against the structures that were actually SOLVED.

WHY. Challenge §30 voided the G4 leg of `g-fab`: its LT-alpha carried 37 residues of signal
peptide AND was missing ~24 residues internally (120 where the mature domain is 144), so the
trimer never assembled and every inter-chain ipSAE read exactly 0.0000. §27 had already
diagnosed and fixed that construct -- by rebuilding it from the observed 1TNR chain -- but the
fix lived in one file and the run read another. **A fix is not applied until every live
consumer of the broken artifact is re-pointed at the fixed one.**

That means the interesting question is not "is G4 fixed" (it is) but "what ELSE is reading an
alignment-derived construct". This script answers it the cheap decisive way: build an index of
every chain sequence in every local observed structure, then ask of each target chain whether
it appears there. A construct that matches a solved chain is trustworthy by construction. One
that does not was written by hand or by alignment, and is where the §30 class of error lives.

VERDICTS
  EXACT      byte-identical to an observed chain -- the standard §27 set
  SUBSEQ     a contiguous slice of an observed chain, or vice versa -- fine, offset reported
  UNVERIFIED matches no observed chain. NOT necessarily wrong: designed binders are supposed
             to be novel. Judge by what the chain IS, which is why --targets-only exists.
  SIGNAL     an N-terminal hydrophobic run consistent with an uncleaved signal peptide. This
             is what G4 had. Reported independently of the match verdict.

usage:
  python3 analysis/02-tnf/construct_audit.py                 # TNF/LT-alpha constructs
  python3 analysis/02-tnf/construct_audit.py --all
  python3 analysis/02-tnf/construct_audit.py --selftest
"""
import glob
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STRUCT_DIRS = ("analysis/02-tnf/structures", "targets")
MIN_LEN = 30                       # below this, a match is not evidence of anything


def observed_chains():
    """{sequence: [labels]} for every chain in every local observed structure."""
    import gemmi
    idx = defaultdict(list)
    for d in STRUCT_DIRS:
        for f in sorted(glob.glob(str(ROOT / d / "*.cif")) + glob.glob(str(ROOT / d / "*.pdb"))):
            try:
                st = gemmi.read_structure(f)
                st.setup_entities()
                st.remove_ligands_and_waters()
                st.remove_alternative_conformations()
            except Exception:
                continue
            for ch in st[0]:
                s = gemmi.one_letter_code([r.name for r in ch]).upper()
                s = re.sub(r"[^A-Z]", "", s)
                if len(s) >= MIN_LEN:
                    idx[s].append(f"{Path(f).stem}:{ch.name}")
    return idx


def read_faa(path):
    """[(chain_id, sequence)] from a Boltz-style >protein|A fasta."""
    out, cid, buf = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if cid is not None:
                out.append((cid, "".join(buf)))
            cid, buf = line.lstrip(">").split("|")[-1].strip(), []
        elif line.strip():
            buf.append(line.strip())
    if cid is not None:
        out.append((cid, "".join(buf)))
    return out


def signal_peptide(seq, win=15, need=11):
    """True if the N-terminal `win` residues contain a hydrophobic run of `need`.

    Calibrated on the G4 failure: ERFLPRTHLLLLGLLLVLLPGAQ... carries LLLLGLLLVLL. This is a
    heuristic and is reported SEPARATELY from the match verdict -- it never silently rejects.
    """
    hydro = set("AVLIMFWCG")
    head = seq[:win + need]
    run = best = 0
    for c in head:
        run = run + 1 if c in hydro else 0
        best = max(best, run)
    return best >= need


def classify(seq, idx):
    if len(seq) < MIN_LEN:
        return "SHORT", ""
    if seq in idx:
        return "EXACT", idx[seq][0]
    for obs, labels in idx.items():
        if seq in obs:
            return "SUBSEQ", f"{labels[0]} [slice @{obs.index(seq)}, obs {len(obs)}]"
        if obs in seq:
            return "SUBSEQ", f"{labels[0]} [contains obs, +{len(seq) - len(obs)} extra]"
    return "UNVERIFIED", ""


def main():
    if "--selftest" in sys.argv:
        return selftest()
    want_all = "--all" in sys.argv
    idx = observed_chains()
    print(f"indexed {len(idx)} distinct observed chain sequences "
          f"from {STRUCT_DIRS[0]} and {STRUCT_DIRS[1]}\n")

    files = sorted(glob.glob(str(ROOT / "targets" / "**" / "*.faa"), recursive=True))
    if not want_all:
        files = [f for f in files
                 if re.search(r"tnf|lt-?a|lta|TNFR", Path(f).name, re.I)
                 or re.search(r"tnf|lta", str(Path(f).parent.name), re.I)]

    # A QUARANTINED file is already known-void (see targets/VOID-*/README.md), so it must not
    # inflate the live-risk count -- a tally that can never go green stops being read. It is
    # still listed, because silently skipping it is how a void construct becomes invisible.
    quarantined = [f for f in files if "VOID" in f]
    files = [f for f in files if "VOID" not in f]

    tally = defaultdict(int)
    flagged = []
    for f in files:
        rows = read_faa(f)
        notes = []
        for cid, seq in rows:
            verdict, where = classify(seq, idx)
            sig = signal_peptide(seq)
            tally[verdict] += 1
            if sig:
                tally["SIGNAL"] += 1
            if verdict == "UNVERIFIED" or sig:
                notes.append((cid, len(seq), verdict, sig, seq[:26]))
        if notes:
            flagged.append((f, notes))

    for f, notes in flagged:
        print(f"{Path(f).relative_to(ROOT)}")
        for cid, n, verdict, sig, head in notes:
            mark = "🔴 SIGNAL" if sig else "  "
            print(f"   {mark} chain {cid}  {n:4d} res  {verdict:10s}  {head}...")
    if quarantined:
        print("quarantined, excluded from the tally below:")
        for f in quarantined:
            print(f"   ⬛ {Path(f).relative_to(ROOT)}")
        print()
    print(f"{dict(tally)}")
    print(f"{len(files)} live files examined, {len(flagged)} carry at least one flagged chain"
          f"{f', {len(quarantined)} quarantined' if quarantined else ''}")
    if tally["SIGNAL"]:
        print("🔴 SIGNAL chains are the §30 failure mode. Verify each against an observed chain.")
    else:
        print("✅ no live construct carries a signal peptide.")


def selftest():
    # The G4 sequence that actually failed, and the 1TNR chain that replaced it.
    g4 = ("ERFLPRTHLLLLGLLLVLLPGAQGPGVGLPSAATAHLAHSTLKPAAHLGDPSSLLWRANTDRAFLQDGSLSNNS"
          "LLTSGIVYSQVVFSGKAYSPKATSSPLAHEVQLFSSQYPFHVPVYPGQEPWLHSMYHGAAQLTQGDQLDLVLSS"
          "TVFFGAFAL")
    good = ("KPAAHLIGDPSKQNSLLWRANTDRAFLQDGFSLSNNSLLVPTSGIYFVYSQVVFSGKAYSPKATSSPLYLAHEV"
            "QLFSSQYPFHVPLLSSQKMVYPGLQEPWLHSMYHGAAFQLTQGDQLSTHTDGIPHLVLSPSTVFFGAFAL")
    assert signal_peptide(g4), "the signal-peptide detector misses the sequence it was built on"
    assert not signal_peptide(good), "the CORRECTED construct must not be flagged"

    idx = {good: ["1TNR:A"]}
    assert classify(good, idx) == ("EXACT", "1TNR:A"), classify(good, idx)
    v, where = classify(g4, idx)
    assert v == "UNVERIFIED", (v, where)          # the broken one must NOT match
    # a slice of an observed chain is SUBSEQ, and the offset is reported
    v, where = classify(good[10:90], idx)
    assert v == "SUBSEQ" and "@10" in where, (v, where)
    # an observed chain plus an N-terminal extension is SUBSEQ, extension counted
    v, where = classify("MKVKL" + good, idx)
    assert v == "SUBSEQ" and "+5 extra" in where, (v, where)
    # MUTATION TEST: a detector that always returns True is useless -- a hydrophilic
    # N-terminus must come back clean, or SIGNAL carries no information.
    assert not signal_peptide("DEKRNQSTDEKRNQSTDEKRNQSTDEKRNQST"), "SIGNAL fires on anything"
    # and a short chain is not judged on a coincidental match
    assert classify("ACDEFG", idx)[0] == "SHORT"

    print("selftest  G4 flagged SIGNAL + UNVERIFIED; corrected construct EXACT and clean;")
    print("          slice->SUBSEQ@10; +5 extension counted; hydrophilic terminus clean (mutation)")
    print("selftest OK")


if __name__ == "__main__":
    main()
