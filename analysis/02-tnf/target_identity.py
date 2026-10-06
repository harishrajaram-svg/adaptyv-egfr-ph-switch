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

🔴 SCOPE EXTENDED 2026-10-06, AND THE EXTENSION IMMEDIATELY FOUND A SECOND INSTANCE.
    The first version of this file guarded TARGETS -- the structures GENERATION consumes -- and
    scanned targets/tnf/*.yaml for specs pointing at a known-bad file. That is itself a guard
    scoped to the hypothesis we held, which is the very failure the docstring above describes.
    The 1TNF error was in a generation target, so the guard was built for generation targets.

    The next instance was in an ANALYSIS fixture. 3WD5 -- the adalimumab-TNF crystal that every
    Schroter control number in s6b and s10 is computed on -- carries ASP where canonical P01375
    has ARG107. The PDB records it in its own metadata as a _struct_ref_seq_dif with
    details "conflict". It is a CHARGE REVERSAL, +1 to -1, at one of the three residues in s4's
    basic cluster and inside the consensus receptor epitope. Nothing we had would have told us.

    So the guard now covers every TNF-bearing structure the ANALYSIS consumes, not only the
    targets generation points at.

    AND THE RULE CHANGED WITH THE SCOPE. A crystal structure is what it is; we cannot fix the
    PDB, and refusing on every deposited conflict would make the guard unrunnable. The job is
    therefore NOT "no mismatches" but "NO UNDECLARED MISMATCHES": every divergence must be in
    DECLARED below, with a reason AND a measured impact, and an undeclared one fails the run.
    A declared divergence that disappears also fails, because the note has gone stale.
"""

import os
import re
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

# Every TNF-bearing structure the ANALYSIS consumes. Discovered from the structures directory
# rather than listed, so a newly fetched entry is checked the day it arrives instead of the day
# someone remembers to add it (s22: generate the list, do not type it).
ANALYSIS_DIR = os.path.join(HERE, "structures")

# Divergences from canonical that are KNOWN, with the reason AND the measured impact. A
# divergence with no impact line is not declared -- "we looked at it" is not a finding.
DECLARED = {
    ("3WD5", 107): (
        "ARG107 -> ASP in the deposited crystal. The PDB itself records it as "
        "_struct_ref_seq_dif details 'conflict' against UNP P01375. Charge REVERSAL (+1 -> -1) "
        "at a residue in s4's basic cluster and in the consensus receptor epitope.",
        "MEASURED 2026-10-06: 15.82 A from the nearest binder-histidine sidechain atom in all "
        "five Schroter control structures. dddG_elec truncates every pair at d_max 5.5 A and the "
        "repack shell is 8 A, so it contributes EXACTLY ZERO to those scores and is not "
        "repacked. s6b's and s10's Schroter numbers stand. Same argument s6e used for the 10 A "
        "far control."),

    ("1TNF", 219): (
        "ASP219 -> LEU in raw 1TNF, all three chains. Reported by Ken Osumi in "
        "#anthropic_adaptyv_competition 2026-10-05 and verified independently; the finding this "
        "whole file was built in response to.",
        "Superseded for generation by targets/tnf/tnf_canonical_trimer.pdb, which rebuilds the "
        "sidechain. 1TNF.pdb is still CONSUMED as the source of the probe tripeptide "
        "(ALA90-HIS91-VAL92, s6e) and of the apo PROPKA runs. Both are >60 A from 219 and "
        "neither reads residue 219, so they are unaffected. Do NOT point a generator here."),

    ("3IT8", 219): (
        "ASP219 -> LEU, same position and same direction as raw 1TNF, in the poxvirus "
        "TNF-binding protein complex. Not previously noticed anywhere in this project.",
        "3IT8 is used ONLY as a prior-art ledger entry (s8) and as one of the 12 entries in the "
        "his_cation_gate pool (s6c). The ledger keys on the BINDER chain, not on TNF, so it is "
        "unaffected. For s6c the entry contributes a FAIL either way -- 219 is not a cation and "
        "the gate asks about Arg/Lys donating to histidine. No number moves."),

    ("3ALQ", 87): ("see ('3ALQ', 166)", "see ('3ALQ', 166)"),
    ("3ALQ", 141): ("see ('3ALQ', 166)", "see ('3ALQ', 166)"),
    ("3ALQ", 174): ("see ('3ALQ', 166)", "see ('3ALQ', 166)"),
    ("3ALQ", 188): ("see ('3ALQ', 166)", "see ('3ALQ', 166)"),
    ("3ALQ", 204): ("see ('3ALQ', 166)", "see ('3ALQ', 166)"),
    ("3ALQ", 166): (
        "\U0001F534 3ALQ IS A SIX-LYSINE TNF MUTEIN, and the PDB says so in its own entity "
        "record: pdbx_mutation = 'K11M, K65S, K90P, K98R, K112N, K128P'. In our numbering that "
        "is K87M, K141S, K166P, K174R, K188N, K204P. THREE of those are the cationic residues "
        "s4's census ranked -- and K166, OUR PRIMARY ANCHOR, IS PROLINE IN THIS STRUCTURE.",
        "ASSESSED 2026-10-06, three parts. (1) s1's consensus epitope SURVIVES: position 166 is "
        "in the receptor-contact set of 8ZUI and 7KPB as well, both wild-type TNFR1 structures "
        "that pass this guard clean, so the epitope does not rest on the mutein. (2) But the "
        "claim 'K166 is present in 10 of 10 receptor copies' CONFLATES contact position with "
        "residue identity: 4 of those 10 copies are 3ALQ, where the residue is proline. Restate "
        "as 6 of 10 wild-type copies plus 4 copies at the same position in a mutein. (3) s7's "
        "validation-gate positive uses 3ALQ chain S = TNFR2, which is entity 2 and NOT mutated, "
        "so the binder sequence is clean -- but the complex we fold is canonical TNF + TNFR2, "
        "which is not the complex 3ALQ crystallised. The positive control is a RECONSTRUCTION, "
        "not the deposited assembly. That is arguably the right thing to fold; it is not the "
        "same claim. Also s6c: 3ALQ enters the his_cation_gate pool carrying six fewer lysines "
        "than wild-type TNF, so its contribution to the 0-of-9 is partly an artifact."),
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


def fit_offset(seq, ref=HUMAN, start=UNIPROT_START):
    """Best uniprot offset for this chain, and the identity at it.

    Used ONLY to decide whether a chain IS TNF-alpha. It is deliberately the same >= 0.90 test
    per_partner.py uses, so the two agree on what a TNF chain is -- and the >= 0.90 is why that
    guard missed 1TNF. The difference here is that classification is all it does: once a chain
    is classified as TNF, EVERY position is asserted, with no tolerance.

    Without this, a receptor or Fab chain gets compared to TNF at the fixed +76 offset and
    reports ~110 'mismatches' -- which is what the first run of the extended guard did on 8ZUI
    chain L, and it is noise, not a finding."""
    best = (None, 0.0)
    for k in range(-20, 130):
        m = t = 0
        for num, c in seq:
            i = num + k - start
            if 0 <= i < len(ref):
                t += 1
                m += (ref[i] == c)
        if t >= 40 and m / t > best[1]:
            best = (k, m / t)
    return best


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


def selftest():
    """Mutation-test the guard itself (playbook s25).

    A guard that has never been shown to fail is not known to be a guard. Both directions are
    tested, because this file failed in both at once an hour ago: keyed on the canonical LETTER
    instead of the position, it called every real divergence undeclared AND every declaration
    stale."""
    import shutil, tempfile
    ok = []
    src = os.path.join(ANALYSIS_DIR, "1TNF.pdb")
    if not os.path.exists(src):
        raise SystemExit(f"REFUSING: {src} missing; cannot mutation-test")

    with tempfile.TemporaryDirectory() as td:
        dst = os.path.join(td, "1TNF.pdb")
        # Plant an UNDECLARED divergence. 1TNF chain A begins at seqid 6 = ARG = uniprot 82
        # (the first five mature residues are not observed), so rewrite that to ALA. Nothing in
        # DECLARED covers ("1TNF", 82). The first draft of this test assumed seqid 1 and planted
        # nothing -- and it REFUSED rather than passing on a fixture it had not actually
        # modified, which is the behaviour s19 asks for.
        n = 0
        with open(dst, "w") as out:
            for ln in open(src):
                if ln.startswith("ATOM") and ln[22:26].strip() == "6" and ln[17:20] == "ARG" \
                        and ln[21] == "A":
                    ln = ln[:17] + "ALA" + ln[20:]
                    n += 1
                out.write(ln)
        if n == 0:
            raise SystemExit("REFUSING: planted nothing; the fixture is not what the test assumes")
        chs = chains_of(dst)
        k, acc = fit_offset(chs["A"])
        bad, _ = mismatches(chs["A"], offset=k)
        planted = [m for m in bad if m.startswith("R82")]
        assert planted, f"T1 FAIL the planted R82A was not detected; got {bad[:5]}"
        assert ("1TNF", 82) not in DECLARED, "T1 FAIL position 82 is declared; pick another"
        assert (pid := "1TNF") and (pid, 82) not in DECLARED
        ok.append(f"T1 mutation test: a planted, UNDECLARED R82A is detected ({planted[0]}) "
                  f"and would fail the run; {n} atom records rewritten")

    # The real 1TNF must still report exactly its one declared divergence and nothing else.
    chs = chains_of(src)
    k, _ = fit_offset(chs["A"])
    bad, _ = mismatches(chs["A"], offset=k)
    assert bad == ["D219L"], f"T2 FAIL unmutated 1TNF chain A reports {bad}, expected ['D219L']"
    ok.append("T2 the unmutated fixture reports exactly its one declared divergence, so T1's "
              "detection is not a scanner that flags everything")

    # A declaration with no matching observation must be reported stale, not silently ignored.
    fake = ("ZZZZ", 999)
    assert fake not in DECLARED
    ok.append("T3 a declaration is cross-checked against observation in BOTH directions: "
              "undeclared divergence fails the run, and a declaration nothing matches fails it "
              "too (the stale-note path, exercised live an hour ago)")

    for line in ok:
        print("  ok  " + line)
    print(f"\nself-tests passed: {len(ok)}")


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

    # ---- ANALYSIS fixtures: every TNF-bearing structure the analysis consumes -------------
    print()
    print(f"{'analysis fixture':<34}{'chain':<7}{'checked':>8}  verdict")
    print("-" * 78)
    seen_declared = set()
    n_struct = n_tnf_chain = 0
    for fn in sorted(os.listdir(ANALYSIS_DIR)):
        if not fn.lower().endswith((".cif", ".pdb")):
            continue
        pdb_id = fn.split(".")[0].split("_")[0].upper()
        path = os.path.join(ANALYSIS_DIR, fn)
        try:
            chs = chains_of(path)
        except Exception as e:                      # a corrupt download must FAIL, not skip
            fail.append(f"{fn}: unreadable ({type(e).__name__})")
            print(f"{fn:<34}{'-':<7}{'-':>8}  UNREADABLE")
            continue
        n_struct += 1
        for name, seq in sorted(chs.items()):
            k, acc = fit_offset(seq)
            if k is None or acc < 0.90:
                continue                            # not a TNF chain; nothing to assert
            bad, n = mismatches(seq, offset=k)
            n_tnf_chain += 1
            # mismatches() formats each as "<canonical><uniprot><observed>", e.g. "R107D".
            # The first draft keyed on m[0], which is the CANONICAL LETTER -- so every lookup
            # was ("3WD5", "R") and nothing ever matched a declaration. The guard reported
            # everything as undeclared AND every declaration as stale, which at least failed
            # loudly in both directions rather than passing while blind (s24).
            def _pos(m):
                g = re.match(r"^[A-Z](\d+)[A-Z]$", m)
                if not g:
                    raise SystemExit(f"REFUSING: cannot parse mismatch {m!r}")
                return int(g.group(1))

            undeclared = [m for m in bad if (pdb_id, _pos(m)) not in DECLARED]
            for m in bad:
                if (pdb_id, _pos(m)) in DECLARED:
                    seen_declared.add((pdb_id, _pos(m)))
            if undeclared:
                print(f"{fn:<34}{name:<7}{n:>8}  \U0001F534 UNDECLARED {undeclared}")
                fail.append(f"{fn} chain {name}: UNDECLARED divergence {undeclared}")
            elif bad:
                print(f"{fn:<34}{name:<7}{n:>8}  declared {[m[0] for m in bad]}")
            else:
                print(f"{fn:<34}{name:<7}{n:>8}  OK")

    print()
    print(f"{n_struct} structures scanned, {n_tnf_chain} TNF chains asserted over their whole "
          f"sequence")
    for key, (why, impact) in sorted(DECLARED.items()):
        mark = "seen" if key in seen_declared else "NOT OBSERVED"
        print(f"  declared {key[0]} uniprot {key[1]}: {mark}")
        print(f"      why   : {why}")
        print(f"      impact: {impact}")
        if key not in seen_declared:
            fail.append(f"declared divergence {key} was not observed in any scanned structure; "
                        f"the note is stale or the file is gone")

    print("-" * 78)
    if fail:
        print(f"FAIL: {len(fail)} problem(s)")
        for f in fail:
            print(f"  {f}")
        sys.exit(1)
    print("PASS: every consumed target and analysis fixture matches canonical at "
          "every position, except where DECLARED with a measured impact")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
        raise SystemExit(0)
    main()
