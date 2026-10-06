#!/usr/bin/env python3
"""
The control s6e named and did not run: K166 -> Gln in silico, probe geometry held FIXED.

WHY THIS EXISTS
    s6e gave dddG_elec its first mechanism-present control and reached two findings. The
    `>= 0` threshold is dead -- a histidine with no cation within 10 A passes it 72% of the
    time, against R108's 70%. But the MAGNITUDE separates at K166 on four of five band legs,
    median +0.627 against the far control's +0.059 on the default leg.

    s6e then named its own weakness, item 4 of "What this changes":

      "The control that would close this properly: mutate K166 -> Gln in silico and re-place
       the same probe at the same geometry. That removes the cation while holding the local
       environment fixed; the far control at THR83 changes both at once, so some of the
       separation above could be environment rather than charge."

    It also recorded a named confound: the far-control patch is NOT the same site between the
    pre- and post-correction runs (A:THR83:CB before, A:GLN97:CG after), so "the correction
    improved the instrument" is not cleanly attributable. THIS control has no such confound.
    Same site, same pose, same scorer. One variable: the charge.

WHAT IS HELD FIXED, AND HOW THAT IS ENFORCED RATHER THAN PROMISED
    The mutation is applied to the FIVE POSE FILES ALREADY ON DISK, not to the target before
    re-placing the probe. pose_anchor_K166_NZ_{1..5}.pdb each hold target chains A,B,C plus the
    probe as chain Z. Mutating the pose means:

      - the probe's coordinates are BIT-IDENTICAL between conditions, by construction, not by
        re-running a placement search and hoping it lands the same way
      - every target heavy atom except A:166's sidechain is unchanged
      - the five poses stay the five poses, so the n=5 spread is comparable pose-for-pose

    T2 asserts all of that atom by atom. A rebuild that quietly moved a hydroxyl would
    otherwise read as a charge effect -- which is exactly the fault that cost s6b a result
    ("two parametrisation runs do not cancel": addHydrogens re-optimised the whole H-bond
    network and left a residual 2.5x the signal).

    T2 FIRED ON THE FIRST RUN, and the fix is why the comparison is now symmetric. Comparing
    the raw parent pose against its PDBFixer-rebuilt mutant failed: addMissingAtoms added an
    OXT to the probe's C-terminal VAL, so the two files differed by an atom that has nothing
    to do with the mutation. That is the OXT fault anchor_dddg.py's own docstring warns about,
    arriving from a new direction. So BOTH conditions now go through the IDENTICAL PDBFixer
    path and differ only in whether a mutation was requested -- one pipeline, two sequences,
    which is the generalisation of s6b's "one geometry, two charge vectors".

    The raw parent is still scored, as a third row, so the rebuild's own effect on the number
    is MEASURED rather than assumed negligible.

WHY GLUTAMINE
    Isosteric to lysine within ~0.3 A of sidechain length, retains H-bond donor and acceptor
    capacity, and carries NO formal charge. It removes the cation and changes as little else as
    a substitution can. s6e specifies Gln by name; this does not reopen the choice.

WHAT EACH OUTCOME MEANS -- written down before the numbers exist (playbook s13)
    K166Q collapses to the far control's level   the K166 separation IS the cation. dddG_elec
                                                 ranks on mechanism at this anchor, and s6e's
                                                 re-scoping to a relative ranker is correct.
    K166Q stays near wild-type K166              the separation is the local ENVIRONMENT, not
                                                 the charge. dddG_elec cannot be a ranker here
                                                 either, and s6e's surviving half dies too.
    K166Q lands in between                       partially charge-driven; report the fraction
                                                 and never quote the raw number alone.

    The sign matters as much as the magnitude: removing a cation should move dddG_elec DOWN,
    because a protonated histidine beside a cation is destabilised and that is the whole
    mechanism. An INCREASE on mutation would mean the score is not reading what we think.

Scoring is bin/dddg_elec.py UNMODIFIED, called as a subprocess -- same as anchor_dddg.py, so
there is no second implementation to drift.

    anchor_dddg_k166q.py --selftest
    anchor_dddg_k166q.py            # build the mutants and score both conditions
"""
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

from third_site_census import tnf_chain_offset                    # noqa: E402

ANCHOR_UNIPROT = 166
ANCHOR_CHAIN = "A"
WT_FROM, MUT_TO = "LYS", "GLN"
N_POSE = 5
LEGS = [("sigmoid", "HIE"), ("sigmoid", "HID"), ("10", "HIE"), ("6", "HIE"), ("80", "HIE")]
MOVE_TOL = 1e-3          # A. PDB files carry 3 decimals, so this is "the same number".
VENV = os.path.join(ROOT, ".venv", "bin", "python")
PY = VENV if os.path.exists(VENV) else sys.executable


def pose_path(i):
    return os.path.join(HERE, f"pose_anchor_K166_NZ_{i}.pdb")


def wt_path(i):
    """The parent pose through the IDENTICAL PDBFixer path, no mutation requested. This is the
    wild-type condition, not the raw file -- see the T2 note in the module docstring."""
    return os.path.join(HERE, f"pose_anchor_K166_rebuilt_{i}.pdb")


def mut_path(i):
    return os.path.join(HERE, f"pose_anchor_K166Q_{i}.pdb")


def anchor_resnum():
    """A:166 in UniProt numbering -> the residue number the POSE FILE actually uses.

    Derived, never hardcoded. 1TNF numbers the mature chain so uniprot = seqid + 76, which
    makes K166 residue 90 -- but playbook s21 is on this page three times over (three live
    numbering schemes, and 'the fourth scheme' broke a gate line in the Mosaic work), so the
    offset is recomputed from the file in hand and the residue identity is asserted."""
    import gemmi
    st = gemmi.read_structure(pose_path(1))
    st.setup_entities()
    st.remove_ligands_and_waters()
    for ch in st[0]:
        if ch.name != ANCHOR_CHAIN:
            continue
        k, acc = tnf_chain_offset(ch)
        if k is None or acc < 0.90:
            raise SystemExit(f"REFUSING: chain {ch.name} does not read as TNF (acc={acc})")
        for r in ch:
            if r.seqid.num + k == ANCHOR_UNIPROT:
                if r.name != WT_FROM:
                    raise SystemExit(
                        f"REFUSING: uniprot {ANCHOR_UNIPROT} on chain {ANCHOR_CHAIN} is "
                        f"{r.name}, not {WT_FROM}. The offset or the target is wrong.")
                return r.seqid.num, k
    raise SystemExit(f"REFUSING: uniprot {ANCHOR_UNIPROT} not found on chain {ANCHOR_CHAIN}")


def heavy_atoms(path):
    """{(chain, resnum, resname, atomname): xyz} over heavy atoms."""
    import gemmi
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    out = {}
    for ch in st[0]:
        for r in ch:
            for at in r:
                if at.element == gemmi.Element("H"):
                    continue
                out[(ch.name, r.seqid.num, r.name, at.name)] = np.array(
                    [at.pos.x, at.pos.y, at.pos.z])
    return out


def build(i, resnum, mutate):
    """Run pose i through PDBFixer. `mutate` is the ONLY difference between the two conditions."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    src = pose_path(i)
    if not os.path.exists(src):
        raise SystemExit(f"REFUSING: {src} is missing. Run anchor_dddg.py first; this control "
                         f"must score the SAME poses, not freshly placed ones.")
    dest = mut_path(i) if mutate else wt_path(i)
    fixer = PDBFixer(filename=src)
    fixer.missingResidues = {}
    if mutate:
        fixer.applyMutations([f"{WT_FROM}-{resnum}-{MUT_TO}"], ANCHOR_CHAIN)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with open(dest, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    _normalise(dest, resnum)
    return dest


def _normalise(dest, resnum):
    """Delete any heavy atom the raw parent pose does not have, except at the mutated site.

    WHY. `addMissingAtoms` completes the probe's C-terminal VAL with an OXT the parent never
    had -- and it places that OXT 0.46 A apart between the two conditions, because the
    mutation perturbs the completion. T2 caught both facts in succession.

    Arguing the OXT away was tempting and wrong. dddG_elec cancels every cross-interface pair
    that touches no histidine to machine precision (its ZERO_TOL guard), so a displaced
    non-histidine atom looks harmless -- but TNF carries H91, H149 and H154, and a binder atom
    paired with a TARGET histidine does not cancel. The exemption would have been an argument
    about which pairs survive, made under deadline, about a score that has already produced one
    clean plausible wrong number per section of this file.

    So the atom set is restored to the parent's instead, which needs no argument: the parent
    atom set is the one s6e's published numbers were computed on, and the only difference
    between the two conditions becomes the mutated sidechain."""
    parent = set(heavy_atoms(pose_path(int(dest.rsplit("_", 1)[1].split(".")[0]))))
    keep, dropped = [], []
    for ln in open(dest).read().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            keep.append(ln)
            continue
        ch, num, rn, an = ln[21], int(ln[22:26]), ln[17:20].strip(), ln[12:16].strip()
        if an.startswith("H") and not an[:1].isdigit():
            continue                                  # heavy atoms only, as everywhere else
        at_site = (ch == ANCHOR_CHAIN and num == resnum)
        if at_site or any(k[0] == ch and k[1] == num and k[3] == an for k in parent):
            keep.append(ln)
        else:
            dropped.append(f"{ch}:{rn}{num}:{an}")
    with open(dest, "w") as fh:
        fh.write("\n".join(keep) + "\n")
    _normalise.last_dropped = dropped
    return dropped


def build_pair(i, resnum):
    return build(i, resnum, False), build(i, resnum, True)


def diff_report(i, resnum, left=None, right=None):
    """What changed between two files. Defaults to rebuilt-wild-type vs mutant, which is the
    comparison the control rests on; pass `left=pose_path(i)` to see the rebuild's own effect."""
    a = heavy_atoms(left or wt_path(i))
    b = heavy_atoms(right or mut_path(i))
    key_at_site = lambda k: k[0] == ANCHOR_CHAIN and k[1] == resnum
    moved, gone, added = [], [], []
    for k, p in a.items():
        if k not in b:
            gone.append(k)
        elif float(np.linalg.norm(p - b[k])) > MOVE_TOL:
            moved.append((k, float(np.linalg.norm(p - b[k]))))
    for k in b:
        if k not in a:
            added.append(k)
    return dict(moved=moved, gone=gone, added=added,
                moved_off_site=[m for m in moved if not key_at_site(m[0])],
                gone_off_site=[k for k in gone if not key_at_site(k)],
                added_off_site=[k for k in added if not key_at_site(k)],
                probe_changed=[m for m in moved if m[0][0] == "Z"]
                              + [k for k in gone if k[0] == "Z"]
                              + [k for k in added if k[0] == "Z"])


def score(path):
    vals = []
    for eps, taut in LEGS:
        r = subprocess.run([PY, os.path.join(ROOT, "bin", "dddg_elec.py"), path,
                            "--binder", "Z", "--target", "A,B,C",
                            "--eps", eps, "--tautomer", taut],
                           capture_output=True, text=True)
        v = None
        for ln in (r.stdout or "").splitlines():
            m = ln.strip()
            if m.startswith("dddG_elec"):
                try:
                    v = float(m.split()[1])
                except (ValueError, IndexError):
                    pass
        if v is None:
            tail = "\n    ".join(((r.stdout or "") + (r.stderr or "")).strip().splitlines()[-6:])
            print(f"    {eps}/{taut} returned no number (rc={r.returncode}):\n    {tail}")
        vals.append(v)
    return vals


def selftest():
    ok = []
    resnum, off = anchor_resnum()
    ok.append(f"T1 uniprot {ANCHOR_UNIPROT} on chain {ANCHOR_CHAIN} resolves to residue "
              f"{resnum} (offset +{off}) and is {WT_FROM}, asserted not assumed")

    build_pair(1, resnum)
    d = diff_report(1, resnum)

    # T2 -- the mutation changed the anchor site and NOTHING ELSE. This is the test the whole
    # control rests on: a rebuild that silently moved a hydroxyl elsewhere would read as a
    # charge effect (s6b bug 2, where re-parametrisation left a residual 2.5x the signal).
    assert not d["moved_off_site"], \
        f"T2 FAIL {len(d['moved_off_site'])} heavy atoms moved away from the site, worst " \
        f"{max(d['moved_off_site'], key=lambda m: m[1])}"
    assert not d["gone_off_site"], f"T2 FAIL atoms deleted off-site: {d['gone_off_site'][:4]}"
    assert not d["added_off_site"], f"T2 FAIL atoms added off-site: {d['added_off_site'][:4]}"
    ok.append(f"T2 off-site heavy atoms are bit-identical between the two conditions: "
              f"0 moved, 0 removed, 0 added (tolerance {MOVE_TOL} A)")

    # T3 -- the probe is untouched. If chain Z moved, the two conditions are not the same pose.
    assert not d["probe_changed"], f"T3 FAIL the probe changed: {d['probe_changed'][:4]}"
    ok.append("T3 the probe (chain Z) is unchanged -- same pose, not a re-placement")

    # T4 -- the site really did become a glutamine, and the cation really is gone.
    b = heavy_atoms(mut_path(1))
    site = {k[3]: v for k, v in b.items() if k[0] == ANCHOR_CHAIN and k[1] == resnum}
    names = {k[2] for k in b if k[0] == ANCHOR_CHAIN and k[1] == resnum}
    assert names == {MUT_TO}, f"T4 FAIL residue is {names}, expected {{{MUT_TO}}}"
    assert "NZ" not in site, "T4 FAIL the lysine NZ is still present"
    assert {"OE1", "NE2"} <= set(site), f"T4 FAIL glutamine amide missing, got {sorted(site)}"
    ok.append(f"T4 the site is {MUT_TO} with its amide (OE1, NE2) and no NZ -- the cation is gone")

    # T5 -- MUTATION TEST on the guard itself. Copy the parent over the mutant and confirm T4
    # fails: otherwise T4 passing is consistent with it inspecting the wrong file.
    keep = mut_path(1) + ".keep"
    shutil.copy(mut_path(1), keep)
    shutil.copy(wt_path(1), mut_path(1))
    try:
        names2 = {k[2] for k in heavy_atoms(mut_path(1))
                  if k[0] == ANCHOR_CHAIN and k[1] == resnum}
        assert names2 == {WT_FROM}, f"T5 FAIL planted parent reads {names2}"
        fired = (names2 != {MUT_TO})
    finally:
        shutil.move(keep, mut_path(1))
    assert fired, "T5 FAIL the residue-identity check cannot tell the files apart"
    ok.append("T5 mutation test: planting the parent over the mutant makes T4's check fail, "
              "so T4 is reading the file it claims to")

    # T6 -- the probe histidine survives as something dddG_elec will score. A mutant that
    # lost the interface histidine would return a clean, plausible, meaningless number.
    his = [k for k in b if k[0] == "Z" and k[2] == "HIS"]
    assert his, "T6 FAIL no histidine left on chain Z in the mutant"
    assert {"ND1", "NE2", "CG", "CD2", "CE1"} <= {k[3] for k in his}, \
        "T6 FAIL the probe imidazole is incomplete in the mutant"
    ok.append(f"T6 the probe histidine survives the rebuild intact ({len(his)} heavy atoms)")

    # T7 -- the rebuild's own footprint, MEASURED. This is what T2 could not see when the
    # comparison was parent-vs-mutant: PDBFixer adds an OXT to the probe's C-terminal VAL.
    # Both conditions now carry it, so it cancels -- but it is recorded, not waved away.
    rb = diff_report(1, resnum, left=pose_path(1), right=wt_path(1))
    assert not rb["added"], f"T7 FAIL normalisation left added atoms: {rb['added'][:4]}"
    off = [m for m in rb["moved"] if not (m[0][0] == ANCHOR_CHAIN and m[0][1] == resnum)]
    assert not off, f"T7 FAIL the rebuild moved off-site atoms vs the raw parent: {off[:3]}"
    ok.append(f"T7 the rebuilt wild type is atom-for-atom the raw parent off-site "
              f"(normalisation dropped {len(getattr(_normalise, 'last_dropped', []))}: "
              f"{sorted(set(getattr(_normalise, 'last_dropped', [])))}), so the rebuild "
              f"contributes nothing of its own")

    for line in ok:
        print("  ok  " + line)
    print(f"\nself-tests passed: {len(ok)}")
    return resnum


def main():
    resnum, off = anchor_resnum()
    print(f"K166 -> Gln control. uniprot {ANCHOR_UNIPROT} = residue {resnum} on chain "
          f"{ANCHOR_CHAIN} (offset +{off}), {WT_FROM} -> {MUT_TO}.")
    print(f"Scoring the SAME {N_POSE} poses as s6e, mutated in place.\n")

    wt, mt, raw = {}, {}, {}
    for i in range(1, N_POSE + 1):
        if not os.path.exists(pose_path(i)):
            print(f"  pose {i}: parent missing, skipped")
            continue
        build_pair(i, resnum)
        d = diff_report(i, resnum)
        bad = d["moved_off_site"] + d["gone_off_site"] + d["added_off_site"]
        if bad or d["probe_changed"]:
            raise SystemExit(f"REFUSING pose {i}: the two conditions differ off-site "
                             f"({len(bad)} items, probe {len(d['probe_changed'])}). They would "
                             f"not differ only by the charge, so the comparison is void.")
        wt[i] = score(wt_path(i))
        mt[i] = score(mut_path(i))
        raw[i] = score(pose_path(i))
        print(f"  pose {i}: {len(d['moved']) + len(d['gone']) + len(d['added'])} heavy atoms "
              f"differ, all at the site; off-site 0; probe unchanged")

    hdr = "".join(f"{e}/{t}"[:12].ljust(12) for e, t in LEGS)
    print(f"\n{'':<16}{hdr}")
    print("-" * (16 + 12 * len(LEGS)))
    for lab, tab in (("K166 wild type, rebuilt", wt), ("K166Q no cation, rebuilt", mt),
                     ("K166 RAW parent (rebuild check)", raw)):
        print(f"{lab}")
        for i, vals in sorted(tab.items()):
            cells = "".join((f"{v:+.3f}" if v is not None else "ERR").ljust(12) for v in vals)
            print(f"  pose {i:<10}{cells}")

    print(f"\n{'leg':<16}{'K166 median':>14}{'K166Q median':>14}{'delta':>10}   reading")
    print("-" * 78)
    verdicts = []
    for j, (e, t) in enumerate(LEGS):
        a = [wt[i][j] for i in wt if wt[i][j] is not None]
        b = [mt[i][j] for i in mt if mt[i][j] is not None]
        if not a or not b:
            print(f"{e+'/'+t:<16}{'ERR':>14}{'ERR':>14}")
            continue
        ma, mb = float(np.median(a)), float(np.median(b))
        overlap = not (max(b) < min(a) or max(a) < min(b))
        read = ("ranges OVERLAP -- not separated" if overlap
                else "separated, and the mutation moves it "
                     + ("DOWN, as removing a cation should" if mb < ma else "UP, which is WRONG"))
        verdicts.append((e + "/" + t, ma, mb, mb - ma, overlap, mb < ma))
        print(f"{e+'/'+t:<16}{ma:>+14.3f}{mb:>+14.3f}{mb-ma:>+10.3f}   {read}")

    sep = [v for v in verdicts if not v[4]]
    right_way = [v for v in sep if v[5]]
    print(f"\n{len(sep)} of {len(verdicts)} legs separate K166 from K166Q; "
          f"{len(right_way)} of those move DOWN on losing the cation.")
    print("A leg that separates in the WRONG direction is evidence against the score reading\n"
          "the mechanism, not for it -- a protonated histidine beside a cation is destabilised,\n"
          "so removing the cation must lower dddG_elec.")

    # How much of any separation could be the rebuild rather than the charge? Measured, because
    # s6e's far control could not separate those two and said so.
    print(f"\n{'leg':<16}{'raw parent':>14}{'rebuilt WT':>14}{'rebuild shift':>15}")
    print("-" * 60)
    for j, (e, t) in enumerate(LEGS):
        r = [raw[i][j] for i in raw if raw[i][j] is not None]
        w = [wt[i][j] for i in wt if wt[i][j] is not None]
        if not r or not w:
            continue
        print(f"{e+'/'+t:<16}{float(np.median(r)):>+14.3f}{float(np.median(w)):>+14.3f}"
              f"{float(np.median(w)) - float(np.median(r)):>+15.3f}")
    print("If a rebuild shift is comparable to the K166 -> K166Q delta above, the mutation\n"
          "result is not attributable and must be reported as such.")
    return verdicts


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        selftest()
        print()
        main()
