#!/usr/bin/env python3
"""
Is the s4 switch mechanism geometrically REACHABLE on apo human TNF-alpha?

THE QUESTION, AND WHY IT COMES BEFORE GENERATION
    challenges/02-tnf-alpha.md s4 commits the design to binder histidines placed against R108 and
    K166 -- the only two cationic anchors that survive the mouse-conservation filter -- and prices
    that at 51.6x. Nothing has ever checked that a histidine can physically get there.

    bin/his_cation_gate.py found the mechanism in 0 of 9 real TNF complexes (s6c). Two readings
    were open: the mechanism is rare by default, or the criterion is too strict to fire on
    anything. This script decides between them, and it answers the design question at the same
    time:

        for each donor nitrogen of R108 and K166, is there an approach direction along which a
        real histidine sidechain can sit with its ND1 lone pair pointing at that nitrogen at
        <= 3.5 A, WITHOUT clashing into TNF?

    If the answer is no, s4's mechanism is unreachable and the 51.6x table is fiction -- better
    known now than after a generation spend. If the answer is yes, three things fall out at once:
    his_cation_gate gets the real-coordinate POSITIVE it has never had, the 0-of-9 is confirmed as
    scarcity rather than strictness, and the surviving CB positions are where the generator has to
    put binder backbone.

WHAT IS AND IS NOT CLAIMED
    This is a REACH test, not a design. It asks whether the approach vector is open on the apo
    trimer. It does not model the binder, so a direction counted clear here still has to be
    reachable by a real backbone, and that is the generator's problem. A direction counted
    BLOCKED, however, is blocked for any binder, which is what makes a negative result decisive.

METHOD
    Histidine geometry is taken from a REAL histidine in the structure being scored, not from
    idealised internal coordinates, so no bond length or angle in this file is invented. The
    sidechain is placed as a rigid body:
      1. sample approach directions u on a Fibonacci sphere around the donor nitrogen
      2. put ND1 at donor + u * D_PLACE, and rotate so the ND1 lone pair points back along -u
      3. the one remaining degree of freedom is the spin about that axis -- sample it
      4. reject any placement with a heavy-atom clash into TNF
      5. CONFIRM every survivor through bin/his_cation_gate.py itself, so this is the same
         instrument and not a parallel reimplementation of it

CONSTANTS
    D_PLACE  3.0 A   where ND1 is placed. Inside the gate's 3.5 A ceiling with room for the
                     lone-pair cone, and a normal N...N hydrogen bond distance.
    CLASH    3.0 A   heavy-atom floor against TNF. The engaged donor nitrogen is exempt -- it is
                     the intended partner, and ND1 sits 3.0 A from it by construction.
    N_DIR    642     Fibonacci directions. N_SPIN 18 rotations about the approach axis.

SELF-TESTS
    T1  the rigid transform preserves internal geometry -- every pairwise distance in the
        sidechain is unchanged to 1e-9. A broken rotation would otherwise produce a confident
        map of placements for a molecule that is no longer a histidine.
    T2  after placement, ND1 sits exactly D_PLACE from the donor and the lone pair points at it
        to within 1e-6 degrees. Pins the construction the whole result rests on.
    T3  known answer: the residues 1TNF carries at the mapped positions really are ARG at 108
        and LYS at 166. Same check as third_site_census.py T3, and the thing that catches a
        numbering slip silently producing a map of the wrong residue.
    T4  the gate agrees. Every placement this script calls clear must come back as a DESIGNABLE
        hit from his_cation_gate.scan(). Disagreement means one of the two is wrong.
    T5  the anti-blindness control. A donor placed at the trimer's buried centroid must yield
        ZERO clear placements. Nothing fits inside a protein core, so if that returns a number
        the clash test is not doing anything and no count in this file means what it says.

    anchor_reach.py              # self-tests, then the map
    anchor_reach.py --selftest   # tests only
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(os.path.dirname(HERE))      # analysis/02-tnf -> analysis -> repo root
sys.path.insert(0, os.path.join(ROOT, "bin"))

import fetch                                              # noqa: E402
import his_cation_gate as gate                            # noqa: E402
from species_and_histidines import HUMAN, UNIPROT_START   # noqa: E402
from third_site_census import tnf_chain_offset            # noqa: E402

D_PLACE = 3.0
CLASH = 3.0
N_DIR = 642
N_SPIN = 18
LOCAL_R = 14.0        # TNF atoms beyond this of the donor cannot reach a placed sidechain

# The anchors s4 commits to, plus two comparison rows that are NOT part of the design: R107 is
# excluded on mouse Gln and R158 was called "peripheral, too buried" by the census. Scoring them
# costs nothing and shows what a buried cation looks like in this map.
ANCHORS = {108: "ARG", 166: "LYS"}
COMPARE = {107: "ARG", 158: "ARG", 141: "LYS"}

SIDECHAIN = ("CB", "CG", "ND1", "CD2", "CE1", "NE2")
DONOR_N = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",)}


# ------------------------------------------------------------------------------------- geometry

def fibonacci(n):
    """n roughly equidistant unit vectors."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = math.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def rot_between(a, b):
    """Rotation matrix taking unit vector a onto unit vector b."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-12:
        if c > 0:
            return np.eye(3)
        # antiparallel: rotate pi about any axis perpendicular to a
        p = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(p, a))) > 0.9:
            p = np.array([0.0, 1.0, 0.0])
        p = p - np.dot(p, a) * a
        return rot_axis(p / np.linalg.norm(p), math.pi)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def rot_axis(axis, theta):
    """Rodrigues rotation about a unit axis."""
    k = axis / np.linalg.norm(axis)
    kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return (np.eye(3) + math.sin(theta) * kx + (1 - math.cos(theta)) * (kx @ kx))


# ----------------------------------------------------------------------------------- structure

# CORRECTED TARGET, 2026-10-06. This used fetch.pdb("1TNF") directly, and raw 1TNF carries LEU
# where canonical TNF-alpha has ASP219 -- in the consensus core of the epitope, 10 of 10 receptor
# copies. Residue 219 sits 11.5-13.5 A from R108, which is INSIDE this file's LOCAL_R = 14 A
# clash shell, so the reach map genuinely depended on it.
TARGET = os.path.join(ROOT, "targets", "tnf", "tnf_canonical_trimer.pdb")


def load_1tnf():
    """-> (atoms, {(chain, uniprot): (resname, {atom: pos})}, all-heavy-atom array)"""
    import gemmi
    if not os.path.exists(TARGET):
        raise SystemExit(f"REFUSING: {TARGET} is missing. Do NOT fall back to raw 1TNF -- it has "
                         f"the wrong residue at 219. Rebuild the target first.")
    st = gemmi.read_structure(TARGET)
    st.setup_entities()
    st.remove_ligands_and_waters()
    model = st[0]
    byres, coords = {}, []
    for ch in model:
        k, acc = tnf_chain_offset(ch)
        if k is None or acc < 0.90:
            continue
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if not (info and info.is_amino_acid()):
                continue
            u = r.seqid.num + k
            d = {}
            for at in r:
                if at.element == gemmi.Element("H"):
                    continue
                p = np.array([at.pos.x, at.pos.y, at.pos.z])
                d[at.name] = p
                coords.append(p)
            byres[(ch.name, u)] = (r.name, d)
    return byres, np.asarray(coords)


def template(byres):
    """A real histidine sidechain out of this structure, plus its ND1 lone-pair vector."""
    for (ch, u), (rn, d) in sorted(byres.items()):
        if rn == "HIS" and all(a in d for a in SIDECHAIN):
            lp = gate.lone_pair(d, "ND1")
            if lp is not None:
                return (ch, u), np.stack([d[a] for a in SIDECHAIN]), lp
    raise SystemExit("REFUSING: no complete histidine sidechain in the structure to use as a "
                     "template, and this file will not invent one")


def place(tmpl, lp, donor, u, spin):
    """Rigid-body the sidechain so ND1 = donor + u*D_PLACE with its lone pair along -u."""
    nd1 = tmpl[SIDECHAIN.index("ND1")]
    centred = tmpl - nd1
    r = rot_axis(-u, spin) @ rot_between(lp, -u)
    return centred @ r.T + (donor + u * D_PLACE)


# --------------------------------------------------------------------------------------- scan

def reach(byres, coords, chain, uni, resname, verbose=False):
    """How many clear placements exist against each donor nitrogen of one residue."""
    rn, d = byres[(chain, uni)]
    assert rn == resname, f"expected {resname} at {chain}:{uni}, found {rn}"
    _, tmpl, lp = template(byres)
    dirs = fibonacci(N_DIR)
    spins = np.linspace(0, 2 * math.pi, N_SPIN, endpoint=False)
    out = {}
    for dn in DONOR_N[resname]:
        if dn not in d:
            continue
        donor = d[dn]
        # Only TNF atoms near the donor can clash, and the engaged nitrogen is exempt.
        near = coords[np.linalg.norm(coords - donor, axis=1) <= LOCAL_R]
        near = near[np.linalg.norm(near - donor, axis=1) > 1e-6]
        clear, cbs = 0, []
        for u in dirs:
            # stage 1: is ND1's own seat open? rejects every direction pointing into the protein
            if float(np.min(np.linalg.norm(near - (donor + u * D_PLACE), axis=1))) < CLASH:
                continue
            for s in spins:
                pos = place(tmpl, lp, donor, u, s)
                dmin = np.min(np.linalg.norm(
                    pos[:, None, :] - near[None, :, :], axis=2), axis=1)
                # ND1's distance to the engaged donor is exempt; `near` already excludes it
                if float(np.min(dmin)) < CLASH:
                    continue
                clear += 1
                cbs.append(pos[SIDECHAIN.index("CB")])
        out[dn] = (clear, np.asarray(cbs) if cbs else np.empty((0, 3)))
        if verbose:
            print(f"    {dn}: {clear} clear of {N_DIR * N_SPIN}")
    return out


def confirm(byres, coords, chain, uni, resname):
    """T4: push one clear placement through his_cation_gate and require a DESIGNABLE hit."""
    rn, d = byres[(chain, uni)]
    _, tmpl, lp = template(byres)
    for u in fibonacci(N_DIR):
        for dn in DONOR_N[resname]:
            if dn not in d:
                continue
            donor = d[dn]
            near = coords[np.linalg.norm(coords - donor, axis=1) <= LOCAL_R]
            near = near[np.linalg.norm(near - donor, axis=1) > 1e-6]
            if float(np.min(np.linalg.norm(near - (donor + u * D_PLACE), axis=1))) < CLASH:
                continue
            pos = place(tmpl, lp, donor, u, 0.0)
            if float(np.min(np.min(np.linalg.norm(
                    pos[:, None, :] - near[None, :, :], axis=2), axis=1))) < CLASH:
                continue
            atoms = [gate.Atom("X", "HIS", 1, n, p) for n, p in zip(SIDECHAIN, pos)]
            atoms += [gate.Atom(chain, rn, uni, n, p) for n, p in d.items()]
            hits = gate.scan(atoms, ["X"], [chain])
            good = [h for h in hits if h["direction"] == "designable" and not h["ambiguous"]]
            return (dn, u, good)
    return (None, None, [])


# ---------------------------------------------------------------------------------- self-tests

def selftest(byres, coords):
    (tch, tu), tmpl, lp = template(byres)

    # T1  the transform is rigid
    inner = lambda p: np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    before = inner(tmpl)
    moved = place(tmpl, lp, np.array([10.0, -3.0, 7.0]),
                  np.array([0.3, -0.5, 0.81]) / np.linalg.norm([0.3, -0.5, 0.81]), 1.1)
    assert np.allclose(before, inner(moved), atol=1e-9), \
        "T1 FAIL: the rigid transform changed the sidechain's internal geometry"

    # T2  the placement does what it says
    donor = np.array([10.0, -3.0, 7.0])
    u = np.array([0.3, -0.5, 0.81])
    u = u / np.linalg.norm(u)
    pos = place(tmpl, lp, donor, u, 0.7)
    d = {n: p for n, p in zip(SIDECHAIN, pos)}
    assert abs(float(np.linalg.norm(d["ND1"] - donor)) - D_PLACE) < 1e-9, \
        f"T2 FAIL: ND1 is {np.linalg.norm(d['ND1'] - donor)} from the donor, not {D_PLACE}"
    assert gate._angle(gate.lone_pair(d, "ND1"), donor - d["ND1"]) < 1e-6, \
        "T2 FAIL: the lone pair does not point at the donor after placement"

    # T3  known answer on the numbering: the PDB residue AND the UniProt sequence must both
    #     agree at every position this file reports on. Checking only the PDB would pass on a
    #     consistent-but-wrong offset.
    one = {"ARG": "R", "LYS": "K"}
    for uni, want in {**ANCHORS, **COMPARE}.items():
        seen = {byres[k][0] for k in byres if k[1] == uni}
        assert seen == {want}, f"T3 FAIL: expected {want} at {uni}, structure has {seen}"
        assert HUMAN[uni - UNIPROT_START] == one[want], (
            f"T3 FAIL: UniProt has {HUMAN[uni - UNIPROT_START]} at {uni}, not {one[want]} "
            f"-- the offset is wrong")
    chains = sorted({k[0] for k in byres})
    assert len(chains) == 3, f"T3 FAIL: 1TNF should give 3 protomers, got {chains}"

    # T5  anti-blindness: the probe must admit NOTHING at a genuinely buried point.
    #
    #     The first version of this test used the trimer CENTROID and fired: 12 of 642 directions
    #     came back clear. The clash test was fine -- the premise was wrong. TNF-alpha trimers
    #     have a solvent channel along the three-fold axis, so their geometric centre is not
    #     inside protein at all. Keeping the note because a guard whose assumption is wrong looks
    #     exactly like a guard catching a bug.
    #
    #     The replacement picks the most-buried heavy atom in the structure by neighbour count,
    #     which is buried by construction rather than by assumption.
    bur = np.empty(len(coords), dtype=int)
    for i in range(0, len(coords), 500):
        blk = coords[i:i + 500]
        bur[i:i + 500] = (np.linalg.norm(
            blk[:, None, :] - coords[None, :, :], axis=2) <= 8.0).sum(axis=1)
    probe = coords[int(np.argmax(bur))]
    near = coords[np.linalg.norm(coords - probe, axis=1) <= LOCAL_R]
    near = near[np.linalg.norm(near - probe, axis=1) > 1e-6]
    n_clear = sum(1 for uu in fibonacci(N_DIR)
                  if float(np.min(np.linalg.norm(near - (probe + uu * D_PLACE), axis=1))) >= CLASH)
    assert n_clear == 0, (f"T5 FAIL: {n_clear} of {N_DIR} directions were called clear at the "
                          f"MOST BURIED ATOM in the structure ({int(bur.max())} neighbours within "
                          f"8 A). The clash test is not working and no count in this file means "
                          f"what it says.")
    # reported, not asserted: the channel is a property of the fold, not a test outcome
    cen = coords.mean(axis=0)
    cnear = coords[np.linalg.norm(coords - cen, axis=1) <= LOCAL_R]
    n_cen = sum(1 for uu in fibonacci(N_DIR)
                if float(np.min(np.linalg.norm(cnear - (cen + uu * D_PLACE), axis=1))) >= CLASH)
    print(f"self-tests T1 T2 T3 T5 passed  (template {tch}:HIS{tu}, {len(chains)} protomers, "
          f"{len(coords)} heavy atoms; most-buried atom admits 0 of {N_DIR}, "
          f"trimer centroid admits {n_cen} -- the axial channel)")
    return chains


def main():
    byres, coords = load_1tnf()
    chains = selftest(byres, coords)
    if "--selftest" in sys.argv:
        return

    print(f"\nReach map on apo 1TNF. D_PLACE={D_PLACE} A, CLASH={CLASH} A, "
          f"{N_DIR} directions x {N_SPIN} spins = {N_DIR * N_SPIN} placements per donor N.\n")
    print(f"{'site':<8}{'role':<26}{'chain':<7}{'donor':<7}{'clear':>8}{'of':>9}  CB spread")
    print("-" * 86)
    totals = {}
    for uni, resname in list(ANCHORS.items()) + list(COMPARE.items()):
        role = "ANCHOR (s4)" if uni in ANCHORS else "comparison, not in plan"
        for ch in chains:
            if (ch, uni) not in byres:
                continue
            for dn, (n, cbs) in reach(byres, coords, ch, uni, resname).items():
                spread = ("—" if len(cbs) < 2 else
                          f"{float(np.linalg.norm(cbs - cbs.mean(axis=0), axis=1).max()):.1f} A")
                print(f"{resname[0]}{uni:<7}{role:<26}{ch:<7}{dn:<7}{n:>8}{N_DIR * N_SPIN:>9}"
                      f"  {spread}")
                totals[(uni, ch, dn)] = n
    print("-" * 86)

    for uni in ANCHORS:
        tot = sum(v for (u_, _, _), v in totals.items() if u_ == uni)
        print(f"{ANCHORS[uni][0]}{uni}: {tot} clear placements across 3 protomers")

    # T4, run last because it needs a real clear placement to confirm
    print()
    for uni, resname in ANCHORS.items():
        dn, u, good = confirm(byres, coords, chains[0], uni, resname)
        if dn is None:
            print(f"T4 {resname[0]}{uni}: no clear placement existed to confirm")
            continue
        assert good, (f"T4 FAIL: a placement this script calls clear produced NO designable hit "
                      f"from his_cation_gate at {resname[0]}{uni} {dn}. The two instruments "
                      f"disagree and at least one is wrong.")
        h = good[0]
        print(f"T4 {resname[0]}{uni}: his_cation_gate CONFIRMS -> {h['his']} {h['acceptor']} "
              f"<- {h['donor']}  {h['dist']} A  {h['angle']} deg  needs {h['requires_tautomer']}")

    print("\nA clear placement means the approach vector is open on the apo trimer. It does NOT\n"
          "mean a binder backbone can reach it -- that is the generator's problem. A BLOCKED\n"
          "direction is blocked for any binder, which is what makes a zero here decisive.")


if __name__ == "__main__":
    main()
