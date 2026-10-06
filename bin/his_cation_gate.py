#!/usr/bin/env python3
"""
his_cation_gate -- Ahn et al.'s SECOND selection criterion, built as a GATE on dddG_elec.

WHY THIS IS A GATE AND NOT A SIBLING FILTER
    Ahn et al. (bioRxiv 2025.09.29.678932) selected designs on TWO criteria, verbatim:
        `ddg elec >= 0` AND "at least one histidine accepting a hydrogen bond from a
        positively charged residue".
    bin/dddg_elec.py implements the first. On 2026-10-05 it failed its known-answer control
    0 of 5 (challenges/02-tnf-alpha.md s6b): all three MEASURED adalimumab pH-switch variants
    of Schroeter et al. (231x, 505x, 785x) scored NEGATIVE, and a scramble negative with five
    histidines placed >12 A away scored -0.129 against wild type's -0.125 -- i.e. the score was
    insensitive to where the histidines went. Every number was dominated by one pre-existing
    residue, H:HIS57, sitting 6.2 A from a GLU.

    That is the diagnosis: dddG_elec computed across an interface with no cationic partner is
    not a weak signal, it is BACKGROUND. So this criterion is not the second half of a pair.
    It is the PRECONDITION for the first one meaning anything, and it runs first.

    playbook s24 at the design level: a guard that refuses to score beats one that scores blind.

WHY IT ENCODES THE MECHANISM (which is why it belongs in the generator too, lesson 12)
    At pH 7.4 a neutral histidine has a lone pair on exactly one ring nitrogen, which ACCEPTS
    an H-bond from an Arg/Lys donor -> the complex is stabilised. Protonate it and two things
    happen at once, both pushing the same way:
        1. that nitrogen now carries a hydrogen, so the H-bond is LOST, and
        2. a +1 imidazolium now sits against a +1 guanidinium/ammonium -> repulsion.
    Bind at neutral, release at acid. That is problem 2's direction. A histidine with no
    cationic partner across the interface has neither effect available to it.

HEAVY ATOMS ONLY, AND WHY
    The usual H-bond criterion uses the D-H...A angle. Crystal structures have no hydrogens,
    and when they are added the positions are whatever the adder chose -- tonight's G2 guard in
    dddg_elec.py fired precisely because OpenMM's Modeller.addHydrogens re-optimised the whole
    H-bond network between two calls and moved hydroxyl hydrogens that were supposed to cancel.
    An angle computed on inferred hydrogens is therefore not reproducible across tools.

    The acceptor-side geometry does NOT need them. The lone pair of an imidazole nitrogen lies
    in the ring plane and bisects the EXTERIOR of its C-N-C angle, so it is fixed by three heavy
    atoms. This file uses:
        distance   donor N  ->  His ND1/NE2, heavy atom to heavy atom
        direction  angle between the acceptor's lone-pair vector and the vector to the donor
    Both are computable from the deposited coordinates alone.

THRESHOLDS -- written down before any structure was scored (playbook s13)
    D_HB    3.5 A   heavy-atom N...N. The conventional protein H-bond donor-acceptor ceiling;
                    reported as a band 3.2 / 3.5 / 3.9 A rather than asserted as the truth.
    CONE    60 deg  OURS, not published. Labelled, and banded 45 / 60 / 90 deg. 90 deg is the
                    weakest defensible form (donor merely on the lone-pair side of the N).
    CUT     4.5 A   heavy-atom contact distance defining "at the interface". Same value as
                    analysis/02-tnf/per_partner.py, so the two agree on what an interface is.

    If the verdict changes across that 3x3 band, the verdict is NOT ESTABLISHED and this file
    says so instead of picking the cell that reads well (playbook s16).

TWO DIRECTIONS, REPORTED SEPARATELY, BECAUSE ONLY ONE OF THEM IS OURS TO DESIGN
    binder-His <- target-cation    the designable one. R108 and K166 on TNF are the anchors
                                   s4 names; the histidine is ours to place.
    target-His <- binder-cation    a different mechanism. It switches on TNF's OWN histidine,
                                   whose pKa and position we do not choose. Reported, never
                                   counted toward the primary verdict.
    Collapsing these would let a design pass on a mechanism it cannot control.

HIS AS A DONOR IS AMBIGUOUS ON PURPOSE
    A protonated His can donate to a neutral His, but then BOTH sides are pH-dependent and the
    pair has no defined direction of switching. Arg and Lys are the primary donor set; His is
    scored into a separate `ambiguous` column and excluded from the verdict. Backbone and
    terminal amines are excluded too -- a terminus is not a designable sidechain.

GUARDS (each one has a mutation in --mutate that must turn it red; playbook s25)
    G1  no histidine at the interface -> REFUSE, with "no histidine present", never "0 hits".
        playbook s10: a zero is not an absence. dddG_elec's own null passes `>= 0` trivially.
    G2  cross-interface only. An intra-chain His-Arg pair is not a binding switch.
    G3  the band is always printed, and a verdict that flips inside it is reported as unset.
    G4  tautomer feasibility. A neutral His has ONE acceptor nitrogen. A hit on ND1 requires
        the HIE tautomer, a hit on NE2 requires HID. A histidine whose BOTH nitrogens are
        engaged by cationic donors has no neutral tautomer that satisfies both -> flagged
        `tautomer_conflict` and not counted.
    G5  the anti-blindness test, which is the one that matters. A filter that returns zero
        everywhere is indistinguishable from a filter that is broken. T6 therefore re-runs
        3WD5 under a deliberately slack criterion (any His nitrogen within 8 A of any cationic
        nitrogen, no angle at all) and REQUIRES a non-zero count. If the slack form also
        returns zero, the zero is a parsing bug and the file fails rather than reporting a
        finding.

KNOWN ANSWER, ALREADY IN HAND
    3WD5 is adalimumab Fab + monomeric TNF. analysis/02-tnf/per_partner.json gives its epitope
    as TNF 96 97 99 141-143 186 187 189 191 216-223: no R108, no K166, no H149, and K141 was
    classed peripheral by the day-1 census. So this gate must return ZERO designable hits on
    3WD5 -- which converts s6b's after-the-fact epitope explanation into an instrument output
    produced BEFORE scoring. That is T5.

USAGE
    his_cation_gate.py --selftest                     # T1-T6, no network, no I/O
    his_cation_gate.py --mutate cone|cross|dist       # each must make a named test go red
    his_cation_gate.py --pdb 3WD5 --binder H,L --target A
    his_cation_gate.py --pool                         # pass rate over the 13 TNF complexes (s14)
"""

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "analysis", "02-tnf"))

# ---------------------------------------------------------------- published / declared constants

D_HB = 3.5                      # A, heavy-atom donor N ... acceptor N
D_BAND = (3.2, 3.5, 3.9)        # reported band, playbook s16
CONE = 60.0                     # deg, OURS -- the paper publishes no angle
CONE_BAND = (45.0, 60.0, 90.0)
CUT = 4.5                       # A, "at the interface"; same as analysis/02-tnf/per_partner.py

SLACK_D = 8.0                   # G5 only: the deliberately blind criterion

# Cationic sidechain donor nitrogens. CZ/CE are carbons and are not donors.
DONORS = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",)}
AMBIGUOUS_DONORS = {"HIS": ("ND1", "NE2")}

# Imidazole ring connectivity: each acceptor nitrogen and its two ring neighbours. The lone
# pair bisects the exterior of that C-N-C angle.
RING = {"ND1": ("CG", "CE1"), "NE2": ("CD2", "CE1")}

# The neutral tautomer each hit requires: a hit on ND1 needs the proton on NE2, i.e. HIE.
REQUIRES = {"ND1": "HIE", "NE2": "HID"}

# Human TNF-alpha mature sequence, for identifying the target side by CONTENT rather than by
# chain letter -- 29 of 69 problem-1 run directories had the target in chain A. Copied from
# analysis/02-tnf/per_partner.py; that file is a script with top-level side effects and cannot
# be imported.
HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTH"
         "TISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")


class Atom:
    """One heavy atom. Deliberately not a gemmi object: the geometry below is then testable
    with hand-built coordinates and no structure file at all."""
    __slots__ = ("chain", "resname", "resnum", "name", "pos")

    def __init__(self, chain, resname, resnum, name, pos):
        self.chain, self.resname, self.resnum = chain, resname, resnum
        self.name, self.pos = name, np.asarray(pos, dtype=float)

    @property
    def key(self):
        return (self.chain, self.resnum)

    def __repr__(self):
        return f"{self.chain}:{self.resname}{self.resnum}:{self.name}"


# ---------------------------------------------------------------------------------- geometry

def _unit(v):
    n = float(np.linalg.norm(v))
    return None if n < 1e-9 else v / n


def _angle(u, v):
    uu, vv = _unit(u), _unit(v)
    if uu is None or vv is None:
        return None
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(uu, vv))))))


def lone_pair(res, n_name, flip=False):
    """In-plane unit vector along the lone pair of His ND1/NE2, from heavy atoms only.
    `flip` is mutation M1: it points the vector INTO the ring, which must break T3."""
    if n_name not in res:
        return None
    a, b = RING[n_name]
    if a not in res or b not in res:
        return None                      # incomplete sidechain: no claim, rather than a guess
    v = res[n_name] - 0.5 * (res[a] + res[b])
    v = _unit(v)
    if v is None:
        return None
    return -v if flip else v


# ------------------------------------------------------------------------------- the criterion

def residues(atoms):
    """{(chain, resnum): (resname, {atomname: pos})}"""
    out = {}
    for a in atoms:
        name, d = out.setdefault(a.key, (a.resname, {}))
        d[a.name] = a.pos
    return out


def interface_his(atoms, side, other, cut=CUT):
    """Histidines on `side` with a sidechain heavy atom within `cut` of any heavy atom on
    `other`. Sidechain only: a histidine touching across the interface through its backbone
    is not positioned to do chemistry with its ring."""
    sc = {"CB", "CG", "ND1", "CD2", "CE1", "NE2"}
    opp = np.array([a.pos for a in atoms if a.chain in other]) if other else np.empty((0, 3))
    hits = set()
    if len(opp) == 0:
        return hits
    for a in atoms:
        if a.chain in side and a.resname == "HIS" and a.name in sc:
            if float(np.min(np.linalg.norm(opp - a.pos, axis=1))) <= cut:
                hits.add(a.key)
    return hits


def scan(atoms, binder, target, d_hb=D_HB, cone=CONE,
         flip=False, allow_intra=False, slack=False):
    """Every cationic-donor -> His-acceptor pair across the interface.

    Returns a list of dicts. `direction` is 'designable' when the histidine is on the binder
    (ours to place) and 'target-his' when it is on TNF (not ours to choose).
    """
    res = residues(atoms)
    sides = {}
    for ch in binder:
        sides[ch] = "binder"
    for ch in target:
        sides[ch] = "target"

    acceptors = []      # (key, resname, n_name, pos, lp)
    donors = []         # (key, resname, atom_name, pos, ambiguous)
    for key, (rn, d) in res.items():
        if key[0] not in sides:
            continue
        if rn == "HIS":
            for n_name in ("ND1", "NE2"):
                if n_name in d:
                    lp = lone_pair(d, n_name, flip=flip)
                    if lp is not None:
                        acceptors.append((key, rn, n_name, d[n_name], lp))
        for table, amb in ((DONORS, False), (AMBIGUOUS_DONORS, True)):
            for name in table.get(rn, ()):
                if name in d:
                    donors.append((key, rn, name, d[name], amb))

    hits = []
    for akey, arn, n_name, apos, lp in acceptors:
        for dkey, drn, dname, dpos, amb in donors:
            if akey == dkey:
                continue
            same = sides[akey[0]] == sides[dkey[0]]
            if same and not allow_intra:
                continue                          # G2
            dist = float(np.linalg.norm(dpos - apos))
            limit = SLACK_D if slack else d_hb
            if dist > limit:
                continue
            ang = _angle(lp, dpos - apos)
            if not slack and (ang is None or ang > cone):
                continue
            hits.append({
                "his": f"{akey[0]}:HIS{akey[1]}", "his_key": akey, "acceptor": n_name,
                "donor": f"{dkey[0]}:{drn}{dkey[1]}:{dname}", "donor_res": drn,
                "dist": round(dist, 2), "angle": None if ang is None else round(ang, 1),
                "requires_tautomer": REQUIRES[n_name],
                "ambiguous": amb,
                "intra": same,
                "direction": "designable" if sides[akey[0]] == "binder" else "target-his",
            })
    return hits


def tautomer_conflicts(hits):
    """G4: a histidine with cationic donors on BOTH ring nitrogens. No neutral tautomer
    satisfies both, so neither hit is a clean switch."""
    seen = {}
    for h in hits:
        if h["ambiguous"]:
            continue
        seen.setdefault(h["his_key"], set()).add(h["acceptor"])
    return {k for k, v in seen.items() if len(v) > 1}


def verdict(atoms, binder, target, d_hb=D_HB, cone=CONE, **kw):
    """The gate itself. Returns (state, payload).

    state is 'REFUSE' | 'PASS' | 'FAIL'. REFUSE is not FAIL: it means the question was not
    asked, because there is no histidine at the interface to ask it about (G1, playbook s10).
    """
    bi = interface_his(atoms, set(binder), set(target))
    ti = interface_his(atoms, set(target), set(binder))
    if not bi and not ti:
        return "REFUSE", {"reason": "no histidine at the interface on either side",
                          "n_his_binder": 0, "n_his_target": 0, "hits": []}

    hits = scan(atoms, binder, target, d_hb=d_hb, cone=cone, **kw)
    conflict = tautomer_conflicts(hits)
    good = [h for h in hits if not h["ambiguous"] and h["his_key"] not in conflict]
    designable = [h for h in good if h["direction"] == "designable"]
    return ("PASS" if designable else "FAIL"), {
        "n_his_binder": len(bi), "n_his_target": len(ti),
        "hits": hits, "designable": designable,
        "target_his": [h for h in good if h["direction"] == "target-his"],
        "ambiguous": [h for h in hits if h["ambiguous"]],
        "tautomer_conflict": sorted(conflict),
    }


def band(atoms, binder, target):
    """G3: the verdict at every cell of the declared 3x3 band."""
    rows = []
    for d in D_BAND:
        for c in CONE_BAND:
            st, pay = verdict(atoms, binder, target, d_hb=d, cone=c)
            rows.append((d, c, st, len(pay.get("designable", []))))
    states = {st for _, _, st, _ in rows}
    return rows, ("ESTABLISHED" if len(states) == 1 else "NOT ESTABLISHED")


# --------------------------------------------------------------------------------- structure I/O

def load(path):
    """Heavy atoms out of a CIF/PDB, plus which chains look like TNF-alpha."""
    import gemmi
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_ligands_and_waters()
    model = st[0]
    atoms, chains = [], {}
    for ch in model:
        seq, n = [], 0
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if not (info and info.is_amino_acid()):
                continue
            n += 1
            seq.append((r.seqid.num, info.one_letter_code.upper()))
            for at in r:
                if at.element == gemmi.Element("H"):
                    continue
                atoms.append(Atom(ch.name, r.name, r.seqid.num, at.name,
                                  (at.pos.x, at.pos.y, at.pos.z)))
        chains[ch.name] = (n, seq)
    return atoms, chains


def tnf_chains(chains):
    """Chains that ARE TNF, by sequence identity, not by letter. Mirrors per_partner.py's
    best_offset: slide the numbering, require >=40 aligned and >=90% identity."""
    out = set()
    for name, (n, seq) in chains.items():
        best = 0.0
        for k in range(-20, 130):
            m = t = 0
            for num, c in seq:
                u = num + k - 77
                if 0 <= u < len(HUMAN):
                    t += 1
                    m += (c == HUMAN[u])
            if t >= 40 and m / t > best:
                best = m / t
        if best >= 0.90:
            out.add(name)
    return out


# ------------------------------------------------------------------------------------ reporting

def report(label, atoms, binder, target):
    st, pay = verdict(atoms, binder, target)
    rows, established = band(atoms, binder, target)
    print(f"\n### {label}   binder={','.join(sorted(binder))}  target={','.join(sorted(target))}")
    if st == "REFUSE":
        print(f"  REFUSE -- {pay['reason']}")
        print("  (this is not a FAIL. dddG_elec's `>= 0` would have passed this pose trivially.)")
        return st, pay
    print(f"  interface histidines: binder {pay['n_his_binder']}, target {pay['n_his_target']}")
    print(f"  VERDICT {st} at D_HB={D_HB} A, cone={CONE} deg   [band: {established}]")
    for tag, key in (("DESIGNABLE  binder-His <- target-cation", "designable"),
                     ("target-his  target-His <- binder-cation", "target_his"),
                     ("ambiguous   His as donor, excluded", "ambiguous")):
        hs = pay.get(key) or []
        print(f"    {tag}: {len(hs)}")
        for h in hs:
            print(f"      {h['his']} {h['acceptor']} <- {h['donor']}  "
                  f"{h['dist']} A  {h['angle']} deg  needs {h['requires_tautomer']}")
    if pay["tautomer_conflict"]:
        print(f"    tautomer_conflict (both N engaged, no neutral form works): "
              f"{pay['tautomer_conflict']}")
    cells = " ".join(f"{d}/{int(c)}:{s[0]}{n}" for d, c, s, n in rows)
    print(f"    band  {cells}")
    if established == "NOT ESTABLISHED":
        print("    the verdict FLIPS inside the declared band -- it is not established (s16)")
    return st, pay


# ----------------------------------------------------------------------------------- self-tests

def _his(chain, num, lp_dir=(1.0, 0.0, 0.0)):
    """A minimal imidazole whose ND1 lone pair points along `lp_dir` from the origin.

    ND1 at the origin; CG and CE1 placed so their midpoint lies at -lp_dir, which is exactly
    the construction lone_pair() inverts. CB/CD2/NE2 are present so interface_his() sees a
    sidechain, and sit behind the ring where they cannot be mistaken for the lone-pair side.
    """
    u = np.asarray(lp_dir, dtype=float)
    u = u / np.linalg.norm(u)
    perp = np.cross(u, [0.0, 0.0, 1.0])
    if np.linalg.norm(perp) < 1e-6:
        perp = np.cross(u, [0.0, 1.0, 0.0])
    perp = perp / np.linalg.norm(perp)
    return [
        Atom(chain, "HIS", num, "ND1", (0.0, 0.0, 0.0)),
        Atom(chain, "HIS", num, "CG", -u * 1.38 + perp * 1.1),
        Atom(chain, "HIS", num, "CE1", -u * 1.38 - perp * 1.1),
        Atom(chain, "HIS", num, "CD2", -u * 2.6 + perp * 2.0),
        Atom(chain, "HIS", num, "NE2", -u * 2.9 - perp * 1.6),
        Atom(chain, "HIS", num, "CB", -u * 2.9 + perp * 2.6),
    ]


def _arg(chain, num, pos):
    p = np.asarray(pos, dtype=float)
    return [Atom(chain, "ARG", num, "NH1", p),
            Atom(chain, "ARG", num, "CZ", p + np.array([0.0, 1.3, 0.0])),
            Atom(chain, "ARG", num, "CB", p + np.array([0.0, 3.0, 0.0]))]


def selftest(flip=False, allow_intra=False, d_override=None):
    """T1-T6. Mutations are injected through the three keyword arguments so that each one
    maps to exactly one named assertion below."""
    d_hb = D_HB if d_override is None else d_override
    kw = dict(flip=flip, allow_intra=allow_intra)
    B, T = ["B"], ["T"]

    # T0  THE FIXTURE ITSELF. T1-T4 are only meaningful if _his really does point ND1's lone
    #     pair along the requested direction -- otherwise they could all pass on coincidental
    #     geometry and the suite would be the "guard that passes while blind" again (s24).
    for want in ((1, 0, 0), (0, 1, 0), (-0.3, 0.5, 0.81)):
        res = residues(_his("B", 1, want))[("B", 1)][1]
        got = lone_pair(res, "ND1")
        u = np.asarray(want, dtype=float)
        u = u / np.linalg.norm(u)
        off = _angle(got, u)
        assert off is not None and off < 1.0, f"T0 FAIL: fixture lone pair off by {off} deg"
        assert _angle(lone_pair(res, "ND1", flip=True), u) > 179.0, "T0 FAIL: flip is not a flip"

    # T1  donor 3.0 A out along the lone pair, on the far side of the interface -> detected
    a = _his("B", 10, (1, 0, 0)) + _arg("T", 20, (3.0, 0.0, 0.0))
    h = scan(a, B, T, d_hb=d_hb, cone=CONE, **kw)
    assert len(h) == 1 and h[0]["direction"] == "designable", f"T1 FAIL: {h}"
    assert h[0]["requires_tautomer"] == "HIE", f"T1 FAIL tautomer: {h}"

    # T2  same direction, 5.0 A -> rejected on DISTANCE. Mutation M3 (--mutate dist) breaks this.
    a = _his("B", 10, (1, 0, 0)) + _arg("T", 20, (5.0, 0.0, 0.0))
    h = scan(a, B, T, d_hb=d_hb, cone=CONE, **kw)
    assert not h, f"T2 FAIL: a donor 5.0 A away was accepted at d_hb={d_hb}: {h}"

    # T3  3.0 A but on the OPPOSITE side of the ring -> rejected on ANGLE.
    #     Mutation M1 (--mutate cone) inverts the lone-pair vector and must break this.
    a = _his("B", 10, (1, 0, 0)) + _arg("T", 20, (-3.0, 0.0, 0.0))
    h = scan(a, B, T, d_hb=d_hb, cone=CONE, **kw)
    assert not h, f"T3 FAIL: a donor behind the ring was accepted: {h}"

    # T4  good geometry but SAME side -> not an interface switch.
    #     Mutation M2 (--mutate cross) must break this.
    a = _his("B", 10, (1, 0, 0)) + _arg("B", 20, (3.0, 0.0, 0.0))
    h = scan(a, B, T, d_hb=d_hb, cone=CONE, **kw)
    assert not h, f"T4 FAIL: an intra-chain pair was counted: {h}"

    # T4b G1: no histidine anywhere -> REFUSE, distinct from FAIL
    a = _arg("B", 20, (0.0, 0.0, 0.0)) + _arg("T", 21, (3.0, 0.0, 0.0))
    st, pay = verdict(a, B, T)
    assert st == "REFUSE", f"T4b FAIL: expected REFUSE with no histidine, got {st}"

    # T4c G4: donors on BOTH ring nitrogens -> tautomer_conflict, and the verdict is FAIL,
    #     not PASS. The second donor is placed along NE2's own lone pair so it is a real hit
    #     rather than a near miss that would make this test vacuous.
    his = _his("B", 10, (1, 0, 0))
    ne2_lp = lone_pair(residues(his)[("B", 10)][1], "NE2")
    ne2_pos = next(x.pos for x in his if x.name == "NE2")
    a = his + _arg("T", 20, (3.0, 0.0, 0.0)) + _arg("T", 21, ne2_pos + ne2_lp * 3.0)
    h = scan(a, B, T, d_hb=D_HB, cone=CONE)
    assert {x["acceptor"] for x in h} == {"ND1", "NE2"}, \
        f"T4c FAIL: the fixture did not engage both nitrogens, so the test is vacuous: {h}"
    assert tautomer_conflicts(h), f"T4c FAIL: both-N engagement was not flagged: {h}"
    st, _ = verdict(a, B, T)
    assert st == "FAIL", f"T4c FAIL: a tautomer-conflicted His gave {st}, expected FAIL"

    if flip or allow_intra or d_override is not None:
        return                      # mutations are only expected to reach T1-T4

    # ------------------------------------------------------------------ structural known answers
    try:
        import fetch
        path = fetch.cif("3WD5")
    except Exception as e:                                   # noqa: BLE001
        print(f"T5/T6 SKIPPED: 3WD5 unavailable ({e})")
        print("self-tests T1-T4c passed")
        return
    atoms, chains = load(path)
    tnf = tnf_chains(chains)
    fab = {c for c, (n, _) in chains.items() if c not in tnf and n >= 15}
    assert tnf and fab, f"T5 FAIL: chain split went wrong (tnf={tnf}, fab={fab})"

    # T5  the known answer from per_partner.json: adalimumab's epitope holds no R108, no K166,
    #     no H149, and K141 is peripheral -> ZERO designable hits.
    st, pay = verdict(atoms, fab, tnf)
    assert st == "FAIL", f"T5 FAIL: 3WD5 gave {st}, expected FAIL (no cation in this epitope)"
    assert not pay["designable"], f"T5 FAIL: unexpected designable hits on 3WD5: {pay['designable']}"

    # T6  G5, THE anti-blindness test. A filter returning zero everywhere looks exactly like a
    #     broken one. Under the slack criterion -- any His N within 8 A of any cationic N, no
    #     angle at all -- 3WD5 MUST produce hits. If it does not, the zero above is a parsing
    #     bug and this file has no business reporting it as a finding.
    slack = scan(atoms, fab, tnf, slack=True)
    assert slack, ("T6 FAIL: the slack 8 A criterion also found nothing on 3WD5. The zero in T5 "
                   "is therefore not evidence about the epitope -- the scan is blind. Do not "
                   "report T5 until this passes.")
    print(f"self-tests T1-T6 passed  (3WD5: tnf={sorted(tnf)} fab={sorted(fab)}, "
          f"{pay['n_his_binder']} binder + {pay['n_his_target']} target interface His, "
          f"0 designable at 3.5 A/60 deg, {len(slack)} under the slack 8 A control)")


POOL = ["3ALQ", "8ZUI", "7KPB", "3WD5", "4G3Y", "5WUX", "5YOY", "5M2I", "5M2J", "5M2M",
        "9BN7", "9DJW"]


def pool():
    """playbook s14: a filter that passes 0% or 100% of a real pool is not a filter. These are
    the TNF complexes the day-1 census already used, so the pass rate costs nothing to get."""
    import fetch
    tally = {"PASS": 0, "FAIL": 0, "REFUSE": 0, "SKIP": 0}
    for pid in POOL:
        try:
            atoms, chains = load(fetch.cif(pid))
            tnf = tnf_chains(chains)
            partner = {c for c, (n, _) in chains.items() if c not in tnf and n >= 15}
            if not tnf or not partner:
                print(f"\n### {pid}  SKIP  tnf={sorted(tnf)} partner={sorted(partner)}")
                tally["SKIP"] += 1
                continue
            st, _ = report(pid, atoms, partner, tnf)
            tally[st] += 1
        except Exception as e:                               # noqa: BLE001
            print(f"\n### {pid}  SKIP  {type(e).__name__}: {e}")
            tally["SKIP"] += 1
    n = tally["PASS"] + tally["FAIL"]
    print("\n" + "=" * 78)
    print(f"POOL  {tally['PASS']} PASS / {tally['FAIL']} FAIL / {tally['REFUSE']} REFUSE "
          f"/ {tally['SKIP']} SKIP")
    if n:
        rate = tally["PASS"] / n
        print(f"pass rate {rate:.0%} of {n} scored complexes")
        if rate in (0.0, 1.0):
            print("playbook s14: at 0% or 100% this is not yet shown to be a FILTER. Report it "
                  "as a mechanism-availability measurement on TNF-alpha interfaces, which is "
                  "what it is, and do not quote it as a selection yield.")
    print("REFUSE is not FAIL: it counts poses the criterion could not be asked about.")


def main():
    av = sys.argv[1:]
    mut = None
    if "--mutate" in av:
        mut = av[av.index("--mutate") + 1]
    if "--selftest" in av or mut:
        kw = {}
        expect = None
        # which assertion each mutation is expected to reach. `cone` names two because
        # inverting the lone pair breaks T1 (the real hit stops being found) before it reaches
        # T3 (the false hit starts being found) -- either is a correct detection.
        if mut == "cone":
            kw, expect = dict(flip=True), "T1 or T3"
        elif mut == "cross":
            kw, expect = dict(allow_intra=True), "T4"
        elif mut == "dist":
            kw, expect = dict(d_override=99.0), "T2"
        elif mut:
            sys.exit(f"unknown mutation {mut!r}; one of cone, cross, dist")
        try:
            selftest(**kw)
        except AssertionError as e:
            if mut:
                print(f"MUTATION {mut}: correctly RED -> {e}")
                return
            raise
        if mut:
            sys.exit(f"MUTATION {mut} PASSED THE SELF-TESTS. {expect} does not actually test "
                     f"what it claims to -- playbook s25. Fix the test, not the mutation.")
        return
    if "--pool" in av:
        pool()
        return
    if "--pdb" in av:
        import fetch
        pid = av[av.index("--pdb") + 1]
        atoms, chains = load(fetch.cif(pid))
        if "--binder" in av and "--target" in av:
            b = set(av[av.index("--binder") + 1].split(","))
            t = set(av[av.index("--target") + 1].split(","))
        else:
            t = tnf_chains(chains)
            b = {c for c, (n, _) in chains.items() if c not in t and n >= 15}
            print(f"chains assigned by sequence: target(TNF)={sorted(t)} binder={sorted(b)}")
        report(pid, atoms, b, t)
        return
    print(__doc__.strip().rsplit("USAGE", 1)[-1].strip())


if __name__ == "__main__":
    main()
