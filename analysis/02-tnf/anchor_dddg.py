#!/usr/bin/env python3
"""
The control dddG_elec never had: a MECHANISM-PRESENT pose on real coordinates.

WHY THIS EXISTS
    bin/dddg_elec.py failed its known-answer control 0 of 5 (s6b). The diagnosis was that the
    adalimumab epitope contains no cation, so the filter was never shown its own mechanism --
    every number was dominated by one pre-existing histidine 6.2 A from a GLU, and a scramble
    negative with five histidines >12 A away was indistinguishable from wild type.

    s6d then showed the mechanism IS reachable: anchor_reach.py found 5,565 clear placements at
    R108 and 5,976 at K166, and his_cation_gate confirmed two of them as designable hits. So the
    fixture for the missing control already existed. This builds it.

    THE QUESTION IS NOT whether the formula works -- dddg_elec selftest T6 already gives +0.201
    for a His+ 3.5 A from a cation in vacuum. It is whether that anchor signal SURVIVES TNF's
    real electrostatic environment, which is exactly what killed the adalimumab case.

    Three outcomes, all decisive:
      strongly positive    the instrument was never broken, it was shown the wrong case
      near zero / negative the background swamps the anchor at OUR anchors too, so dddG_elec can
                           never rank our designs and must be retired by name with a date (s27)
      positive but small   it ranks weakly, and the band says how weakly

THE PROBE, AND WHY IT IS A TRIPEPTIDE
    ALA90-HIS91-VAL92 lifted intact from 1TNF chain A and rigid-body placed. Three residues, not
    one, because a SINGLE-residue chain is both N- and C-terminal at once and ff14SB treats a
    terminal histidine as a different residue carrying OXT -- that exact fault already fired
    dddg_elec's net-charge guard once ("REFUSING: H:HIS219 net charge protonated -0.0000"). With
    flanking residues the histidine is internal and parametrises as standard HIP, net +1.

    The probe's own charged termini do NOT contaminate the result, and this is structural rather
    than lucky: dddG_elec is the difference between two charge vectors that differ ONLY on
    histidine atoms. Every non-histidine charge, termini included, cancels EXACTLY -- which is
    what dddg_elec's ZERO_TOL guard asserts to machine precision.

    No bond length, angle or charge in this file is invented. The probe is real coordinates, the
    placement is rigid, and the scoring is bin/dddg_elec.py unmodified, called as a subprocess so
    there is no second implementation to drift.

THE MATCHED FAR CONTROL, WHICH IS THE POINT
    A positive number at the anchor means nothing on its own -- "any histidine near a protein
    scores positive" would look identical. So the same probe is also placed at an interface patch
    with NO cationic nitrogen within FAR_MIN of its ND1. Same probe, same chemistry, same
    scorer, same band; only the partner changes. If the anchor and the far control score alike,
    the filter is reading background and not mechanism.

CONSTANTS
    D_PLACE 3.0 A    ND1 to the engaged donor nitrogen
    CLASH   3.0 A    heavy-atom floor against TNF, the engaged donor exempt
    FAR_MIN 10.0 A   how far the far control's ND1 must be from every cationic nitrogen.
                     NOT the census's 12 A: measured on this structure, NO TNF heavy atom is
                     more than 14.8 A from a cationic nitrogen and only 48 of 3,552 clear 12 A,
                     so a 12 A control is not constructible here. 10 A is the right bar anyway,
                     because dddG_elec truncates every pair at d_max = 5.5 A -- a partner at
                     10 A contributes EXACTLY ZERO to the score, by the scorer's own definition.
                     That TNF's surface is this cation-dense is itself worth knowing: "put the
                     histidine away from cations" is nearly impossible on this target, which is
                     also why s6b's scramble negative was never a clean negative.
    CONTACT 4.5 A    the far control must still TOUCH TNF, or there is no interface to score
"""

import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import fetch                                                      # noqa: E402
import his_cation_gate as gate                                     # noqa: E402
from anchor_reach import (fibonacci, rot_axis, rot_between,        # noqa: E402
                          N_DIR, N_SPIN, LOCAL_R)
from third_site_census import tnf_chain_offset                     # noqa: E402

D_PLACE = 3.0
CLASH = 3.0
FAR_MIN = 10.0
CONTACT = 4.5
PROBE_RESNUMS = (14, 15, 16)          # ALA90 HIS91 VAL92 in 1TNF's own numbering
N_POSE = 5                            # playbook s18: n >= 5 poses, never one
CATION_N = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",), "HIS": ("ND1", "NE2")}
VENV = os.path.join(ROOT, ".venv", "bin", "python")
PY = VENV if os.path.exists(VENV) else sys.executable


def load():
    import gemmi
    st = gemmi.read_structure(fetch.pdb("1TNF"))
    st.setup_entities()
    st.remove_ligands_and_waters()
    while len(st) > 1:
        del st[1]
    model = st[0]
    probe, coords, cations, byres = [], [], [], {}
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
                if ch.name == "A" and r.seqid.num in PROBE_RESNUMS:
                    probe.append((r.seqid.num, r.name, at.name, p))
            byres[(ch.name, u)] = (r.name, d)
            for an in CATION_N.get(r.name, ()):
                if an in d:
                    cations.append(d[an])
    return st, probe, np.asarray(coords), np.asarray(cations), byres


def placed(probe, lp_ref, aim, u, spin):
    """Rigid-body the whole tripeptide so its ND1 sits at aim + u*D_PLACE with the ND1 lone
    pair pointing back along -u. Returns the transformed atom list."""
    nd1 = next(p for n, rn, an, p in probe if rn == "HIS" and an == "ND1")
    r = rot_axis(-u, spin) @ rot_between(lp_ref, -u)
    shift = aim + u * D_PLACE
    return [(n, rn, an, (p - nd1) @ r.T + shift) for n, rn, an, p in probe]


def best_placement(probe, lp_ref, coords, aim, donor_pos, require_far=None, cations=None,
                   n_want=1, min_sep_deg=30.0):
    """Up to n_want clear placements, highest clearance first and mutually at least
    min_sep_deg apart in approach direction.

    n_want > 1 exists because playbook s18 says n >= 5 poses, and s6 already established on
    this target that a single pose cannot rank anything. The angular separation is what stops
    the five "poses" being five views of the same one.

    Deterministic by construction, so the reported number is not whichever direction the loop
    happened to reach first."""
    near = coords[np.linalg.norm(coords - aim, axis=1) <= LOCAL_R + 8.0]
    if donor_pos is not None:
        near = near[np.linalg.norm(near - donor_pos, axis=1) > 1e-6]
    found = []
    for u in fibonacci(N_DIR):
        for s in np.linspace(0, 2 * math.pi, N_SPIN, endpoint=False):
            pl = placed(probe, lp_ref, aim, u, s)
            pos = np.array([p for _, _, _, p in pl])
            dmin = float(np.min(np.linalg.norm(
                pos[:, None, :] - near[None, :, :], axis=2)))
            if dmin < CLASH:
                continue
            nd1 = next(p for _, rn, an, p in pl if rn == "HIS" and an == "ND1")
            if require_far is not None:
                if float(np.min(np.linalg.norm(cations - nd1, axis=1))) < require_far:
                    continue
                if float(np.min(np.linalg.norm(coords - nd1, axis=1))) > CONTACT:
                    continue
            found.append((dmin, u, pl))
    found.sort(key=lambda x: -x[0])
    keep, dirs = [], []
    cos_lim = math.cos(math.radians(min_sep_deg))
    for dmin, u, pl in found:
        if all(float(np.dot(u, v)) < cos_lim for v in dirs):
            keep.append((pl, dmin))
            dirs.append(u)
            if len(keep) >= n_want:
                break
    if n_want == 1:
        return keep[0] if keep else (None, -1.0)
    return keep


def write_pose(st, pl, path):
    import gemmi
    out = gemmi.Structure()
    out.spacegroup_hm = "P 1"
    m = gemmi.Model("1")
    for ch in st[0]:
        m.add_chain(ch.clone())
    z = gemmi.Chain("Z")
    # gemmi's Chain.add_residue() appends a COPY. The first version of this function created
    # each Residue, added it to the chain, then added atoms to the ORIGINAL object -- so chain Z
    # was written with zero atoms, an empty chain was dropped, and dddg_elec refused with
    # "chains ['Z'] not in structure". Build each residue complete, THEN add it.
    grouped = {}
    for n, rn, an, pp in pl:
        grouped.setdefault(n, (rn, []))[1].append((an, pp))
    for n in PROBE_RESNUMS:
        if n not in grouped:
            continue
        rn, ats = grouped[n]
        r = gemmi.Residue()
        r.name = rn
        r.seqid = gemmi.SeqId(str(PROBE_RESNUMS.index(n) + 1))
        for an, pp in ats:
            a = gemmi.Atom()
            a.name = an
            a.element = gemmi.Element("N" if an.startswith("N") else
                                      "O" if an.startswith("O") else "C")
            a.pos = gemmi.Position(*pp)
            a.occ = 1.0
            a.b_iso = 20.0
            r.add_atom(a)
        z.add_residue(r)
    m.add_chain(z)
    out.add_model(m)
    out.setup_entities()
    out.write_pdb(path)
    # s24: assert the file on disk actually carries the probe. The silent-empty-chain bug above
    # produced a perfectly valid PDB that simply did not contain the thing being tested.
    back = gemmi.read_structure(path)
    nz = sum(len(r) for ch in back[0] if ch.name == "Z" for r in ch)
    if nz != len(pl):
        raise SystemExit(f"REFUSING: wrote {path} but chain Z holds {nz} atoms, expected "
                         f"{len(pl)}. The probe is not in the file.")
    return path


def score(path, legs):
    rows = []
    for eps, taut in legs:
        r = subprocess.run([PY, os.path.join(ROOT, "bin", "dddg_elec.py"), path,
                            "--binder", "Z", "--target", "A,B,C",
                            "--eps", eps, "--tautomer", taut],
                           capture_output=True, text=True)
        val = None
        for ln in (r.stdout or "").splitlines():
            m = ln.strip()
            # the line reads "dddG_elec     +0.127   n_his_interface=1   -> PASS" -- capital G,
            # and the number is the SECOND field, not the last
            if m.startswith("dddG_elec"):
                try:
                    val = float(m.split()[1])
                except (ValueError, IndexError):
                    pass
        rows.append((eps, taut, val, r.returncode,
                     (r.stdout or "") + (r.stderr or "")))
    return rows


def main():
    st, probe, coords, cations, byres = load()
    lp_ref = gate.lone_pair(
        {an: p for _, rn, an, p in probe if rn == "HIS"}, "ND1")
    assert lp_ref is not None, "no lone pair on the probe histidine"
    print(f"probe: {'-'.join(rn for n, rn, an, p in probe if an == 'CA')} "
          f"from 1TNF A:{PROBE_RESNUMS}, {len(probe)} heavy atoms; "
          f"{len(coords)} TNF heavy atoms, {len(cations)} cationic nitrogens")

    poses = []
    for label, (ch, uni, an) in {"anchor_R108_NH2": ("B", 108, "NH2"),
                                 "anchor_K166_NZ": ("A", 166, "NZ")}.items():
        rn, d = byres[(ch, uni)]
        got = best_placement(probe, lp_ref, coords, d[an], d[an], n_want=N_POSE)
        if not got:
            print(f"{label}: NO clear placement for the tripeptide")
            continue
        for i, (pl, clr) in enumerate(got):
            poses.append((f"{label}#{i + 1}",
                          write_pose(st, pl, os.path.join(HERE, f"pose_{label}_{i + 1}.pdb")),
                          clr, f"{rn}{uni}:{an} on chain {ch}"))

    # matched far control: same probe, an interface patch with no cation within FAR_MIN
    done = False
    for (ch, uni), (rn, d) in sorted(byres.items()):
        if done or rn in CATION_N:
            continue
        for an, p in d.items():
            if an in ("N", "CA", "C", "O"):
                continue
            if float(np.min(np.linalg.norm(cations - p, axis=1))) < FAR_MIN:
                continue
            got = best_placement(probe, lp_ref, coords, p, None, require_far=FAR_MIN,
                                 cations=cations, n_want=N_POSE)
            if got:
                for i, (pl, clr) in enumerate(got):
                    poses.append((f"far_control#{i + 1}",
                                  write_pose(st, pl,
                                             os.path.join(HERE, f"pose_far_control_{i + 1}.pdb")),
                                  clr, f"aimed at {ch}:{rn}{uni}:{an}, "
                                       f"no cation within {FAR_MIN} A"))
                done = True
                break

    legs = [("sigmoid", "HIE"), ("sigmoid", "HID"), ("10", "HIE"), ("6", "HIE"), ("80", "HIE")]
    print(f"\n{'pose':<20}{'clearance':>10}  partner")
    print("-" * 78)
    for label, path, clr, what in poses:
        print(f"{label:<20}{clr:>9.2f}A  {what}")

    print(f"\n{'pose':<20}" + "".join(f"{e}/{t:<4}"[:12].ljust(12) for e, t in legs))
    print("-" * 78)
    results = {}
    for label, path, clr, what in poses:
        rows = score(path, legs)
        vals = [v for _, _, v, _, _ in rows]
        results[label] = vals
        cells = "".join((f"{v:+.3f}" if v is not None else "ERR").ljust(12) for v in vals)
        print(f"{label:<20}{cells}")
        for eps, taut, v, rc, log in rows:
            if v is None:
                print(f"    {eps}/{taut} returned no number (rc={rc}):")
                print("    " + "\n    ".join(log.strip().splitlines()[-6:]))

    # playbook s18: report the spread across poses, never a single number. s6 established on
    # this target that one pose cannot rank anything.
    groups = {}
    for lab, v in results.items():
        groups.setdefault(lab.split("#")[0], []).append(v)

    print(f"\nMEDIAN [min..max] across poses, per leg")
    print("-" * 78)
    for cond, vv in groups.items():
        print(f"{cond}  (n={len(vv)})")
        for i, (e, t) in enumerate(legs):
            col = [x[i] for x in vv if x[i] is not None]
            if col:
                print(f"    {e+'/'+t:<14}{np.median(col):+.3f}  "
                      f"[{min(col):+.3f} .. {max(col):+.3f}]"
                      + ("   SIGN UNSTABLE" if min(col) < 0 <= max(col) else ""))
            else:
                print(f"    {e+'/'+t:<14}ERR")

    # The filter is `dddG_elec >= 0`, so the only thing it actually consumes is the SIGN.
    print(f"\nWhat the `>= 0` filter would do, over all {len(legs)} legs x {N_POSE} poses:")
    for cond, vv in groups.items():
        allv = [x for row in vv for x in row if x is not None]
        pos = sum(1 for x in allv if x >= 0)
        print(f"  {cond:<20}{pos:>3}/{len(allv)} PASS")
    print("\nIf the far control passes at a rate like the anchors', the filter is reading\n"
          "background and not mechanism -- which is the whole question this file was built to\n"
          "answer, and it is answerable only because the far control exists.")
    return results


if __name__ == "__main__":
    main()
