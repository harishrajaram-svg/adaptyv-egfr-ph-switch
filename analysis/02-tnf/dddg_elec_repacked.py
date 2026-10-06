#!/usr/bin/env python3
"""
Ahn's ACTUAL protocol: dddG_elec with the pose REPACKED at each protonation state.

WHY THIS EXISTS, AND WHY IT COULD NOT EXIST BEFORE TODAY
    bin/dddg_elec.py's docstring says Rosetta/PyRosetta are "out of scope on licence grounds by
    the organisers' own rule". That reading is explicitly wrong for non-commercial use -- Amir
    Shanehsazzadeh, #design-methods, 2026-10-06 08:56 EDT. So the reimplementation was not
    wasted (it produced s6b, s6e and s9) but it was not necessary either.

    Ahn et al. REPACK the pose at both protonation states. Ours is RIGID, and dddg_elec.py's own
    docstring calls the relaxed leg "the honest one": repacking lets a His+ ROTATE AWAY from a
    cation instead of being frozen against it, which is the single most plausible way the rigid
    number is wrong.

THIS IS A NEW LEG, NOT AN EDIT. bin/dddg_elec.py is untouched -- 7 self-tests, 2 mutation tests,
and every published figure in s6b/s6e/s9 came out of it. Keeping them separate means the two
legs can be COMPARED, which is the point, and a bug here cannot retroactively move an old number.

HOW THE STATES ARE SET
    Rosetta's doubly-protonated histidine is a BASE residue type, HIS_P, not a variant -- so
    add_variant_type_to_pose_residue(..., PROTONATED) fails with "Unable to find desired residue
    'HIS' with variant 'PROTONATED'", which cost twenty minutes. The residue is REPLACED with the
    HIS_P type instead, coordinates copied. HIS_P is only in the residue type set under
    -pH_mode, which is why Ahn's Methods mention pH mode at all.

WHICH HISTIDINES SWITCH: ALL OF THEM IN THE SCORED CHAINS, matching the rigid tool exactly.
    bin/dddg_elec.py's prepare() substitutes neutral charges for EVERY histidine in the kept
    chains (its by_res loop), not only the interface ones -- interface_his() is reporting and the
    G1 condition, not the partition's domain. On a probe pose that is the probe His PLUS TNF's
    own H91/H149/H154 in all three protomers: 10 histidines, not 1. Switching a subset here
    would make the two legs incomparable, which is s16 (the reference state is a result).

WHAT IS SCORED
    The cross-interface fa_elec sum, read off Rosetta's own energy graph per residue pair. That
    is the same quantity the rigid tool computes: for a rigid separation the intramolecular terms
    cancel exactly, so ddG_elec reduces to the cross-interface pair sum. No pose surgery, so no
    chance of a split that silently drops a chain.

    REPACKING USES ref2015, SCORING USES fa_elec ALONE. Repacking on electrostatics only would
    place sidechains by a criterion no one uses; Ahn repack with Rosetta's standard function.
    Reported separately so the choice is visible.

CONSTANTS
    SHELL 8.0 A   repack sidechains within this of a shell CENTRE, everything else fixed.
                  LABELLED OURS -- Ahn publish no shell.

                  TWO SHELLS ARE REPORTED, because the first one measured was not attributable.
                  Centring on all ten switched histidines made 264 of 459 residues repackable
                  and moved sidechains 1.195 A RMSD -- the trimer carries H91/H149/H154 in each
                  protomer, so "local" around all of them is most of the protein, and any
                  dddG change would be global rearrangement rather than the probe's response.
                    repacked_local  shell around the BINDER histidine only. The designable one,
                                    and the one P2 is about: can the probe His+ rotate away.
                    repacked_all    shell around every switched histidine. The faithful reading
                                    of "repack the pose", carried as a sensitivity leg (s16).
    Backbone is NEVER moved. A backbone change is a different pose, not a response.

    dddg_elec_repacked.py --selftest
    dddg_elec_repacked.py <pose.pdb> --binder Z --target A,B,C [--rigid-only]
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

SHELL = 8.0
HIS_TYPES = ("HIS", "HIS_D", "HIS_P")
INIT = "-mute all -ex1 -ex2aro -pH_mode true -value_pH 0 -ignore_unrecognized_res false"
_READY = {"done": False}


def init():
    if _READY["done"]:
        return
    import pyrosetta
    pyrosetta.init(INIT, silent=True)
    _READY["done"] = True


def _ts():
    from pyrosetta.rosetta.core.chemical import ChemicalManager
    return ChemicalManager.get_instance().residue_type_set("fa_standard")


def load(path):
    init()
    import pyrosetta
    if not os.path.exists(path):
        raise SystemExit(f"REFUSING: {path} does not exist")
    return pyrosetta.pose_from_pdb(path)


def chain_res(pose, chains):
    pi = pose.pdb_info()
    return [i for i in range(1, pose.total_residue() + 1) if pi.chain(i) in chains]


def histidines(pose, scored):
    """Every histidine in the scored chains, in any protonation state."""
    return [i for i in scored if pose.residue(i).name3() == "HIS"
            or pose.residue(i).name().split(":")[0] in HIS_TYPES]


def set_state(pose, sites, state):
    """Replace each site's residue type with HIS (neutral, NE2-H) or HIS_P (doubly protonated)."""
    from pyrosetta.rosetta.core.pose import \
        replace_pose_residue_copying_existing_coordinates as repl
    ts = _ts()
    rt = ts.name_map(state)
    for i in sites:
        repl(pose, i, rt)


def net_charge(pose, i):
    r = pose.residue(i)
    return sum(r.atomic_charge(k) for k in range(1, r.natoms() + 1))


def cross_fa_elec(pose, A, B):
    """Cross-interface fa_elec, summed off Rosetta's energy graph. Returns (energy, n_pairs)."""
    from pyrosetta.rosetta.core.scoring import ScoreType, ScoreFunction
    sf = ScoreFunction()
    sf.set_weight(ScoreType.fa_elec, 1.0)
    sf(pose)
    eg = pose.energies().energy_graph()
    tot, n = 0.0, 0
    Bs = set(B)
    for i in A:
        for j in Bs:
            e = eg.find_energy_edge(i, j)
            if e is not None:
                tot += e.fill_energy_map()[ScoreType.fa_elec]
                n += 1
    return tot, n


def repack(pose, sites, shell=SHELL):
    """Repack sidechains within `shell` of any switched histidine. Backbone fixed.

    Returns the number of residues allowed to move, so a run that repacked nothing is visible
    rather than silently labelled 'repacked'."""
    from pyrosetta import get_fa_scorefxn
    from pyrosetta.rosetta.core.pack.task import TaskFactory, operation
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    import numpy as np

    centres = []
    for i in sites:
        r = pose.residue(i)
        for k in range(1, r.natoms() + 1):
            if not r.atom_type(k).is_hydrogen():
                centres.append(np.array([r.xyz(k).x, r.xyz(k).y, r.xyz(k).z]))
    centres = np.asarray(centres)

    movable = []
    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        pts = np.array([[r.xyz(k).x, r.xyz(k).y, r.xyz(k).z]
                        for k in range(1, r.natoms() + 1)
                        if not r.atom_type(k).is_hydrogen()])
        if pts.size and float(np.min(np.linalg.norm(
                pts[:, None, :] - centres[None, :, :], axis=2))) <= shell:
            movable.append(i)

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())
    task = tf.create_task_and_apply_taskoperations(pose)
    for i in range(1, pose.total_residue() + 1):
        if i not in movable:
            task.nonconst_residue_task(i).prevent_repacking()
    PackRotamersMover(get_fa_scorefxn(), task).apply(pose)
    return len(movable)


def sidechain_rmsd(a, b, residues):
    import numpy as np
    d = []
    for i in residues:
        ra, rb = a.residue(i), b.residue(i)
        if ra.name() != rb.name():
            continue
        for k in range(1, min(ra.natoms(), rb.natoms()) + 1):
            if ra.atom_type(k).is_hydrogen() or ra.atom_name(k).strip() in ("N", "CA", "C", "O"):
                continue
            pa, pb = ra.xyz(k), rb.xyz(k)
            d.append((pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2)
    return float(np.sqrt(np.mean(d))) if d else 0.0


def score(path, binder, target, rigid_only=False, shell=SHELL):
    """Both legs on one pose. Returns a dict; raises rather than returning a silent zero."""
    scored = sorted(set(binder) | set(target))
    out = {"path": os.path.basename(path)}
    base = load(path)
    A = chain_res(base, binder)
    B = chain_res(base, target)
    if not A or not B:
        raise SystemExit(f"REFUSING {path}: binder side {len(A)} residues, target side {len(B)}")
    sites = histidines(base, chain_res(base, scored))
    if not sites:
        raise SystemExit(f"REFUSING {path}: no histidine in chains {scored}. dddG_elec is "
                         f"undefined without one -- it is a histidine filter.")
    out["n_his"] = len(sites)
    out["his"] = [f"{base.pdb_info().chain(i)}:{base.pdb_info().number(i)}" for i in sites]

    binder_sites = [i for i in sites if base.pdb_info().chain(i) in set(binder)]
    out["n_his_binder"] = len(binder_sites)
    centres = {"repacked_local": binder_sites or sites, "repacked_all": sites}

    legs = {}
    leg_names = ("rigid",) if rigid_only else ("rigid", "repacked_local", "repacked_all")
    for state in ("HIS", "HIS_P"):
        for leg in leg_names:
            p = base.clone()
            set_state(p, sites, state)
            q = [net_charge(p, i) for i in sites]
            want = 1.0 if state == "HIS_P" else 0.0
            bad = [(s, v) for s, v in zip(out["his"], q) if abs(v - want) > 1e-3]
            if bad:
                raise SystemExit(f"REFUSING {path}: net charge wrong in state {state}: {bad[:3]}")
            moved = None
            if leg != "rigid":
                before = p.clone()
                n_mov = repack(p, centres[leg], shell)
                moved = (n_mov, sidechain_rmsd(before, p, range(1, p.total_residue() + 1)))
            e, n = cross_fa_elec(p, A, B)
            legs[(leg, state)] = (e, n, moved)

    for leg in leg_names:
        en, nn, _ = legs[(leg, "HIS")]
        ep, np_, mv = legs[(leg, "HIS_P")]
        out[leg] = {"neutral": round(en, 4), "protonated": round(ep, 4),
                    "dddG": round(ep - en, 4), "n_pairs": nn,
                    "repacked_residues": mv[0] if mv else None,
                    "sidechain_rmsd": round(mv[1], 3) if mv else None}
    return out


def selftest():
    ok = []
    k166 = os.path.join(HERE, "pose_anchor_K166_rebuilt_1.pdb")
    k166q = os.path.join(HERE, "pose_anchor_K166Q_1.pdb")
    for f in (k166, k166q):
        if not os.path.exists(f):
            raise SystemExit(f"REFUSING: {f} missing. Run anchor_dddg_k166q.py first.")

    a = score(k166, "Z", "A,B,C")
    b = score(k166q, "Z", "A,B,C")

    # T1 -- ALL histidines switch, not a subset, and the count is what the structure implies:
    # the probe His plus TNF's H91/H149/H154 in three protomers.
    assert a["n_his"] == 10, f"T1 FAIL {a['n_his']} histidines switched, expected 10: {a['his']}"
    ok.append(f"T1 all {a['n_his']} histidines in the scored chains switch, matching the rigid "
              f"tool's domain, not just the interface one: {a['his']}")

    # T2 -- net charge per state, the same guard the rigid tool carries. score() raises on
    # failure, so reaching here IS the assertion; stated so it is visible in the output.
    ok.append("T2 every switched histidine is net +1 protonated and net 0 neutral, asserted "
              "per site per state (score() refuses otherwise)")

    # T3 -- KNOWN ANSWER, cross-implementation. s9 established with bin/dddg_elec.py that K166
    # scores ABOVE K166Q on all five band legs. Real Rosetta must agree on the ORDERING. It
    # need not agree on magnitude -- different dielectric model, different weights.
    assert a["rigid"]["dddG"] > b["rigid"]["dddG"], \
        f"T3 FAIL rigid: K166 {a['rigid']['dddG']} is not above K166Q {b['rigid']['dddG']}"
    assert a["rigid"]["dddG"] > 0 > b["rigid"]["dddG"], \
        f"T3 FAIL rigid signs: K166 {a['rigid']['dddG']}, K166Q {b['rigid']['dddG']}"
    ok.append(f"T3 known answer, CROSS-IMPLEMENTATION: real Rosetta's rigid leg reproduces s9's "
              f"ordering and both signs — K166 {a['rigid']['dddG']:+.3f} vs K166Q "
              f"{b['rigid']['dddG']:+.3f}")

    # T4 -- the repack actually moved sidechains. Without this, "repacked" is a label.
    r = a["repacked_local"]
    assert r["repacked_residues"] and r["repacked_residues"] > 1, \
        f"T4 FAIL only {r['repacked_residues']} residues were repackable"
    assert r["sidechain_rmsd"] > 0.01, \
        f"T4 FAIL sidechain RMSD {r['sidechain_rmsd']} — nothing moved, so this is not a leg"
    ok.append(f"T4 the PROBE-LOCAL repack moved something: {r['repacked_residues']} residues "
              f"in the {SHELL} A shell around the binder histidine, sidechain RMSD "
              f"{r['sidechain_rmsd']} A. The all-histidine shell reaches "
              f"{a['repacked_all']['repacked_residues']} residues at "
              f"{a['repacked_all']['sidechain_rmsd']} A, which is most of the trimer and is "
              f"why the two legs are reported separately rather than averaged")

    # T5 -- MUTATION TEST. A pose with no histidine must REFUSE, not return 0.0. s10: a zero is
    # not an absence, and this exact failure mode (fail-open to a clean zero) is s19.
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        stripped = os.path.join(td, "nohis.pdb")
        with open(stripped, "w") as fh:
            for ln in open(k166):
                if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() in ("HIS", "HIP"):
                    continue
                fh.write(ln)
        fired = False
        try:
            score(stripped, "Z", "A,B,C")
        except SystemExit as e:
            fired = "no histidine" in str(e) or "net charge" in str(e) or "REFUSING" in str(e)
        assert fired, "T5 FAIL a histidine-free pose did not refuse"
    ok.append("T5 mutation test: a histidine-stripped pose REFUSES rather than returning 0.0")

    for line in ok:
        print("  ok  " + line)
    print(f"\nself-tests passed: {len(ok)}")
    return a, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pose", nargs="?")
    ap.add_argument("--binder", default="Z")
    ap.add_argument("--target", default="A,B,C")
    ap.add_argument("--rigid-only", action="store_true")
    ap.add_argument("--shell", type=float, default=SHELL)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest or not a.pose:
        selftest()
        return
    r = score(a.pose, a.binder.replace(",", ""), a.target.replace(",", ""),
              a.rigid_only, a.shell)
    print(f"{r['path']}  {r['n_his']} histidines switched")
    for leg in ("rigid", "repacked_local", "repacked_all"):
        if leg in r:
            d = r[leg]
            extra = (f"  [{d['repacked_residues']} res repacked, "
                     f"sc-RMSD {d['sidechain_rmsd']} A]" if d["repacked_residues"] else "")
            print(f"  {leg:<15} neutral {d['neutral']:+9.4f}  protonated {d['protonated']:+9.4f}"
                  f"  dddG {d['dddG']:+9.4f}  ({d['n_pairs']} pairs){extra}")


if __name__ == "__main__":
    main()
