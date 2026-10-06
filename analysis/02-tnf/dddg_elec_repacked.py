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
# The two CHARGE states, and the residue types that realise each. HIS (NE2-H, "HIE") and HIS_D
# (ND1-H, "HID") are both net-neutral tautomers -- so locking the neutral leg to the literal
# name "HIS" is wrong, and the guard caught it: the packer returned HIS_D and the run REFUSED.
# The invariant is the NET CHARGE, not the name. Letting Rosetta choose between the two neutral
# tautomers is a feature, not slack: the rigid tool has to GUESS the tautomer and carries
# HIE/HID as two separate band legs (s16). Here the packer picks, and which one it picks is
# reported.
STATE_TYPES = {"neutral": ("HIS", "HIS_D"), "protonated": ("HIS_P",)}
STATE_CHARGE = {"neutral": 0.0, "protonated": 1.0}
SET_TYPE = {"neutral": "HIS", "protonated": "HIS_P"}
INIT = ("-mute all -ex1 -ex2aro -pH_mode true -value_pH 0 -ignore_unrecognized_res false "
        "-constant_seed")
TRIALS = 5          # playbook s18: n >= 5, because the packer is STOCHASTIC -- see set_seed()
BASE_SEED = 20261006
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


def set_seed(k):
    """Pin Rosetta's RNG for one packer trajectory.

    🔴 THE PACKER IS STOCHASTIC AND THIS WAS MEASURED, NOT ASSUMED. Three repacks of the
    SAME pose in the SAME protonation state gave cross-interface energies -4.6869, -4.3795 and
    -4.4282 -- a spread of 0.31 on a quantity whose repacked dddG values sit at -0.24. So a
    single repack pair carries roughly 0.3-0.4 of packer noise, which is LARGER than most of the
    differences the sweep was reporting. Every repacked number from a single trajectory is n=1 on
    a random process, which is s18 and the exact trap s6e fell into at n=1.

    Trials are paired: trial k repacks the neutral and protonated states under the SAME seed, so
    the difference is not also differencing two independent draws."""
    from pyrosetta.rosetta.numeric.random import rg
    rg().set_seed("mt19937", BASE_SEED + k)


def repack(pose, sites, shell=SHELL, lock=None, state=None):
    """Repack sidechains within `shell` of any switched histidine. Backbone fixed.

    Returns the number of residues allowed to move, so a run that repacked nothing is visible
    rather than silently labelled 'repacked'.

    🔴 `lock` + `state` ARE NOT OPTIONAL IN PRACTICE, AND THE FIRST VERSION OMITTED THEM.
    Without them the packer QUIETLY UNDOES THE PROTONATION ASSIGNMENT. `RestrictToRepacking`
    forbids changing the amino acid, but HIS / HIS_D / HIS_P are three residue types of the same
    amino acid, so under -pH_mode the packer is free to pick among them -- that freedom is what
    pH mode is for. Measured on pose_anchor_K166_rebuilt_1: starting from all-neutral, the packer
    re-protonated residue 458 (the probe) back to HIS_P, so BOTH legs converged on the same pose
    and dddG came out EXACTLY +0.0000 on four of five poses.

    That zero read as "the anchor signal vanishes under repacking", which is a dramatic
    conclusion and was false. It is s17 (know the floor your method returns when nothing
    happens), s10 (a zero is not an absence) and s19 (fail closed) arriving together. The tell
    was that the neutral and protonated CROSS-INTERFACE ENERGIES were identical to four
    decimals, not that dddG was zero -- which is why both are now reported, not just the
    difference.

    `restrict_restypes` pins each switched site to the single chosen type, so the packer may
    rotate the sidechain and may not re-title it."""
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
    if lock:
        from pyrosetta.rosetta.utility import vector1_std_string
        if state is None:
            raise SystemExit("REFUSING: repack() got `lock` without `state`; pinning the "
                             "protonation state is the whole point of the argument")
        allowed = vector1_std_string()
        for t in STATE_TYPES[state]:
            allowed.append(t)
        for i in lock:
            task.nonconst_residue_task(i).restrict_restypes(allowed)
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


def score(path, binder, target, rigid_only=False, shell=SHELL, trials=TRIALS,
          legs_wanted=None):
    """Every leg on one pose. Raises rather than returning a silent zero.

    Repacked legs run `trials` PAIRED trajectories: trial k repacks the neutral and the
    protonated state under the SAME seed, so the difference is not also differencing two
    independent draws from a stochastic packer. Reported as median [min..max] over trials,
    never as a single number -- s18."""
    import numpy as np
    scored = sorted(set(binder) | set(target))
    out = {"path": os.path.basename(path), "trials": trials, "shell_A": shell}
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

    leg_names = ("rigid",) if rigid_only else ("rigid", "repacked_local", "repacked_all")
    if legs_wanted:
        leg_names = tuple(l for l in leg_names if l in legs_wanted)

    def one(leg, state, trial):
        p = base.clone()
        set_state(p, sites, SET_TYPE[state])
        want = STATE_CHARGE[state]
        bad = [(s_, v) for s_, v in zip(out["his"], (net_charge(p, i) for i in sites))
               if abs(v - want) > 1e-3]
        if bad:
            raise SystemExit(f"REFUSING {path}: net charge wrong in state {state}: {bad[:3]}")
        moved = None
        if leg != "rigid":
            before = p.clone()
            set_seed(trial)
            n_mov = repack(p, centres[leg], shell, lock=sites, state=state)
            bad_q = [(f"{p.pdb_info().chain(i)}:{p.pdb_info().number(i)}",
                      p.residue(i).name(), round(net_charge(p, i), 3))
                     for i in sites if abs(net_charge(p, i) - want) > 1e-3]
            if bad_q:
                raise SystemExit(f"REFUSING {path}: after repacking, {bad_q[:4]} are not in the "
                                 f"{state} charge state. The two legs would be the same pose and "
                                 f"dddG would be a spurious zero.")
            out.setdefault("tautomers", {}).setdefault(f"{leg}/{state}", set()).update(
                p.residue(i).name().split(":")[0] for i in sites)
            moved = (n_mov, sidechain_rmsd(before, p, range(1, p.total_residue() + 1)))
        e, n = cross_fa_elec(p, A, B)
        return e, n, moved

    for leg in leg_names:
        n_t = 1 if leg == "rigid" else trials
        ds, neus, pros, rms, npair, nmov = [], [], [], [], None, None
        for t in range(1, n_t + 1):
            en, nn, _ = one(leg, "neutral", t)
            ep, np_, mv = one(leg, "protonated", t)
            ds.append(ep - en); neus.append(en); pros.append(ep)
            npair = nn
            if mv:
                nmov, r = mv
                rms.append(r)
        out[leg] = {
            "dddG": round(float(np.median(ds)), 4),
            "dddG_min": round(min(ds), 4), "dddG_max": round(max(ds), 4),
            "dddG_trials": [round(x, 4) for x in ds],
            "sign_unstable": bool(min(ds) < 0 <= max(ds)),
            "neutral": round(float(np.median(neus)), 4),
            "protonated": round(float(np.median(pros)), 4),
            "n_pairs": npair, "n_trials": n_t,
            "repacked_residues": nmov,
            "sidechain_rmsd": round(float(np.median(rms)), 3) if rms else None,
        }
    if "tautomers" in out:
        out["tautomers"] = {k: sorted(v) for k, v in out["tautomers"].items()}
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

    # T6 -- THE PACKER MUST NOT RE-TITRATE. Mutation test on the guard that was missing: run a
    # repack WITHOUT the lock and confirm the protonation state is broken, then WITH it and
    # confirm it holds. Without this test the exact-zero bug comes straight back.
    init()
    base = load(k166)
    sites = histidines(base, chain_res(base, "ZABC"))
    zsites = [i for i in sites if base.pdb_info().chain(i) == "Z"]
    unlocked = base.clone()
    set_state(unlocked, sites, "HIS")
    repack(unlocked, zsites)
    broke = [i for i in sites if abs(net_charge(unlocked, i)) > 1e-3]
    assert broke, ("T6 FAIL an unlocked repack did NOT change any protonation state, so the lock "
                   "cannot be shown to be doing anything")
    locked = base.clone()
    set_state(locked, sites, "HIS")
    repack(locked, zsites, lock=sites, state="neutral")
    held = [i for i in sites if abs(net_charge(locked, i)) > 1e-3]
    taut = sorted({locked.residue(i).name().split(":")[0] for i in sites})
    assert not held, f"T6 FAIL the lock did not hold at {held[:3]}"
    ok.append(f"T6 mutation test on the protonation lock: WITHOUT it the packer re-titrates "
              f"{len(broke)} site(s) (residue {broke[0]}, the probe) to net +1 and both legs "
              f"collapse to the same pose; WITH it all {len(sites)} sites stay net-neutral, and "
              f"Rosetta picks the tautomer itself: {taut}")

    # T7 -- THE PACKER IS STOCHASTIC, AND THE SEED MAKES IT REPRODUCIBLE. Both halves matter.
    # Without the first, a single trajectory looks like a measurement. Without the second,
    # nothing in this file reproduces in a clean clone (s26).
    base2 = load(k166)
    sites2 = histidines(base2, chain_res(base2, "ZABC"))
    z2 = [i for i in sites2 if base2.pdb_info().chain(i) == "Z"]
    A2, B2 = chain_res(base2, "Z"), chain_res(base2, "ABC")

    def one_repack(seed):
        q = base2.clone()
        set_state(q, sites2, "HIS_P")
        set_seed(seed)
        repack(q, z2, lock=sites2, state="protonated")
        return round(cross_fa_elec(q, A2, B2)[0], 4)

    same = [one_repack(1) for _ in range(3)]
    many = [one_repack(k) for k in range(1, 9)]
    assert len(set(same)) == 1, f"T7 FAIL the same seed gave {same} \u2014 not reproducible"
    spread = max(many) - min(many)
    assert len(set(many)) > 1, (f"T7 FAIL eight seeds all gave {many[0]} \u2014 either the seed is "
                                f"not reaching the packer or the {TRIALS}-trial spread is theatre")
    ok.append(f"T7 packer noise is REAL and SEEDED: one seed reproduces exactly ({same[0]}), "
              f"eight seeds span {spread:.4f} on the cross-interface energy "
              f"({len(set(many))} distinct values). Seeds 1\u20133 alone all landed in the same "
              f"basin, so three trials would have read as deterministic \u2014 which is why the "
              f"spread is measured over {TRIALS} and not 3")

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
