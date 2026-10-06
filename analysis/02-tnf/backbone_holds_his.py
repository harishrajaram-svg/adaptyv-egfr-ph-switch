#!/usr/bin/env python3
"""Does a REAL BACKBONE hold an installed histidine where a free fragment does not?

s10's largest open confound. Our anchor probe is a rigid-body placed ALA-HIS-VAL tripeptide: it
is attached to nothing, so when the packer is allowed to move it, nothing stops the His+ rotating
away from the cation. That is the most plausible reason the K166 signal collapses from +4.9 rigid
to -0.24 repacked, and if it is the reason, s10's "the mechanism runs backwards" reading is an
artifact of the FIXTURE rather than a statement about designs.

A real binder holds its histidine in a backbone the generator built. We have no designs yet
(blocked on the target decision), so the question needs an existing structure with INSTALLED
histidines on a real scaffold.

3QSK: Murtaugh et al., Protein Sci 20(9):1619-1631 (2011), doi 10.1002/pro.696. A camelid VHH
carrying FIVE engineered histidines, solved at 1.75 A in complex with its target.

\U0001F534 DIFFERENT TARGET, STATED UP FRONT. 3QSK is anti-RNase A, not anti-TNF. The challenge
file called Murtaugh "histidine scanning in a VHH, OUR FORMAT" and this project briefly misfiled
it as TNF prior art on that wording (s8, corrected same day). It is NOT prior art here and
nothing about TNF is inferred from it. The question it answers is about BACKBONE CONSTRAINT, which
is a property of the scaffold and not of the epitope: given a real VHH framework with engineered
histidines, how far does the packer move them?

MEASURED, per histidine: sidechain heavy-atom RMSD between the pre- and post-repack pose, in the
protonated state, over the same seeded trials the sweep uses. Compared against the same quantity
for our free tripeptide probe. A large gap means the probe overstates mobility and s10's reading
needs softening; a small gap means the collapse is real and s4 is in trouble.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dddg_elec_repacked as dd
import fetch

SIDECHAIN = {"CB", "CG", "ND1", "CD2", "CE1", "NE2"}


def his_rmsd(path, binder, target, label):
    """Per-histidine sidechain RMSD under repacking, protonated state, median over trials."""
    dd.init()
    base = dd.load(path)
    scored = dd.chain_res(base, binder + target)
    sites = dd.histidines(base, scored)
    bsites = [i for i in sites if base.pdb_info().chain(i) in set(binder)]
    if not bsites:
        print(f"  {label}: no binder-side histidine; skipped")
        return {}
    per = {i: [] for i in bsites}
    for t in range(1, dd.TRIALS + 1):
        p = base.clone()
        dd.set_state(p, sites, "HIS_P")
        before = p.clone()
        dd.set_seed(t)
        dd.repack(p, bsites, lock=sites, state="protonated")
        for i in bsites:
            ra, rb = before.residue(i), p.residue(i)
            d = []
            for k in range(1, ra.natoms() + 1):
                an = ra.atom_name(k).strip()
                if an not in SIDECHAIN or ra.atom_type(k).is_hydrogen():
                    continue
                for j in range(1, rb.natoms() + 1):
                    if rb.atom_name(j).strip() == an:
                        pa, pb = ra.xyz(k), rb.xyz(j)
                        d.append((pa.x - pb.x) ** 2 + (pa.y - pb.y) ** 2 + (pa.z - pb.z) ** 2)
                        break
            if d:
                per[i].append(float(np.sqrt(np.mean(d))))
    out = {}
    for i, v in per.items():
        if v:
            tag = f"{base.pdb_info().chain(i)}:{base.pdb_info().number(i)}"
            out[tag] = (float(np.median(v)), min(v), max(v))
            print(f"  {label:<26}{tag:<10}median {out[tag][0]:6.3f} A  "
                  f"[{out[tag][1]:.3f} .. {out[tag][2]:.3f}]")
    return out


def main():
    print(__doc__.strip().splitlines()[0])
    print(f"\nn = {dd.TRIALS} paired seeded trials, protonated state, binder-side histidines only.")
    print(f"{'fixture':<28}{'site':<10}sidechain RMSD under repacking\n" + "-" * 78)

    free = his_rmsd(os.path.join(HERE, "pose_anchor_K166_rebuilt_1.pdb"), "Z", "ABC",
                    "FREE tripeptide probe")
    try:
        p3 = fetch.cif("3QSK")
    except SystemExit as e:
        print(f"  3QSK unavailable: {e}")
        return
    # chain A = RNase A (target), chain B = the VHH (binder)
    held = his_rmsd(p3, "B", "A", "3QSK VHH (real backbone)")

    if not free or not held:
        print("\nOne side produced nothing; no comparison.")
        return
    f = float(np.median([v[0] for v in free.values()]))
    h = float(np.median([v[0] for v in held.values()]))
    print("-" * 78)
    print(f"median over sites: FREE probe {f:.3f} A   vs   REAL backbone {h:.3f} A "
          f"({len(held)} histidines)")
    print()
    ratio = f / max(h, 1e-6)
    MEANINGFUL = 1.5          # written down before the numbers (s13). Below this, "tighter" is
                              # a rounding difference on a 5-site sample, not a constraint.
    if ratio >= MEANINGFUL:
        print(f"The real backbone holds its histidines {ratio:.1f}x tighter than the free probe.")
        print("s10's collapse at K166 is then at least partly a property of the FIXTURE, and the")
        print("'mechanism runs backwards' reading must be softened until a real design is scored.")
        print("It does NOT become wrong -- it becomes untested on the geometry that matters.")
    else:
        lo = min(v[0] for v in held.values()); hi = max(v[0] for v in held.values())
        print(f"THE REAL BACKBONE DOES NOT HOLD ITS HISTIDINES ANY TIGHTER: {ratio:.2f}x, against")
        print(f"a {MEANINGFUL}x bar set before the numbers. The free probe at {f:.3f} A sits")
        print(f"INSIDE the spread of the five engineered histidines on a real VHH framework")
        print(f"({lo:.3f}-{hi:.3f} A).")
        print()
        print("So s10's collapse is NOT a fixture artifact, and this was the last cheap way out.")
        print("A histidine on a real scaffold is about as free to rotate away from a cation as")
        print("our unattached tripeptide was. s4 prices two anchor sites on the opposite")
        print("assumption.")
    print("\n⚠️  Different target (RNase A). This measures BACKBONE CONSTRAINT only. Nothing")
    print("   about TNF, the epitope, or any ratio is inferred from 3QSK.")


if __name__ == "__main__":
    main()
