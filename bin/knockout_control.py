#!/usr/bin/env python3
"""Causal knockout control for a pH switch: delete the carboxylate, change nothing else.

A positive pH-gate result on a generated complex is not evidence on its own -- a pKa
shift can come from the target's own repacking (see bin/ph_gate_target_his.py). This
isolates the claim "the pinned carboxylate causes the shift" by truncating that one
side chain to Ala in place and re-running PROPKA on otherwise identical coordinates.

Measured 2026-10-03 on g532mimic_short_14 (canonical H433 = mature 409):
    free (binder absent)        6.37
    WT complex (Asp21 present)  7.27   dPKa +0.90   ratio 2.48x
    D21A complex                5.18   dPKa -1.19   ratio 0.65x
    => causal effect of one carboxylate: +2.09 pKa units.
The D21A row matters as much as the WT row: it reproduces the week's negative result
(burial alone suppresses the switch and acid weakens binding), which is why the two
findings are the same system and not a contradiction.

Usage:
    python3 bin/knockout_control.py <complex.cif> <binder_resnum> <target_his_resnum>
e.g.  python3 bin/knockout_control.py runs/.../rank07_..._14.cif 21 409
"""
import os, subprocess, sys, tempfile
import gemmi

PH_LO, PH_HI = 6.5, 7.4
BINDER_CHAIN, TARGET_CHAIN = "A", "B"
KEEP = ("N", "CA", "C", "O", "CB")   # Asp/Glu -> Ala: drop the carboxylate, keep the stub


def ratio(free, bound):
    K = lambda ph: (1 + 10 ** (bound - ph)) / (1 + 10 ** (free - ph))
    return K(PH_LO) / K(PH_HI)


def propka_his(st, wd, tag, chain=TARGET_CHAIN):
    p = os.path.join(wd, tag + ".pdb")
    st.write_pdb(p)
    subprocess.run([sys.executable, "-m", "propka", tag + ".pdb"], capture_output=True, cwd=wd)
    f = os.path.join(wd, tag + ".pka")
    out = {}
    if os.path.exists(f):
        for ln in open(f):
            q = ln.split()
            if len(q) > 3 and q[0] == "HIS" and q[2] == chain and int(q[1]) not in out:
                try: out[int(q[1])] = float(q[3].rstrip("*"))
                except ValueError: pass
    return out


def build(cif, mutate_at, keep_binder):
    st = gemmi.read_structure(cif); st.setup_entities(); st.remove_ligands_and_waters()
    m = st[0]
    if not keep_binder:
        for i in range(len(m) - 1, -1, -1):
            if m[i].name != TARGET_CHAIN: del m[i]
    elif mutate_at is not None:
        ch = [c for c in m if c.name == BINDER_CHAIN][0]
        for r in ch:
            if r.seqid.num == mutate_at:
                # Delete in place, in REVERSE index order. Collecting the wanted
                # Atom objects first and re-adding them holds references into the
                # residue's atom vector while that vector is emptied -- the
                # references dangle and the re-added coordinates are not
                # guaranteed to be the originals. Never hold a reference across
                # the deletion.
                for i in range(len(r) - 1, -1, -1):
                    if r[i].name not in KEEP:
                        del r[i]
                r.name = "ALA"
    st.setup_entities()
    return st


def main():
    if len(sys.argv) != 4:
        print(__doc__); sys.exit(2)
    cif, binder_res, his = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    with tempfile.TemporaryDirectory() as wd:
        free = propka_his(build(cif, None, False), wd, "apo").get(his)
        wt   = propka_his(build(cif, None, True),  wd, "wt").get(his)
        mut  = propka_his(build(cif, binder_res, True), wd, "ko").get(his)
    if None in (free, wt, mut):
        print(f"PROPKA did not report HIS {his} in every leg (free={free} wt={wt} ko={mut})")
        sys.exit(1)
    print(f"HIS {his} pKa_free (binder absent)      : {free:.2f}")
    print(f"HIS {his} pKa_bound, WT                 : {wt:.2f}   dPKa {wt-free:+.2f}   ratio {ratio(free,wt):.2f}x")
    print(f"HIS {his} pKa_bound, residue {binder_res}->Ala      : {mut:.2f}   dPKa {mut-free:+.2f}   ratio {ratio(free,mut):.2f}x")
    print(f"\nCAUSAL EFFECT of residue {binder_res}: {wt-mut:+.2f} pKa units")
    if wt - mut < 0.3:
        print("VERDICT: the carboxylate is NOT the cause -- look for target repacking.")
    else:
        print("VERDICT: the carboxylate is causal.")


if __name__ == "__main__":
    main()
