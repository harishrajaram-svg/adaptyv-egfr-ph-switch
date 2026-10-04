#!/usr/bin/env python3
"""Paired mechanism-B gate: pKa_free from EACH complex's own target coordinates.

USE THIS, not an apo-reference version, for any TARGET-histidine (mechanism B) gate.
bin/ph_gate.py takes a single --his and a CONSTANT --free-pka; that is fine for a fixed
experimental structure but wrong for generated complexes. bin/ph_gate_mechA.py is for
BINDER histidines. This one is for target histidines on generated complexes.

Edit FREE/CANON and the glob at the bottom for a new arm.

The first pass took pKa_free from a single apo structure, reasoning that the target is
rigid. It is not: BoltzGen repacks the target per design, and H370's own suppressor
ASP344 sits 2.4-3.8 A from its ring. A D344 rotamer shift then shows up as a "switch"
in designs whose binder is 27 A away with zero atoms within 6 A of the histidine.
Deleting the binder from the SAME file holds every other coordinate fixed, so the
difference is the binder and nothing else.
"""
import glob, os, sys, tempfile, subprocess, csv
from concurrent.futures import ProcessPoolExecutor
import gemmi

PH_LO, PH_HI = 6.5, 7.4
RATIO_BAR = 1.20
CANON = {346: "H370", 409: "H433"}

def ratio(free, bound):
    K = lambda ph: (1 + 10 ** (bound - ph)) / (1 + 10 ** (free - ph))
    return K(PH_LO) / K(PH_HI)

def pkas(st, wd, tag):
    p = os.path.join(wd, tag + ".pdb")
    st.write_pdb(p)
    subprocess.run([sys.executable, "-m", "propka", tag + ".pdb"], capture_output=True, cwd=wd)
    f = os.path.join(wd, tag + ".pka")
    out = {}
    if os.path.exists(f):
        for ln in open(f):
            q = ln.split()
            if len(q) > 3 and q[0] == "HIS" and q[2] == "B":
                try: out[int(q[1])] = float(q[3])
                except ValueError: pass
    return out

def one(cif):
    with tempfile.TemporaryDirectory() as wd:
        st = gemmi.read_structure(cif); st.setup_entities()
        st.remove_ligands_and_waters(); st.setup_entities()
        bound = pkas(st, wd, "cpx")
        m = st[0]                                  # same file, binder removed
        for i in range(len(m) - 1, -1, -1):
            if m[i].name != "B": del m[i]
        st.setup_entities()
        free = pkas(st, wd, "apo")
    return os.path.basename(cif)[:-4], free, bound

if __name__ == "__main__":
    jobs = []
    for a in ["g532mimic", "g532mimic_short", "g532mimic_tiny", "tinyHis"]:
        for f in sorted(glob.glob(f"runs/g532mimic/gm_{a}/final_ranked_designs/final_30_designs/rank*.cif")):
            jobs.append((a, f))
    rows = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for (a, _), (name, free, bound) in zip(jobs, ex.map(one, [f for _, f in jobs])):
            for num, lbl in CANON.items():
                if num not in free or num not in bound: continue
                r = ratio(free[num], bound[num])
                rows.append(dict(arm=a, design=name, his=lbl,
                                 pka_free=round(free[num], 2), pka_bound=round(bound[num], 2),
                                 delta=round(bound[num] - free[num], 2), ratio=round(r, 3),
                                 passes=r >= RATIO_BAR))
            print(f"  {name}", file=sys.stderr)
    with open("analysis/01-egfr/phgate_g532mimic_paired.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} rows")
