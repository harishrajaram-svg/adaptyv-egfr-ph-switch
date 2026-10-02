#!/usr/bin/env python3
"""pH-switch gate for Adaptyv problem 1 (EGFR, pH 6.5 vs 7.4).

The only filter in our pipeline that can see the objective the organizers rank FIRST.
No structure predictor takes pH as input, so this stands in for one.

Physics
-------
Binding and protonation are thermodynamically linked:

    K_obs(pH)  proportional to  (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH))

    pKa_bound > pKa_free  ->  protonation favours binding  ->  acid STRENGTHENS  (wanted)
    pKa_bound < pKa_free  ->  protonation opposes binding   ->  acid WEAKENS     (fatal)

Burying a histidine at an interface desolvates it and LOWERS its pKa by default, which
is why naive designs land on the wrong side. Measured on 6ARU: cetuximab drives H433
from 6.22 to 5.10, giving 0.72x -- binding 1.4x weaker at pH 6.5.

Usage
-----
    python3 bin/ph_gate.py complex.pdb [--target-chain A] [--his 409]
    python3 bin/ph_gate.py designs/*.pdb --tsv results.tsv
"""
import argparse, subprocess, sys, os, tempfile, shutil

FREE_PKA_H433 = 6.22      # PROPKA3 on 6ARU chain A alone; see analysis/01-egfr/
PH_LO, PH_HI  = 6.5, 7.4
MATURE_HIS    = 409       # canonical H433 = mature 409 (PDB numbering)


def ratio(pka_free, pka_bound, ph_lo=PH_LO, ph_hi=PH_HI):
    """K(ph_lo)/K(ph_hi). >1 = stronger in acid."""
    K = lambda ph: (1 + 10 ** (pka_bound - ph)) / (1 + 10 ** (pka_free - ph))
    return K(ph_lo) / K(ph_hi)


RATIO_BAR = 1.20          # below this the "switch" is inside PROPKA's own noise
CONTACT_CUTOFF = 5.0      # A, binder heavy atom to His heavy atom


def his_contacts(pdb_path, chain, his):
    """Binder heavy atoms within CONTACT_CUTOFF of the target histidine.

    0 means the design does not engage the switch at all -- see assess().
    """
    import gemmi
    st = gemmi.read_structure(pdb_path)
    st.setup_entities()
    st.remove_hydrogens()
    mdl = st[0]
    tgt = [c for c in mdl if c.name == chain]
    if not tgt:
        return 0
    res = [r for r in tgt[0] if r.seqid.num == his]
    if not res:
        return 0
    ns = gemmi.NeighborSearch(st, CONTACT_CUTOFF + 1.0).populate()
    n = 0
    for at in res[0]:
        for m in ns.find_atoms(at.pos, '\0', radius=CONTACT_CUTOFF):
            if m.to_cra(mdl).chain.name != chain:
                n += 1
    return n


def run_propka(pdb_path, workdir):
    """Run PROPKA3, return {(resname,resnum,chain): pKa}."""
    local = os.path.join(workdir, os.path.basename(pdb_path))
    shutil.copy(pdb_path, local)
    r = subprocess.run([sys.executable, "-m", "propka", local],
                       capture_output=True, text=True, cwd=workdir)
    pka_file = os.path.splitext(local)[0] + ".pka"
    if not os.path.exists(pka_file):
        raise RuntimeError(f"propka produced no output for {pdb_path}\n{r.stderr[-500:]}")
    out, in_summary = {}, False
    for line in open(pka_file):
        if line.startswith("SUMMARY"):
            in_summary = True; continue
        if in_summary:
            if line.startswith("-----"): break
            f = line.split()
            if len(f) >= 4:
                try: out[(f[0], int(f[1]), f[2])] = float(f[3])
                except ValueError: pass
    return out


def assess(pdb_path, chain="A", his=MATURE_HIS, free=FREE_PKA_H433):
    with tempfile.TemporaryDirectory() as wd:
        pkas = run_propka(pdb_path, wd)
    bound = pkas.get(("HIS", his, chain))
    if bound is None:
        return dict(pdb=pdb_path, error=f"HIS {his} chain {chain} not found in PROPKA output")
    r = ratio(free, bound)
    # A design that never TOUCHES the histidine leaves its pKa untouched and so satisfies
    # `bound > free` trivially, at ratio ~1.01, with no switch at all. Measured 2026-10-02:
    # 9 of 23 gate passers had zero binder atoms within 5 A of H433. Engagement is now
    # required, and the ratio has to clear a margin rather than merely exceed 1.0.
    n_contacts = his_contacts(pdb_path, chain, his)
    engaged = n_contacts > 0
    real = engaged and r >= RATIO_BAR
    return dict(pdb=pdb_path, pka_free=free, pka_bound=bound,
                delta=bound - free, ratio=r, contacts=n_contacts,
                passes=real,
                verdict=("PASS - real switch" if real else
                         ("trivial - H433 never contacted" if not engaged else
                          ("weak - ratio below %.2f" % RATIO_BAR if r > 1.0 else "FAIL - acid weakens"))))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbs", nargs="+", help="complex PDB(s): target chain + binder chain")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--his", type=int, default=MATURE_HIS, help="His residue number, PDB numbering")
    ap.add_argument("--free-pka", type=float, default=FREE_PKA_H433)
    ap.add_argument("--tsv", help="also write results here")
    a = ap.parse_args()

    rows = []
    print(f"{'design':<34} {'pKa_bound':>9} {'delta':>7} {'ratio':>7} {'cont':>5}  verdict")
    print("-" * 86)
    for p in a.pdbs:
        try:
            r = assess(p, a.target_chain, a.his, a.free_pka)
        except Exception as e:
            print(f"{os.path.basename(p):<34} {'ERROR':>9}  {e}"); continue
        if "error" in r:
            print(f"{os.path.basename(p):<34} {'n/a':>9}  {r['error']}"); continue
        rows.append(r)
        print(f"{os.path.basename(p):<34} {r['pka_bound']:>9.2f} {r['delta']:>+7.2f} "
              f"{r['ratio']:>6.2f}x {r['contacts']:>5d}  {r['verdict']}")

    if rows:
        n = sum(1 for r in rows if r["passes"])
        print("-" * 78)
        print(f"{n}/{len(rows)} REAL switches (H433 contacted AND ratio >= {RATIO_BAR}).")
        print(f"  of the rest: {sum(1 for r in rows if r['contacts'] == 0)} never contact H433, "
              f"{sum(1 for r in rows if r['contacts'] > 0 and 1.0 < r['ratio'] < RATIO_BAR)} engage but too weak, "
              f"{sum(1 for r in rows if r['ratio'] <= 1.0)} go the wrong way.")
        print(f"Reference: cetuximab scores 5.10 / 0.72x -- a FAIL, and it is a real drug.")
    if a.tsv and rows:
        with open(a.tsv, "w") as fh:
            fh.write("design\tpka_free\tpka_bound\tdelta\tratio_6.5_over_7.4\t"
                     "contacts\tpasses\tverdict\n")
            for r in rows:
                fh.write(f"{r['pdb']}\t{r['pka_free']}\t{r['pka_bound']:.2f}\t"
                         f"{r['delta']:+.2f}\t{r['ratio']:.3f}\t{r['contacts']}\t"
                         f"{r['passes']}\t{r['verdict']}\n")
        print(f"wrote {a.tsv}")

if __name__ == "__main__":
    main()
