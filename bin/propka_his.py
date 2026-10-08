#!/usr/bin/env python3
"""Structure-based histidine pKa for problem 2's designs. A SCREEN, NOT THE pH GATE.

🔴 READ THIS FIRST. This file reports the pKa of each histidine in ONE structure and how
much it titrates between the two assay pH values. That is NOT the pH-selectivity quantity. The
mechanism depends on the pKa SHIFT ON BINDING and on its DIRECTION:

    K(pH) = (1 + 10^(pKa_bound - pH)) / (1 + 10^(pKa_free - pH)),   ratio = K(6.0)/K(7.4)

Problem 2 binds at 7.4 and must be silent at 6.0, so protonation must WEAKEN the complex:
pKa_bound < pKa_free, ratio < 1. The opposite sign designs for low-pH binding, which is
problem 1's brief and this problem's antithesis. Nothing in this file can tell those apart.

bin/ph_gate_multisite.py is the gate. It computes both legs by deleting each partner in place,
composes over ALL sites rather than the helpful ones, and carries five guards written after an
earlier version manufactured a false claim by selecting on the outcome. Use it with --problem 2.

Keep this file for what it is good for: a cheap first look at whether a design has any histidine
whose local environment lets it titrate in the assayed range at all. If it does not, there is no
mechanism for the gate to score.

WHY. express_qc_p2.py computes net charge from a Bjellqvist/EMBOSS table whose His value is
5.98 -- the MODEL COMPOUND pKa, i.e. the free amino acid in water. Inside a folded protein the
value shifts: this project's own pre-registration records histidines losing 1.5-2.5 pKa units on
burial (outbox/PREREGISTRATION.md:165). A design's pH mechanism only works if its histidines
actually titrate BETWEEN the two assay values, 6.0 and 7.4. A His whose in-situ pKa is 4.2 is
decoration; one at 6.5 is a switch.

PROPKA 3.5.1 is pre-registered (PREREGISTRATION.md:51) and was used on problem 1 (bound pKa
8.89/9.11/8.98/8.98, HANDOFF.md:468). It was never applied to problem 2's designs. This does that.

HONEST LIMIT, from our own pre-registration (:210): the binder is unrelaxed and PROPKA on a
buried histidine is its hardest case. So a single value is a screen, not a measurement; what is
reportable is the DISTRIBUTION and how many histidines land inside the assayed interval.

    propka_his.py --selftest
    propka_his.py --dir runs/p2-esmfold2-out [--chain D]
"""
import glob
import os
import subprocess
import sys
import tempfile

ASSAY_LO, ASSAY_HI = 6.0, 7.4
MODEL_COMPOUND_HIS = 5.98        # what express_qc_p2.py assumes

# WHY NOT "IS THE pKa INSIDE 6.0-7.4". That bracket is the wrong test, and writing it exposed
# the reason: 5.98 is BELOW 6.0, so the model-compound histidine -- the canonical pH-switch
# residue, the whole reason this mechanism uses His -- would be classified as not titrating.
# What matters is how much the protonation state actually CHANGES between the two assay values:
#
#     f(pH) = 1 / (1 + 10^(pH - pKa))          fraction protonated
#     delta = f(6.0) - f(7.4)
#
# That is exact arithmetic once the pKa is known, it peaks at pKa 6.7 (the interval's midpoint)
# where delta = 0.667, and it falls off smoothly in both directions instead of at a cliff.
# Thresholds are declared HERE, before any design is scored (playbook s13):
SWITCH_DELTA = 0.25             # a real protonation switch across the assayed interval
WEAK_DELTA = 0.10               # measurable but small
MAX_DELTA = 0.667               # the ceiling, at pKa 6.7


def frac_protonated(pka, ph):
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


def delta_protonation(pka):
    """How much of the histidine flips between pH 6.0 and pH 7.4. Exact, given the pKa."""
    return frac_protonated(pka, ASSAY_LO) - frac_protonated(pka, ASSAY_HI)


def classify(pka):
    d = delta_protonation(pka)
    return ("switch" if d >= SWITCH_DELTA else "weak" if d >= WEAK_DELTA else "inert"), d


def cif_to_pdb(cif_path, out_pdb):
    """gemmi, because PROPKA 3.5 takes PDB and ESMFold2 writes mmCIF."""
    import gemmi

    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    st.write_pdb(out_pdb)
    return out_pdb


def run_propka(pdb_path, workdir):
    """[(chain, resnum, pKa)] for every HIS PROPKA reports."""
    r = subprocess.run([sys.executable, "-m", "propka", os.path.basename(pdb_path)],
                       cwd=workdir, capture_output=True, text=True)
    pka_file = os.path.join(workdir, os.path.basename(pdb_path)[:-4] + ".pka")
    if not os.path.exists(pka_file):
        return None, (r.stderr or r.stdout)[-400:]
    out, in_summary = [], False
    for line in open(pka_file):
        if line.startswith("SUMMARY OF THIS PREDICTION"):
            in_summary = True
            continue
        if in_summary:
            if line.startswith("-") or not line.strip():
                if out:
                    break
                continue
            f = line.split()
            if len(f) >= 4 and f[0] == "HIS":
                try:
                    out.append((f[2], int(f[1]), float(f[3])))
                except ValueError:
                    pass
    return out, None


def summarise(rows):
    """rows = [(design, chain, resnum, pka)] -> what the mechanism actually depends on."""
    import statistics

    if not rows:
        return {}
    pkas = [r[3] for r in rows]
    deltas = [delta_protonation(p) for p in pkas]
    sw = [r for r in rows if classify(r[3])[0] == "switch"]
    wk = [r for r in rows if classify(r[3])[0] == "weak"]
    inert = [r for r in rows if classify(r[3])[0] == "inert"]
    return {
        "n_his": len(rows),
        "n_designs": len({r[0] for r in rows}),
        "pka_min": min(pkas), "pka_median": statistics.median(pkas), "pka_max": max(pkas),
        "delta_max_observed": max(deltas), "delta_median": statistics.median(deltas),
        "his_switch": len(sw), "his_weak": len(wk), "his_inert": len(inert),
        "designs_with_a_switch_his": len({r[0] for r in sw}),
    }


def selftest():
    import gemmi  # noqa: F401
    import propka  # noqa: F401

    # T1 -- the model-compound His is BELOW the lower assay bound, which is why the bracket
    # test is wrong; it is nonetheless a strong switch by the delta metric.
    assert MODEL_COMPOUND_HIS < ASSAY_LO, "5.98 is below 6.0; that is the point"
    d598 = delta_protonation(MODEL_COMPOUND_HIS)
    assert d598 >= SWITCH_DELTA, d598
    # T2 -- the ceiling really is at the interval midpoint
    best = max((delta_protonation(p), p) for p in [x / 100 for x in range(400, 1000)])
    assert abs(best[1] - 6.7) < 0.06, best
    assert abs(best[0] - MAX_DELTA) < 0.01, best
    # T3 -- monotone fall-off either side of the midpoint
    assert delta_protonation(6.7) > delta_protonation(5.5) > delta_protonation(4.0)
    assert delta_protonation(6.7) > delta_protonation(8.0) > delta_protonation(10.0)
    # T4 -- classification on a known spread
    rows = [("d1", "D", 10, 4.2), ("d1", "D", 20, 6.5), ("d2", "D", 5, 8.1),
            ("d2", "D", 9, 7.0), ("d3", "D", 1, 10.0)]
    s = summarise(rows)
    assert s["n_his"] == 5 and s["n_designs"] == 3, s
    assert s["his_switch"] == 2, s          # 6.5 and 7.0
    assert s["his_inert"] == 2, s           # 4.2 and 10.0
    assert s["his_weak"] == 1, s            # 8.1
    assert s["designs_with_a_switch_his"] == 2, s
    # T5 -- MUTATION: a buried His is the failure this exists to catch. Our own
    # pre-registration records 1.5-2.5 pKa units lost on burial, so shift 6.5 down by 2.0
    # and it must stop being a switch.
    assert classify(6.5)[0] == "switch"
    assert classify(6.5 - 2.0)[0] == "inert", classify(4.5)
    # T6 -- MUTATION: the delta must be computed across BOTH assay values. Using one pH
    # alone cannot distinguish a switch from a permanently charged residue.
    assert frac_protonated(10.0, ASSAY_LO) > 0.99 and frac_protonated(10.0, ASSAY_HI) > 0.99
    assert delta_protonation(10.0) < 0.01
    print(f"  ok  propka 3.5.1 and gemmi import")
    print(f"  ok  model-compound His {MODEL_COMPOUND_HIS} is BELOW {ASSAY_LO} yet switches "
          f"(delta {d598:.3f}) -- the bracket test would have called it inert")
    print(f"  ok  delta peaks at pKa {best[1]:.2f} with {best[0]:.3f}, falls off both sides")
    print(f"  ok  classification at thresholds {SWITCH_DELTA}/{WEAK_DELTA}")
    print(f"  ok  MUTATION: a 2.0-unit burial shift turns a switch (6.5) inert (4.5)")
    print(f"  ok  MUTATION: pKa 10 is >99% protonated at BOTH values, delta <0.01")
    print("\nself-tests passed: 6")


def main(d, chain=None):
    cifs = sorted(glob.glob(os.path.join(d, "**", "*.cif"), recursive=True))
    if not cifs:
        sys.exit(f"no .cif under {d}")
    rows, failed = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for cif in cifs:
            name = os.path.basename(cif).split("_seed")[0]
            pdb = os.path.join(tmp, f"{name}.pdb")
            try:
                cif_to_pdb(cif, pdb)
            except Exception as e:
                failed.append((name, f"cif->pdb: {e}")); continue
            his, err = run_propka(pdb, tmp)
            if his is None:
                failed.append((name, err)); continue
            for ch, num, pka in his:
                if chain and ch != chain:
                    continue
                rows.append((name, ch, num, pka))
    print(f"{len(cifs)} structures, {len(failed)} failed, "
          f"{len(rows)} histidines scored" + (f" (chain {chain} only)" if chain else ""))
    for n, e in failed[:5]:
        print(f"  FAILED {n}: {str(e)[:120]}")
    s = summarise(rows)
    if not s:
        sys.exit("no histidines scored")
    print()
    for k, v in s.items():
        print(f"  {k:32}{v if isinstance(v, int) else round(v, 2)}")
    print(f"\n  express_qc_p2.py assumes every His sits at {MODEL_COMPOUND_HIS} "
          f"(delta {delta_protonation(MODEL_COMPOUND_HIS):.3f}); ceiling is {MAX_DELTA} at 6.7")
    print("\ndesign\tchain\tresnum\tpka\tdelta_protonation\tclass")
    for n, ch, num, pka in sorted(rows, key=lambda r: -delta_protonation(r[3])):
        cls, d = classify(pka)
        print(f"{n}\t{ch}\t{num}\t{pka:.2f}\t{d:.3f}\t{cls}")


if __name__ == "__main__":
    if len(sys.argv) == 1 or "--selftest" in sys.argv:
        selftest()
    else:
        a = sys.argv
        main(a[a.index("--dir") + 1], a[a.index("--chain") + 1] if "--chain" in a else None)
