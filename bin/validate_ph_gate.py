#!/usr/bin/env python3
"""Validate the pH gate against an experimentally measured pH-dependent interaction.

WHY (2026-10-03). The gate -- PROPKA3 plus the thermodynamic linkage below -- decided the
entire problem-1 result, and it had NEVER been run against a real measurement. Everything
reported as "no design has a defensible pH switch" rests on an instrument nobody validated.
Lesson 5 of lessons-problem-1 says to name the instrument for every ranked objective; that
was applied to the BINDING instrument and never to this one.

POSITIVE CONTROL: IgG Fc / FcRn, the canonical histidine-driven pH switch. Acid STRENGTHENS
binding (the receptor grabs IgG in the acidified endosome and releases it at serum pH).
  1I1A  rat FcRn + rat IgG2a Fc, 2.8 A. Chains A=FcRn alpha, B=b2m, C=Fc (functional half),
        D=Fc (deliberately disabled triple mutant). Fc chain C carries His310/433/435/436.
  6WNA  human FcRn + human IgG1 Fc (YTE), 2.4 A. Chains A/B/H. Residue 436 is TYR in human,
        so the human interface has two interfacial His where the rat one has three.

Measured truth over OUR 6.5 -> 7.4 window: ~41x (interpolated from the nine-point
switchSENSE ladder in Mevada et al., mAbs 2024;16:2361585, PMID 38849969 -- K_D ~12 nM at
pH 6.5 vs 494 nM at pH 7.4). The canonical "~100x" figure is for 6.0 -> 7.0, a different
window. Direction is unambiguous in every source.

WHAT PASSING MEANS -- and it cannot mean matching 41x. Our own formula saturates: with
pKa_free = 6.5 the ratio cannot exceed 4.47 however high pKa_bound goes, and the absolute
ceiling per titratable site over a 0.9 pH-unit window is 10^0.9 = 7.94. Reproducing 41x
needs >=2 sites and realistically 3. So the gate is judged on SIGN, RESIDUE RANKING, and
order of magnitude -- never on matching the measured fold-change.

Usage:  validate_ph_gate.py <dir-with-1I1A.pdb-and-6WNA.pdb>
"""
import os, shutil, subprocess, sys, tempfile
from pathlib import Path

PH_LO, PH_HI = 6.5, 7.4


def ratio(pka_free, pka_bound, ph_lo=PH_LO, ph_hi=PH_HI):
    """Identical to bin/ph_gate.py and bin/ph_gate_mechA.py. Do not reimplement."""
    K = lambda ph: (1 + 10 ** (pka_bound - ph)) / (1 + 10 ** (pka_free - ph))
    return K(ph_lo) / K(ph_hi)


def propka(pdb_path):
    """Run PROPKA3, return {(resname, resnum, chain): pKa}. Mirrors ph_gate_mechA.py."""
    with tempfile.TemporaryDirectory() as wd:
        local = os.path.join(wd, os.path.basename(pdb_path))
        shutil.copy(pdb_path, local)
        subprocess.run([sys.executable, "-m", "propka", local],
                       capture_output=True, text=True, cwd=wd)
        pka = os.path.splitext(local)[0] + ".pka"
        if not os.path.exists(pka):
            return {}
        out, on = {}, False
        for line in Path(pka).read_text().splitlines():
            if line.startswith("SUMMARY OF THIS PREDICTION"):
                on = True; continue
            if on:
                f = line.split()
                if len(f) >= 4 and f[0].isalpha() and f[1].isdigit():
                    try: out[(f[0], int(f[1]), f[2])] = float(f[3])
                    except ValueError: pass
                elif line.startswith("-" * 10):
                    break
        return out


def subset(src, dst, chains):
    keep = [l for l in Path(src).read_text().splitlines()
            if l.startswith("ATOM") and l[21] in chains]          # HETATM stripped
    Path(dst).write_text("\n".join(keep) + "\nEND\n")
    return len(keep)


CASES = [
    dict(pdb="1I1A", role="POSITIVE CONTROL (rat)", fc="C",
         complex="ABCD", free_rec="AB", free_lig="CD",
         expect="sign >1; ranking 436 > 310 > 435 >> 433; every FcRn His exactly 0.00"),
    dict(pdb="6WNA", role="DIAGNOSTIC (human)", fc="H",
         complex="ABH", free_rec="AB", free_lig="H",
         expect="His435 is stabilised by cation-pi with Trp131, which PROPKA has no term "
                "for; a wrong-sign result here is the gate's blind spot, not a code bug"),
]
HIS = [310, 433, 435, 436]


def main(d):
    d = Path(d)
    for c in CASES:
        src = d / f"{c['pdb']}.pdb"
        if not src.exists():
            print(f"missing {src}"); continue
        print(f"\n{'='*78}\n{c['pdb']} -- {c['role']}\n  expect: {c['expect']}\n{'='*78}")
        with tempfile.TemporaryDirectory() as wd:
            wd = Path(wd)
            subset(src, wd / "cplx.pdb", set(c["complex"]))
            subset(src, wd / "lig.pdb",  set(c["free_lig"]))
            subset(src, wd / "rec.pdb",  set(c["free_rec"]))
            pc, pl = propka(wd / "cplx.pdb"), propka(wd / "lig.pdb")
            pr = propka(wd / "rec.pdb")     # the receptor leg; without this the
                                            # receptor-side check below is vacuous
        if not pc or not pl:
            print("  PROPKA produced no output"); continue
        rows, prod = [], 1.0
        for n in HIS:
            k = ("HIS", n, c["fc"])
            if k not in pc or k not in pl:
                rows.append((n, None, None, None, None)); continue
            f, b = pl[k], pc[k]
            r = ratio(f, b)
            rows.append((n, f, b, b - f, r))
            if n != 433: prod *= r
        print(f"  {'His':>5}{'pKa_free':>10}{'pKa_bound':>11}{'dpKa':>8}{'ratio':>8}")
        for n, f, b, dd, r in rows:
            if f is None: print(f"  {n:>5}{'not present / not titrated':>37}"); continue
            print(f"  {n:>5}{f:>10.2f}{b:>11.2f}{dd:>+8.2f}{r:>8.2f}")
        rec = [(k, pc[k] - pr[k]) for k in pc if k[0] == "HIS" and k[2] in c["free_rec"]
               and k in pr]
        mx = max((abs(v) for _, v in rec), default=0.0)
        print(f"\n  receptor-side His (chains {c['free_rec']}): {len(rec)} found, "
              f"largest |dpKa| = {mx:.2f}")
        print(f"  combined ratio over the interfacial His (433 excluded): {prod:.2f}x")
        print(f"  measured truth over 6.5->7.4: ~41x   |   our per-site ceiling: "
              f"{10**0.9:.2f}x, formula saturates at {ratio(6.5,30):.2f}x")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
