#!/usr/bin/env python3
"""Expression / developability QC for Adaptyv submissions. Sequence-only, no dependencies.

WHY THIS MATTERS AND WHY IT WAS MISSING
---------------------------------------
Adaptyv expresses designs by CELL-FREE protein synthesis and screens ~375 designs per
problem. A design that does not express cannot be tested, whatever its ipSAE or pH score.
Nothing in this project had ever checked it.

Cell-free-specific liabilities, in rough order of how often they bite:
  * free cysteines -- standard CFPS is reducing; an ODD cysteine count is the worst case
    (unpaired thiol, scrambling risk). Even counts can still mispair.
  * strong net positive charge -- binds ribosomes/nucleic acids in lysate, kills yield
  * long hydrophobic runs -- aggregation / inclusion-body-like precipitation
  * polybasic clusters -- nonspecific binding, heparin-like behaviour
  * Asn/Asp motifs (NG/NS/NT, DG/DP) -- deamidation and isomerization on storage

Charge is reported at BOTH pH 7.4 and pH 6.5 because the competition objective is a
pH switch: a design with more titratable charge in that window has more to work with.
"""
import csv, os
from pathlib import Path, sys

KD = dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,
          L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
PKA = dict(Cterm=3.65, D=3.9, E=4.07, C=8.18, Y=10.46, H=6.04, Nterm=8.2, K=10.54, R=12.48)
NEG, POS = "DECY", "HKR"

def charge(seq, pH):
    z = 0.0
    z += 1/(1+10**(pH-PKA["Nterm"]))
    z -= 1/(1+10**(PKA["Cterm"]-pH))
    for aa in seq:
        if aa in POS: z += 1/(1+10**(pH-PKA[aa]))
        elif aa in NEG: z -= 1/(1+10**(PKA[aa]-pH))
    return z

def pI(seq):
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo+hi)/2
        if charge(seq, mid) > 0: lo = mid
        else: hi = mid
    return (lo+hi)/2

def max_hydrophobic_run(seq, thresh=1.5):
    best = cur = 0
    for aa in seq:
        cur = cur+1 if KD.get(aa,0) >= thresh else 0
        best = max(best, cur)
    return best

def polybasic(seq, win=6, need=4):
    return sum(1 for i in range(len(seq)-win+1)
               if sum(c in "KR" for c in seq[i:i+win]) >= need)

def motifs(seq):
    d = {m: 0 for m in ("NG","NS","NT","DG","DP")}
    for i in range(len(seq)-1):
        p = seq[i:i+2]
        if p in d: d[p] += 1
    sequons = sum(1 for i in range(len(seq)-2)
                  if seq[i]=="N" and seq[i+1]!="P" and seq[i+2] in "ST")
    return d, sequons

def assess(name, seq):
    cys = seq.count("C")
    dm, seq_n = motifs(seq)
    z74, z65 = charge(seq,7.4), charge(seq,6.5)
    flags = []
    if cys % 2 == 1: flags.append("ODD-CYS")
    elif cys: flags.append(f"{cys}cys")
    if z74 > 4: flags.append("v.CATIONIC")
    elif z74 > 2: flags.append("cationic")
    if max_hydrophobic_run(seq) >= 6: flags.append("HYDROPHOBIC-RUN")
    if polybasic(seq): flags.append("POLYBASIC")
    if sum(dm.values()) >= 8: flags.append("deamid-heavy")
    return dict(name=name, aa=len(seq), cys=cys, pI=pI(seq), z74=z74, z65=z65,
                dz=z65-z74, gravy=sum(KD.get(a,0) for a in seq)/len(seq),
                hrun=max_hydrophobic_run(seq), pb=polybasic(seq),
                deamid=sum(dm.values()), sequon=seq_n,
                flags=",".join(flags) if flags else "clean")

DEFAULT_OUT = "analysis/01-egfr/express_qc.tsv"


def read_seqs(src):
    """Accept FASTA (.faa/.fa/.fasta) or CSV with name,sequence columns."""
    if src.lower().endswith((".faa", ".fa", ".fasta")):
        name, buf = None, []
        for line in open(src):
            line = line.strip()
            if line.startswith(">"):
                if name and buf: yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
        if name and buf: yield name, "".join(buf)
    else:
        for r in csv.DictReader(open(src)):
            if r.get("sequence"): yield r["name"], r["sequence"].strip()


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = DEFAULT_OUT
    for a in sys.argv[1:]:
        if a.startswith("--out="): out = a.split("=", 1)[1]
    if not args:
        raise SystemExit("usage: express_qc.py <seqs.faa|seqs.csv> ... [--out=PATH]\n"
                         f"       default --out is {DEFAULT_OUT}")

    # THE CANONICAL RECORD IS NOT A SCRATCH FILE.
    #
    # `analysis/01-egfr/express_qc.tsv` is the submission's expression-QC record and it is
    # cited in METHODS limitation 35. Until 2026-10-05 any ad-hoc run wrote there by
    # default, so clearing a swap candidate silently replaced the artifact describing the
    # shipped slate -- three separate reviewers did exactly that in one afternoon and each
    # had to restore it from git. Writing the default path now requires the submission CSV
    # as the input, which is the only input that should ever produce it.
    SUB = "submissions/01-egfr.csv"
    if os.path.abspath(out) == os.path.abspath(DEFAULT_OUT) and \
            not any(os.path.abspath(a) == os.path.abspath(SUB) for a in args):
        raise SystemExit(
            f"refusing to overwrite {DEFAULT_OUT} from {args!r}.\n"
            f"  That file is the submission's QC record, cited in METHODS limitation 35.\n"
            f"  Pass --out=PATH to score anything else, or run it on {SUB} to regenerate it.")

    rows = [assess(n, s) for src in args for n, s in read_seqs(src)]
    # Guarded because an empty input used to crash on rows[0] AFTER printing a
    # header, which read as "0/0 clean" -- a pass, not a failure (2026-10-02).
    if not rows:
        raise SystemExit("no sequences parsed -- check the input format "
                         "(FASTA, or CSV with name,sequence columns)")

    hdr = f"{'design':34} {'aa':>4} {'cys':>4} {'pI':>5} {'z7.4':>6} {'z6.5':>6} {'dz':>5} {'GRAVY':>6} {'hrun':>5} {'deam':>5}  flags"
    print(hdr); print("-"*len(hdr))
    for r in sorted(rows, key=lambda x: (x["flags"]=="clean", x["name"])):
        print(f"{r['name'][:34]:34} {r['aa']:>4} {r['cys']:>4} {r['pI']:>5.2f} "
              f"{r['z74']:>+6.1f} {r['z65']:>+6.1f} {r['dz']:>+5.1f} {r['gravy']:>+6.2f} "
              f"{r['hrun']:>5} {r['deamid']:>5}  {r['flags']}")
    n_clean = sum(1 for r in rows if r["flags"]=="clean")
    print(f"\n{n_clean}/{len(rows)} carry no flagged liability.")

    Path(out).parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(f"written: {out}")


main()
