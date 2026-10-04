#!/usr/bin/env python3
"""Pose check for the G532 ladder -- run BEFORE any pKa calculation.

PK, 2026-10-03: "A failure to recover this control could arise from the predicted pose or
the protonation model. Assess those separately before using the pH gate to discard
candidates." This script is the pose half. If the predicted complex does not put the
antibody carboxylates near the target histidines, a pH-gate miss says nothing about PROPKA.

The published mechanism (Liu et al., Mol Ther Oncolytics 2022, PMID 36458200) is a
single-light-chain bidentate arrangement, all canonical numbering:
    LCDR1 E32        <-> EGFR H433   (mature H409)
    LCDR2 D52 / D53  <-> EGFR H370   (mature H346)
There is NO solved complex; the paper's geometry is manual docking, and the authors say so.
So this checks whether ESMFold2 independently reproduces the arrangement.

Residues are located by SEQUENCE MOTIF, never by residue number. Three numbering schemes
are live on this target (canonical, mature, renumbered) and conflating them has already
produced one false result in this project.

Chains: A = EGFR domain III (mature 311-514), B = VH, C = VL.

Usage: g532_pose_check.py <dir-of-folded-cif-or-pdb>
"""
import glob, math, sys
from pathlib import Path

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
ACID = {('ASP','OD1'),('ASP','OD2'),('GLU','OE1'),('GLU','OE2')}
RING = {'ND1','CD2','CE1','NE2','CG'}
# motif -> (offset of the residue of interest within the motif, label)
LCDR1 = "GGNNIGE"       # offset 7 is Kabat L31, offset 8 is L32.
                        # Must NOT include the K: G532Ctrl carries K31N, so a "GGNNIGEK"
                        # motif silently skips that row -- the one row the ladder needs most.
LCDR2 = "YD"            # LCDR2 starts YD; the two following positions are L52/L53
# EGFR domain III motifs around each histidine, from the mature 311-514 sequence
H370_MOTIF = "ISGDLHILP"   # the H is at offset 5  -> mature H346 = canonical H370
H433_MOTIF = "RTKQHGQF"    # the H is at offset 4  -> mature H409 = canonical H433


def load(p):
    ch = {}
    txt = Path(p).read_text().splitlines()
    if p.endswith(".cif"):
        import subprocess, tempfile, os
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        subprocess.run([sys.executable, "bin/cif2pdb.py", p, tmp],
                       capture_output=True, check=False)
        txt = Path(tmp).read_text().splitlines()
    for l in txt:
        if l.startswith("ATOM"):
            c, n = l[21], int(l[22:26])
            ch.setdefault(c, {}).setdefault(n, []).append(
                (l[12:16].strip(), l[17:20].strip(),
                 (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return ch


def seq_of(res):
    return "".join(AA3.get(res[n][0][1], "X") for n in sorted(res))


def find(res, motif, off):
    s = seq_of(res)
    i = s.find(motif)
    if i < 0:
        return None
    return sorted(res)[i + off]


def atoms(res, n, names):
    return [a[2] for a in res.get(n, []) if a[0] in names]


def mind(a, b):
    return min((math.dist(p, q) for p in a for q in b), default=float("inf"))


def main(d):
    files = sorted(glob.glob(f"{d}/**/*.cif", recursive=True)) or \
            sorted(glob.glob(f"{d}/**/*.pdb", recursive=True))
    if not files:
        sys.exit(f"no structures under {d}")
    print(f"{'structure':<34}{'E32-H433':>10}{'D52-H370':>10}{'D53-H370':>10}  verdict")
    for f in files:
        ch = load(f)
        A, C = ch.get("A"), ch.get("C")
        if not A or not C:
            print(f"{Path(f).stem[:33]:<34}  chains A/C missing ({sorted(ch)})"); continue
        h433 = find(A, H433_MOTIF, 4)
        h370 = find(A, H370_MOTIF, 5)
        l31 = find(C, LCDR1, 7)   # Kabat L31
        l32 = l31 + 1 if l31 else None
        yd = find(C, LCDR2, 0)
        l52, l53 = (yd + 2, yd + 3) if yd else (None, None)
        if None in (h433, h370, l32, l52):
            print(f"{Path(f).stem[:33]:<34}  motif not found "
                  f"(h433={h433} h370={h370} l32={l32} l52={l52})"); continue
        d1 = mind(atoms(C, l32, {o for _, o in ACID}), atoms(A, h433, RING))
        d2 = mind(atoms(C, l52, {o for _, o in ACID}), atoms(A, h370, RING))
        d3 = mind(atoms(C, l53, {o for _, o in ACID}), atoms(A, h370, RING))
        # A "-" means the position simply is not a carboxylate in THIS variant -- by design,
        # for G532V/G532Ctrl (E32->H, D52/D53->S/T) and for G5V2 at 32 (Y32). Reporting that
        # as "pose failed" conflates two different things and makes the table misleading.
        present = [x for x in (d1, d2, d3) if x < 90]
        if not present:
            v = "n/a -- no carboxylate at either proposed site in this variant (expected)"
        else:
            best = min(present)
            v = ("POSE OK -- a carboxylate within 4.0 A of a target His" if best <= 4.0
                 else "pose MARGINAL (4-8 A)" if best <= 8.0
                 else "POSE WRONG -- nearest carboxylate %.1f A from its target His; the "
                      "published salt bridge is not reproduced, so a gate result on this "
                      "structure is uninformative" % best)
        fmt = lambda x: f"{x:.2f}" if x < 90 else "  -"
        print(f"{Path(f).stem[:33]:<34}{fmt(d1):>10}{fmt(d2):>10}{fmt(d3):>10}  {v}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "runs/g532")
