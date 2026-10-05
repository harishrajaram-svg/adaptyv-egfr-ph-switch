#!/usr/bin/env python3
"""
Known-answer control for bin/dddg_elec.py, on OUR TARGET, with MEASURED outcomes.

Why this exists
    dddG_elec is the selection filter for problem 2's first-ranked objective. playbook s8/s10:
    an instrument that has not been run against a known answer is not an instrument. Problem 1's
    most valuable single result was Boltz-2 FAILING its validation gate -- it built, ran, emitted
    a PAE matrix and produced numbers, and only a known answer exposed that it scored one of the
    tightest complexes known identically to random shuffles.

The known answer
    Schroeter et al., "A generic approach to engineer antibody pH-switches using combinatorial
    histidine scanning libraries and yeast display", mAbs 7(1):138-151 (2015).
    PMID 25523975, PMC4622719, doi 10.4161/19420862.2014.985993.  Retrieved via PubMed.
    Combinatorial histidine libraries into adalimumab CDRs, selected against trimeric human TNF
    -- our target, our switch direction (hold at 7.4, release at 6.0).

    Measured kd(pH 6.0)/kd(pH 7.4), i.e. fold-enhanced antigen release:
        adalimumab WT   9x      (not 1x -- s17: know what your method returns for "nothing")
        PSV#1         231x      VH S100bH S100cH        VL R90H N92H T97H
        PSV#2         785x      VH S100bH S100cH        VL Q89H R90H N92H
        PSV#3         505x      VH L98H                 VL Y32H L33H
    Kabat numbering, as the paper states for all residue positions.

PRE-REGISTERED PREDICTIONS -- written into the plan file before this script was run (s13)
    P1  (decisive) all three PSVs score dddG_elec > 0 AND above WT. If they do not separate
        from WT, the instrument does not work and nothing downstream may use it.
    P2  WT is not zero. Its measured floor is 9x. The paper attributes that partly to TNF His73
        -- which is OUR H149 -- plus adalimumab's own CDR His35/His56.
    P3  the rank order WITHIN the PSVs will probably NOT reproduce, and that is expected rather
        than a failure: the paper attributes PSV#1/#2 largely to INTRAmolecular paratope
        distortion and only PSV#3 to INTERmolecular repulsion ("the side chains of Leu-98 and
        Tyr-32 ... project deep into binding pockets on TNF"). dddG_elec is interface-only by
        construction, so it should rank PSV#3 highest while measurement ranks PSV#2 highest.
        This bounds what the filter can ever see and is stated before looking.
    P4  the scramble negative collapses. Same number of histidines, placed away from the
        interface: a score that does not care where they go is measuring composition.

Self-tests (abort before any result is printed)
    T1  every Kabat position maps to a 3WD5 residue of the expected IDENTITY. This is the check
        that catches a bad Kabat->author offset, which would otherwise mutate the wrong residues
        and produce a clean, plausible, wrong table. Determined by identity matching, not assumed.
    T2  each built mutant differs from WT at exactly the intended positions and nowhere else.
    T3  the published TNF histidines appear where our +76 numbering says (His15/73/78 = H91/H149/
        H154), independently reproducing s3's offset from this structure.

Run it:
    .venv/bin/python analysis/02-tnf/schroter_control.py
Depends: bin/dddg_elec.py, pdbfixer, openmm, numpy. 3WD5.pdb is fetched if absent.
"""

import os
import random
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "bin"))

import dddg_elec as dd  # noqa: E402

STRUCT = os.path.join(HERE, "structures")
WT = os.path.join(STRUCT, "3WD5.pdb")
BINDER, TARGET = ("H", "L"), ("A",)

# ---- Kabat -> 3WD5 author numbering -------------------------------------------------------
# 3WD5's LIGHT chain is Kabat-numbered already (offset 0, all six positions verified by identity).
# Its HEAVY chain is SEQUENTIAL through CDR-H3, which in Kabat carries insertion codes. The
# adalimumab H3 reads VSYLSTASSLDY at author 99-110, so Kabat 95..102 with four insertions maps
# 95->99 ... 98->102, 100a->105, 100b->106, 100c->107, 100d->108, 101->109, 102->110, i.e. +4.
# The offset is DETERMINED by requiring the expected residue identity at every substituted
# position (T1) -- exactly how the +76 TNF offset was fixed in third_site_census.py, and the
# reason a hardcoded guess is not acceptable here. CDR-H1 is offset 0 (Kabat H35 = author 35)
# and CDR-H2 is +1 (the paper's His56 is author His57), both confirmed against the two
# wild-type CDR histidines the paper names.
KABAT = {
    ("H", "S100b"): ("H", 106, "SER"),
    ("H", "S100c"): ("H", 107, "SER"),
    ("H", "L98"):   ("H", 102, "LEU"),
    ("L", "Q89"):   ("L", 89,  "GLN"),
    ("L", "R90"):   ("L", 90,  "ARG"),
    ("L", "N92"):   ("L", 92,  "ASN"),
    ("L", "T97"):   ("L", 97,  "THR"),
    ("L", "Y32"):   ("L", 32,  "TYR"),
    ("L", "L33"):   ("L", 33,  "LEU"),
}

VARIANTS = {
    "adalimumab_WT": ([], 9.0),
    "PSV1": ([("H", "S100b"), ("H", "S100c"), ("L", "R90"), ("L", "N92"), ("L", "T97")], 231.0),
    "PSV2": ([("H", "S100b"), ("H", "S100c"), ("L", "Q89"), ("L", "R90"), ("L", "N92")], 785.0),
    "PSV3": ([("H", "L98"), ("L", "Y32"), ("L", "L33")], 505.0),
}

# the paper's two wild-type CDR histidines, used to fix the heavy-chain offsets (T1)
WT_CDR_HIS = {("H", 35): "CDR-H1, Kabat H35", ("H", 57): "CDR-H2, Kabat H56"}
TNF_HIS = {15: "H91", 73: "H149", 78: "H154"}   # +76 -> our numbering (s3)


def fetch(path, pdb_id):
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"fetching {url}")
        urllib.request.urlretrieve(url, path)
    return path


def residues(pdb_path):
    """{(chain, resseq): resname} from ATOM records. gemmi would also do; stdlib is enough here."""
    out = {}
    for line in open(pdb_path):
        if line.startswith("ATOM"):
            out[(line[21], int(line[22:26]))] = line[17:20].strip()
    return out


def build(name, subs, out_path):
    """Write a mutant PDB. Mutations are applied with PDBFixer, per chain."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename=WT)
    fixer.missingResidues = {}
    by_chain = {}
    for key in subs:
        ch, num, old = KABAT[key]
        by_chain.setdefault(ch, []).append(f"{old}-{num}-HIS")
    for ch, muts in by_chain.items():
        fixer.applyMutations(muts, ch)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with open(out_path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    return out_path


def scramble(n, seed=0):
    """P4 negative: n histidines at random positions FAR from the interface.

    Composition-matched in the only sense that matters for an interface-electrostatics score:
    the same number of histidines, placed where they cannot touch the other chain. The
    ranking-instrument gate used composition-matched sequence shuffles for the same purpose --
    "a scorer returning zero for everything would look identical on the negatives alone".
    """
    atoms, xyz, qp, qn = dd.prepare(WT, BINDER + TARGET, "HIE")
    iface = dd.interface_his(atoms, xyz, BINDER, TARGET)
    import numpy as np
    ch = np.array([a[0] for a in atoms])
    tmask = np.isin(ch, list(TARGET))
    # exclude chain termini: a terminal histidine is a different ff14SB residue (it carries OXT
    # and its own charge set), which the +1/0 net-charge guard in prepare() correctly refuses.
    # The first scramble draw picked ILE219, chain H's C-terminus, and the guard fired.
    ends = set()
    for c in BINDER:
        nums = sorted({int(a[1]) for a in atoms if a[0] == c})
        ends |= {(c, str(nums[0])), (c, str(nums[-1]))}
    cand = []
    for (c, r, rn, an) in atoms:
        if c not in BINDER or an != "CB" or rn in ("GLY", "PRO", "CYS") or (c, str(r)) in ends:
            continue
        sel = np.array([a[0] == c and a[1] == r for a in atoms])
        d = np.linalg.norm(xyz[sel][:, None, :] - xyz[tmask][None, :, :], axis=-1).min()
        if d > 12.0:                                   # safely outside any interface
            cand.append((c, int(r), rn))
    random.Random(seed).shuffle(cand)
    picked = cand[:n]
    return picked, len(cand), len(iface)


def main():
    fetch(WT, "3WD5")
    wt_res = residues(WT)

    # ---- self-tests -------------------------------------------------------
    for key, (ch, num, want) in KABAT.items():
        got = wt_res.get((ch, num))
        assert got == want, (f"T1 FAIL: Kabat {key[0]}:{key[1]} -> {ch}:{num} is {got}, "
                             f"paper requires {want}. The Kabat->author offset is wrong.")
    for (ch, num), where in WT_CDR_HIS.items():
        assert wt_res.get((ch, num)) == "HIS", \
            f"T1 FAIL: the paper's wild-type {where} is not HIS at {ch}:{num}"
    for num, ours in TNF_HIS.items():
        assert wt_res.get(("A", num)) == "HIS", \
            f"T3 FAIL: TNF His{num} ({ours}) absent -- the +76 offset does not reproduce here"

    built = {}
    for name, (subs, _) in VARIANTS.items():
        path = os.path.join(STRUCT, f"3WD5_{name}.pdb")
        built[name] = build(name, subs, path)
        got = residues(path)
        changed = {k for k in wt_res if wt_res[k] != got.get(k)}
        want = {KABAT[s][:2] for s in subs}
        assert changed == want, (f"T2 FAIL: {name} changed {sorted(changed)}, "
                                 f"intended {sorted(want)}")
        assert all(got[k] == "HIS" for k in want), f"T2 FAIL: {name} did not install HIS"

    n_scr = len(VARIANTS["PSV1"][0])
    picked, n_cand, n_iface_wt = scramble(n_scr)
    assert len(picked) == n_scr, f"T2 FAIL: only {len(picked)} scramble sites of {n_scr}"
    for (c, r, rn) in picked:
        KABAT[("scr", f"{rn}{r}")] = (c, r, rn)
    scr_path = os.path.join(STRUCT, "3WD5_scramble.pdb")
    build("scramble", [("scr", f"{rn}{r}") for (c, r, rn) in picked], scr_path)
    built["scramble_neg"] = scr_path

    print(f"self-tests T1 T2 T3 passed ({len(KABAT) - n_scr} Kabat positions identity-checked, "
          f"{len(VARIANTS)} variants built, {n_scr} scramble sites from {n_cand} candidates "
          f">12 A from the interface, WT interface histidines={n_iface_wt})\n")

    # ---- the control, across the s16 band ---------------------------------
    legs = [("sigmoid", "HIE"), ("sigmoid", "HID"), ("10", "HIE"), ("6", "HIE"), ("80", "HIE")]
    rows = {}
    for name, path in built.items():
        rows[name] = {}
        for eps, taut in legs:
            r = dd.score_pose(path, BINDER, TARGET, eps, taut)
            rows[name][(eps, taut)] = r
            if abs(r["decomposition_residual"]) > 1e-6:
                raise SystemExit(f"G2: decomposition residual {r['decomposition_residual']:.2e} "
                                 f"on {name} -- the per-site sum is not the score")

    hdr = f"{'variant':<16}{'measured':>10}  " + "".join(
        f"{e + '/' + t:>14}" for e, t in legs) + f"{'n_his':>7}"
    print(hdr)
    print("-" * len(hdr))
    order = ["adalimumab_WT", "PSV1", "PSV3", "PSV2", "scramble_neg"]
    for name in order:
        meas = VARIANTS.get(name, ([], None))[1]
        m = f"{meas:>9.0f}x" if meas else "       neg"
        cells = "".join(f"{rows[name][k]['dddg_elec']:>+14.3f}" for k in [(e, t) for e, t in legs])
        n = rows[name][legs[0]]["n_his_interface"]
        print(f"{name:<16}{m}  {cells}{n:>7}")

    print("\nPer-site decomposition, primary leg (sigmoid/HIE) -- G2, all sites, never a subset:")
    for name in order:
        r = rows[name][("sigmoid", "HIE")]
        print(f"  {name}")
        if not r["per_site"]:
            print("    (no interface histidine)")
        for (c, rr), s in sorted(r["per_site"].items()):
            cc = (f"{s['counter_charge']} @ {s['counter_charge_dist']:.2f} A"
                  if s["counter_charge_dist"] is not None else "none in range")
            print(f"    {c}:HIS{rr:<5} {s['side']:<7} dddG {s['dddg']:+8.3f}  "
                  f"iface {s['min_dist_to_other_side']:.2f} A  cc {cc}"
                  f"{'  [' + s['flag'] + ']' if s['flag'] else ''}")

    # ---- the pre-registered predictions, scored ---------------------------
    print("\nPRE-REGISTERED PREDICTIONS, scored on the primary leg (sigmoid/HIE):")
    prim = {k: rows[k][("sigmoid", "HIE")]["dddg_elec"] for k in built}
    wt = prim["adalimumab_WT"]
    psv = {k: prim[k] for k in ("PSV1", "PSV2", "PSV3")}
    p1a = all(v > 0 for v in psv.values())
    p1b = all(v > wt for v in psv.values())
    print(f"  P1 all PSVs > 0        : {'PASS' if p1a else 'FAIL'}  "
          f"({', '.join(f'{k} {v:+.3f}' for k, v in psv.items())})")
    print(f"  P1 all PSVs > WT       : {'PASS' if p1b else 'FAIL'}  (WT {wt:+.3f})")
    print(f"  P2 WT is not zero      : WT = {wt:+.3f} on a measured 9x floor")
    best = max(psv, key=psv.get)
    print(f"  P3 PSV3 ranks highest  : {'as predicted' if best == 'PSV3' else 'NO -- ' + best}"
          f"   (measurement ranks PSV2 highest)")
    print(f"  P4 scramble collapses  : scramble {prim['scramble_neg']:+.3f} vs "
          f"best PSV {max(psv.values()):+.3f}")
    # s14: a filter that passes everything is not a filter
    npass = sum(1 for k in built if rows[k][("sigmoid", "HIE")]["verdict"] == "pass")
    print(f"  s14 pass rate          : {npass}/{len(built)} clear dddG_elec >= 0")

    # ---- why it fails: is our mechanism even present at this epitope? -----
    import json
    pp = json.load(open(os.path.join(HERE, "per_partner.json")))
    epi = sorted(set(pp["3WD5"]["H"]["total"]) | set(pp["3WD5"]["L"]["total"]))
    anchors = {108: "R108 (anchor)", 166: "K166 (anchor)", 149: "H149", 141: "K141 (peripheral)",
               107: "R107 (mouse Gln)", 204: "K204 (peripheral)", 158: "R158 (peripheral)"}
    print("\nMECHANISM AVAILABILITY at the control's epitope (per_partner.json, 4.5 A):")
    print(f"  adalimumab contacts on TNF (UniProt numbering): {epi}")
    for num, label in sorted(anchors.items()):
        print(f"    {label:<22} {'PRESENT' if num in epi else 'absent'}")
    print("  s4's design mechanism is a binder histidine against R108 and K166. Neither is in")
    print("  this epitope, so the control cannot test the mechanism the design is built on --")
    print("  it tests whether dddG_elec recovers a measured switch driven by something else.")

    print("\nBand across all five legs (s16 -- if the band spans the order, the order is not "
          "established):")
    for name in order:
        vals = [rows[name][k]["dddg_elec"] for k in [(e, t) for e, t in legs]]
        print(f"  {name:<16} {min(vals):+.3f} .. {max(vals):+.3f}   "
              f"sign {'stable' if min(vals) * max(vals) > 0 else 'FLIPS'}")


if __name__ == "__main__":
    main()
