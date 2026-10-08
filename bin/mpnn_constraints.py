#!/usr/bin/env python3
"""Per-position sequence constraints for SolubleMPNN: keep acids off the interface, put
histidines where they can reach a cation.

WHY THIS EXISTS. The direction-aware linkage gate run on 2026-10-08 found that 139 of 219
moving sites across the 35 Mosaic designs push pH-selectivity the WRONG WAY, and that the
whole-molecule product points the wrong way for 24 of 35. Problem 2 needs acid to WEAKEN
binding (ratio < 1); those designs mostly get acid TIGHTENING it.

The hypothesis -- stated as a hypothesis, not a finding -- is composition. These designs are
acid-rich, and both generators we have produce acid-rich sequences (the Genie3 seed sequences
measured 30-39% glutamate). An interface where BOTH sides carry carboxylates is mildly
self-repelling at pH 7.4; drop to 6.0, the acids protonate, the repulsion eases and the complex
gets TIGHTER. That is an inverted switch built by accident, and no amount of histidine placement
fixes it if the acids dominate.

So this emits two constraints for ProteinMPNN/SolubleMPNN:

  1. ANTI-INVERSION. At binder positions facing the target, omit ASP and GLU.
  2. MECHANISM. At binder positions that could reach a target cation, allow only HIS.

THRESHOLDS, declared before use (playbook s13), and taken from this project's own numbers or
derived from its own structures -- not invented here:

  IFACE_CA  8.00 A   his_cation_gate.py's SLACK_D, its deliberately generous "at the interface"
  REACH_CA  8.74 A   DERIVED: max CA->CB 1.55 + max CB->ring N 3.69 + D_HB 3.5, measured over
                     1539 CA-CB bonds and 222 histidines in this project's own ESMFold2 output.
                     It is the furthest a histidine placed at that CA could possibly put a ring
                     nitrogen within hydrogen-bonding distance of a cation, assuming an ideal
                     rotamer. Necessary, not sufficient -- a screen for candidate positions.

CATIONS are ARG NE/NH1/NH2 and LYS NZ, the same donor set his_cation_gate.py uses.

THE TRANSFER ASSUMPTION, STATED. Positions are computed from the Genie 3 COMPLEX, where the
target is present, but SolubleMPNN designs on the FOLDED MONOMER, where it is not. The binder is
residues 1..N in both, so the sets transfer by index. That is only sound if the fold reproduces
the generated backbone; `backbone_rmsd()` reports it and the caller must look.

    mpnn_constraints.py --selftest
    mpnn_constraints.py --complex <genie3.pdb> --binder A --out-dir <dir> [--monomer <folded>]
"""
import json
import math
import os
import sys

IFACE_CA = 8.00
REACH_CA = 8.74
D_HB = 3.5                      # his_cation_gate.py
MAX_CA_CB = 1.55                # measured, n=1539
MAX_CB_RING_N = 3.69            # measured, n=222
CATIONS = {("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ")}
ACIDS = "DE"
ALL_AA = "ACDEFGHIKLMNPQRSTVWY"


def read_atoms(path):
    """[(chain, resnum, resname, atomname, (x,y,z))] for every ATOM record."""
    out = []
    for l in open(path):
        if l.startswith("ATOM"):
            out.append((l[21], int(l[22:26]), l[17:20].strip(), l[12:16].strip(),
                        (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    return out


def positions(complex_pdb, binder_chain):
    """(interface, cation_reach) -- 1-based binder residue numbers, as ProteinMPNN wants."""
    at = read_atoms(complex_pdb)
    binder = [(n, p) for c, n, rn, an, p in at if c == binder_chain and an == "CA"]
    if not binder:
        sys.exit(f"REFUSE: no CA atoms on binder chain {binder_chain} in {complex_pdb}")
    tgt_all = [p for c, n, rn, an, p in at if c != binder_chain]
    tgt_cat = [p for c, n, rn, an, p in at if c != binder_chain and (rn, an) in CATIONS]
    if not tgt_all:
        sys.exit(f"REFUSE: no target atoms outside chain {binder_chain}")
    iface = sorted(n for n, p in binder
                   if any(math.dist(p, q) <= IFACE_CA for q in tgt_all))
    reach = sorted(n for n, p in binder
                   if any(math.dist(p, q) <= REACH_CA for q in tgt_cat))
    return iface, reach


def backbone_rmsd(a_pdb, b_pdb, chain_a, chain_b=None):
    """CA RMSD between two structures WITHOUT superposition -- they share a frame only if the
    monomer was not re-centred. Reported so the caller can see whether index transfer is sound;
    a large value means the fold moved and the position sets may not apply."""
    ca = lambda p, c: {n: q for ch, n, rn, an, q in read_atoms(p) if ch == c and an == "CA"}
    A, B = ca(a_pdb, chain_a), ca(b_pdb, chain_b or chain_a)
    common = sorted(set(A) & set(B))
    if not common:
        return None, 0
    d2 = [math.dist(A[n], B[n]) ** 2 for n in common]
    return math.sqrt(sum(d2) / len(d2)), len(common)


def emit(name, binder_chain, iface, reach, out_dir, force_his=True):
    """Write ProteinMPNN omit_AA_jsonl. Format: {name: {chain: [[[1-based pos], 'AAs'], ...]}}"""
    os.makedirs(out_dir, exist_ok=True)
    items = []
    acid_only = [n for n in iface if n not in set(reach)]
    if acid_only:
        items.append([acid_only, ACIDS])
    if force_his and reach:
        # allow ONLY histidine: omit the other nineteen
        items.append([list(reach), "".join(c for c in ALL_AA if c != "H")])
    elif reach:
        items.append([list(reach), ACIDS])
    path = os.path.join(out_dir, "omit_AA.jsonl")
    with open(path, "a") as fh:
        fh.write(json.dumps({name: {binder_chain: items}}) + "\n")
    return path, items


def selftest():
    # the derived reach must equal its three measured parts, or the docstring is lying
    assert abs(REACH_CA - (MAX_CA_CB + MAX_CB_RING_N + D_HB)) < 0.01, REACH_CA
    assert IFACE_CA == 8.00, "IFACE_CA is his_cation_gate's SLACK_D"
    # synthetic complex: binder CAs marching away from a single target ARG NH1 at the origin
    import tempfile
    def line(ch, n, rn, an, x):
        return (f"ATOM  {n:5d}  {an:<3s} {rn} {ch}{n:4d}    "
                f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00\n")
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        for i, x in enumerate([2.0, 7.9, 8.6, 20.0], start=1):
            fh.write(line("A", i, "GLY", "CA", x))
        fh.write(line("B", 500, "ARG", "NH1", 0.0))
        fh.write(line("B", 501, "GLY", "CA", 0.0))
        tmp = fh.name
    iface, reach = positions(tmp, "A")
    # 2.0 and 7.9 are inside both cutoffs; 8.6 is outside IFACE 8.0 but inside REACH 8.74;
    # 20.0 is outside both. That gap is the point: a position can be reachable by a histidine
    # sidechain while its CA is not itself "at the interface".
    assert iface == [1, 2], iface
    assert reach == [1, 2, 3], reach
    # MUTATION: a target with NO cation must yield an empty reach set, not fall back to iface
    with open(tmp) as f:
        body = f.read().replace("ARG NH1", "ALA CB ").replace(" NH1", " CB ")
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        fh.write(body.replace("ARG", "ALA")); tmp2 = fh.name
    i2, r2 = positions(tmp2, "A")
    assert r2 == [], r2
    assert i2 == [1, 2], i2
    # emission: acid omission only where there is no histidine, and His forced where there is
    with tempfile.TemporaryDirectory() as d:
        _, items = emit("x", "A", iface, reach, d)
        by = {it[1]: it[0] for it in items}
        # every interface position here is also a His site, so there is NO acid item at all.
        # Writing an empty [] would hand ProteinMPNN a no-op constraint and make the file's
        # item count a lie about how many positions were constrained.
        assert ACIDS not in by, f"an empty acid constraint must not be written: {by}"
        his = [k for k in by if "H" not in k]
        assert len(his) == 1 and by[his[0]] == [1, 2, 3], by
        assert len(his[0]) == 19, f"forcing His must omit nineteen, got {len(his[0])}"
        # MUTATION: a non-reachable interface position must get the acid omission
        _, items2 = emit("y", "A", [1, 2, 3, 9], [1], d)
        by2 = {it[1]: it[0] for it in items2}
        assert by2[ACIDS] == [2, 3, 9], by2
    os.unlink(tmp); os.unlink(tmp2)
    print(f"  ok  REACH_CA {REACH_CA} = CA-CB {MAX_CA_CB} + CB-ringN {MAX_CB_RING_N} + D_HB {D_HB}")
    print(f"  ok  IFACE_CA {IFACE_CA} is his_cation_gate's own SLACK_D")
    print(f"  ok  a CA can be reachable ({REACH_CA}) without being at the interface ({IFACE_CA})")
    print(f"  ok  MUTATION: no cation on the target -> empty reach set, no fallback")
    print(f"  ok  MUTATION: forcing His omits exactly nineteen residues")
    print(f"  ok  MUTATION: acid omission covers interface positions that are NOT His sites")
    print("\nself-tests passed: 6")


if __name__ == "__main__":
    a = sys.argv
    if len(a) == 1 or "--selftest" in a:
        selftest()
    else:
        cx = a[a.index("--complex") + 1]
        bc = a[a.index("--binder") + 1] if "--binder" in a else "A"
        od = a[a.index("--out-dir") + 1]
        name = os.path.basename(cx)[:-4]
        iface, reach = positions(cx, bc)
        path, items = emit(name, bc, iface, reach, od)
        print(f"{name}: {len(iface)} interface positions, {len(reach)} cation-reachable")
        print(f"  acid-omitted : {[n for n in iface if n not in set(reach)]}")
        print(f"  His-forced   : {reach}")
        print(f"  -> {path}")
        if "--monomer" in a:
            r, n = backbone_rmsd(cx, a[a.index("--monomer") + 1], bc)
            print(f"  backbone CA RMSD vs the folded monomer: "
                  f"{'n/a' if r is None else f'{r:.2f} A over {n} residues'}")
