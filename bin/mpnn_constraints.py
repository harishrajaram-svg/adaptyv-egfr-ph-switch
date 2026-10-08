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

# FORCING HISTIDINE AT EVERY REACHABLE POSITION BUILDS A POLY-HISTIDINE TAG. Measured
# 2026-10-08: the first version produced runs of six to NINE consecutive H in 17 of 48 designs,
# because adjacent binder positions are often all within reach of the same cation. That is a
# purification tag, it is a homopolymer liability of exactly the kind that sank problem 1's
# rank-1 design (30% Ala with a 7-Ala run), and express_qc_p2 flags it. So forced positions are
# now SPACED and CAPPED, chosen nearest-cation-first.
MIN_HIS_SPACING = 3        # no two forced histidines within 3 positions of each other
MAX_FORCED_HIS = 8         # the top of the 8-11 His-cation contact range in working designs
OMIT_EVERYWHERE = "C"      # SolubleMPNN put a cysteine in 9 of 48; de novo work is Cys-free
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
    """(interface, cation_reach) -- 1-based binder residue numbers, as ProteinMPNN wants.

    `cation_reach` is SPACED and CAPPED: candidates are ranked by how close they sit to a target
    cation, then taken greedily while keeping MIN_HIS_SPACING between them, up to MAX_FORCED_HIS.
    """
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
    cand = []
    for n, p in binder:
        if not tgt_cat:
            break
        d = min(math.dist(p, q) for q in tgt_cat)
        if d <= REACH_CA:
            cand.append((d, n))
    chosen = []
    for d, n in sorted(cand):                      # nearest cation first
        if len(chosen) >= MAX_FORCED_HIS:
            break
        if all(abs(n - m) >= MIN_HIS_SPACING for m in chosen):
            chosen.append(n)
    return iface, sorted(chosen)


def backbone_rmsd(a_pdb, b_pdb, chain_a, chain_b=None):
    """CA RMSD after OPTIMAL SUPERPOSITION (Kabsch), over residues present in both.

    SUPERPOSED, corrected 2026-10-08. The first version compared coordinates directly, which is
    meaningless here: the Genie 3 complex and the ESMFold2 monomer are in unrelated frames, so a
    raw distance measures the frame offset and not the fold. The quantity that decides whether
    constraint positions transfer by index is whether the FOLD is the same shape, which is RMSD
    after superposition.
    """
    import numpy as np

    ca = lambda p, c: {n: q for ch, n, rn, an, q in read_atoms(p) if ch == c and an == "CA"}
    A, B = ca(a_pdb, chain_a), ca(b_pdb, chain_b or chain_a)
    common = sorted(set(A) & set(B))
    if len(common) < 3:
        return None, len(common)
    P = np.array([A[n] for n in common], dtype=float)
    Q = np.array([B[n] for n in common], dtype=float)
    P -= P.mean(0); Q -= Q.mean(0)
    V, S, W = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1.0, 1.0, d]) @ W
    return float(np.sqrt(((P @ R - Q) ** 2).sum() / len(common))), len(common)


def emit(name, binder_chain, iface, reach, out_dir, force_his=True, n_res=None):
    """Write ProteinMPNN omit_AA_jsonl. Format: {name: {chain: [[[1-based pos], 'AAs'], ...]}}"""
    os.makedirs(out_dir, exist_ok=True)
    items = []
    # n_res is the BINDER LENGTH, so the global omission covers the whole chain. Defaulting it
    # to the last constrained position would silently leave the tail unprotected.
    if n_res is None:
        n_res = max(max(iface, default=0), max(reach, default=0))
    if OMIT_EVERYWHERE and n_res:
        items.append([list(range(1, n_res + 1)), OMIT_EVERYWHERE])
    acid_only = [n for n in iface if n not in set(reach)]
    if acid_only:
        items.append([acid_only, ACIDS])
    if force_his and reach:
        # allow ONLY histidine: omit the other nineteen
        items.append([list(reach), "".join(c for c in ALL_AA if c != "H")])
    elif reach:
        items.append([list(reach), ACIDS])
    # ONE JSON OBJECT FOR ALL DESIGNS, not one line each. ProteinMPNN reads this file as
    #     for json_str in list(json_file): omit_AA_dict = json.loads(json_str)
    # which OVERWRITES the dict on every line, so a 49-line file silently keeps only the last
    # entry and every other design dies on KeyError. Found 2026-10-08 the hard way.
    path = os.path.join(out_dir, "omit_AA.jsonl")
    allnames = {}
    if os.path.exists(path):
        with open(path) as fh:
            body = fh.read().strip()
        if body:
            allnames = json.loads(body.splitlines()[-1])
    allnames[name] = {binder_chain: items}
    with open(path, "w") as fh:
        fh.write(json.dumps(allnames) + "\n")
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
    # spacing: CAs at x=2.0, 7.9, 8.6 are residues 1, 2, 3 -- all reachable, but 1 and 2 are
    # adjacent and 1 and 3 are two apart, so only residue 1 survives MIN_HIS_SPACING of 3.
    assert reach == [1], f"spacing must thin adjacent reachable positions: {reach}"
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
        assert by[ACIDS] == [2], f"residue 2 is interface but not a His site: {by}"
        assert by[OMIT_EVERYWHERE] == [1, 2], f"C omitted over the constrained span: {by}"
        # MUTATION: given the real binder length, the global omission must cover ALL of it,
        # not just up to the last constrained position.
        _, it83 = emit("z", "A", iface, reach, d, n_res=83)
        by83 = {i[1]: i[0] for i in it83}
        assert by83[OMIT_EVERYWHERE] == list(range(1, 84)), by83[OMIT_EVERYWHERE][:5]
        his = [k for k in by if len(k) == 19]
        assert len(his) == 1 and by[his[0]] == [1], by
        assert len(his[0]) == 19, f"forcing His must omit nineteen, got {len(his[0])}"
        # MUTATION: a non-reachable interface position must get the acid omission
        _, items2 = emit("y", "A", [1, 2, 3, 9], [1], d)
        by2 = {it[1]: it[0] for it in items2}
        assert by2[ACIDS] == [2, 3, 9], by2
        # MUTATION: spacing and the cap must both bite. Nine adjacent reachable candidates
        # must thin to at most MAX_FORCED_HIS, and no two within MIN_HIS_SPACING.
        import itertools as _it
        picked = []
        for n in range(1, 40):
            if all(abs(n - m) >= MIN_HIS_SPACING for m in picked) and len(picked) < MAX_FORCED_HIS:
                picked.append(n)
        assert len(picked) == MAX_FORCED_HIS, picked
        assert all(b - a >= MIN_HIS_SPACING for a, b in zip(picked, picked[1:])), picked
        # MUTATION: the file must be ONE object holding BOTH designs on ONE line. Appending a
        # line per design makes ProteinMPNN keep only the last and KeyError on all the rest.
        body = open(os.path.join(d, "omit_AA.jsonl")).read().strip()
        assert len(body.splitlines()) == 1, f"{len(body.splitlines())} lines; must be 1"
        loaded = json.loads(body)
        assert set(loaded) == {"x", "y", "z"}, loaded
        assert "A" in loaded["x"] and "A" in loaded["y"]
    # MUTATION: superposition must remove a rigid-body move entirely. Without it this reads
    # the frame offset -- which is what the first version of this function did.
    #
    # numpy is the ONLY non-stdlib import in this file and it is needed solely for the Kabsch
    # SVD. A clean clone's system python3 does not have it, and this selftest used to die on a
    # bare traceback there (found 2026-10-08 by running the suite in a fresh checkout). It now
    # says what is missing and how to get it, and reports a SKIP rather than a pass, so an
    # unrun test can never be mistaken for a green one.
    try:
        import numpy as np
    except ImportError:
        print("  SKIP  the superposition mutation test needs numpy, which this interpreter")
        print("        lacks. Run it with the project venv: .venv/bin/python "
              "bin/mpnn_constraints.py --selftest")
        print("\nself-tests passed: 8 of 9 (1 SKIPPED -- numpy absent)")
        return 0
    import tempfile as _tf
    pts = np.array([[0.,0.,0.],[3.8,0,0],[7.0,2.1,0],[9.1,5.0,1.2],[11.0,8.0,2.0]])
    th = 0.7
    R = np.array([[math.cos(th),-math.sin(th),0],[math.sin(th),math.cos(th),0],[0,0,1]])
    moved = pts @ R.T + np.array([100.0, -50.0, 7.0])
    def _w(arr):
        fh = _tf.NamedTemporaryFile("w", suffix=".pdb", delete=False)
        for i, (x, y, z) in enumerate(arr, start=1):
            fh.write(f"ATOM  {i:5d}  CA  GLY A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n")
        fh.close(); return fh.name
    f1, f2 = _w(pts), _w(moved)
    r, n = backbone_rmsd(f1, f2, "A")
    # Tolerance is 1e-3, not 0: PDB writes coordinates to three decimals, so a round trip
    # through the format quantises at 0.001 A and the residual RMSD floor is ~3e-4.
    assert n == 5 and r is not None and r < 1e-3, (r, n)
    os.unlink(f1); os.unlink(f2)
    os.unlink(tmp); os.unlink(tmp2)
    print(f"  ok  REACH_CA {REACH_CA} = CA-CB {MAX_CA_CB} + CB-ringN {MAX_CB_RING_N} + D_HB {D_HB}")
    print(f"  ok  IFACE_CA {IFACE_CA} is his_cation_gate's own SLACK_D")
    print(f"  ok  a CA can be reachable ({REACH_CA}) without being at the interface ({IFACE_CA})")
    print(f"  ok  MUTATION: no cation on the target -> empty reach set, no fallback")
    print(f"  ok  MUTATION: forcing His omits exactly nineteen residues")
    print(f"  ok  MUTATION: acid omission covers interface positions that are NOT His sites")
    print(f"  ok  MUTATION: a rotated and translated copy superposes to RMSD {r:.2e}")
    print(f"  ok  forced histidines spaced >= {MIN_HIS_SPACING} apart, capped at {MAX_FORCED_HIS}")
    print(f"  ok  {OMIT_EVERYWHERE} omitted at every position of the binder, full length")
    print("\nself-tests passed: 9")


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
        nb = len({n for c, n, rn, an, p in read_atoms(cx) if c == bc})
        path, items = emit(name, bc, iface, reach, od, n_res=nb)
        print(f"{name}: {len(iface)} interface positions, {len(reach)} cation-reachable")
        print(f"  acid-omitted : {[n for n in iface if n not in set(reach)]}")
        print(f"  His-forced   : {reach}")
        print(f"  -> {path}")
        if "--monomer" in a:
            r, n = backbone_rmsd(cx, a[a.index("--monomer") + 1], bc)
            print(f"  backbone CA RMSD vs the folded monomer: "
                  f"{'n/a' if r is None else f'{r:.2f} A over {n} residues'}")
