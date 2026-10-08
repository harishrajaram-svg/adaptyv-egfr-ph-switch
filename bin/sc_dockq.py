#!/usr/bin/env python3
"""CA-based DockQ between a DESIGNED complex and its PREDICTED co-fold.

WHY THIS EXISTS, AND IT SHOULD HAVE EXISTED FIRST. Every margin in this work asks whether the
predicted interface scores better than a shuffled control. None of them asks whether the binder
ended up WHERE IT WAS DESIGNED TO GO. A design can score well with the predictor parking the
binder on some other surface it happens to like, and our instrument cannot tell that apart from
success at the intended epitope.

The organisers' own protocol ranks on `ipSAE` AND `sc_DockQ` together, 4:1, at every promotion
step, and defines sc_DockQ as DockQ(designed_complex, prediction). We computed the first and
never the second.

DockQ, Basu & Wallner 2016:
    DockQ = ( Fnat + 1/(1+(LRMS/8.5)^2) + 1/(1+(iRMS/1.5)^2) ) / 3
with the receptor superposed, LRMS the ligand RMSD after that superposition, iRMS the interface
RMSD, and Fnat the fraction of designed interface contacts the prediction preserves.

CA-ONLY, AND SAID OUT LOUD. The Genie 3 backbone carries one atom per residue, so every term here
is computed on CA positions and the contact definition is CA-CA within CONTACT_CA. That is a
coarser DockQ than the all-atom original and its absolute values are not comparable to published
all-atom numbers. It is used only to compare designs with each other on one consistent basis.

SYMMETRY. The target is a homotrimer, so the designed and predicted target chains can correspond
under any permutation. The protocol says to take "the best symmetric relabeling (the one that
maximizes DockQ)"; this tries all permutations of the target chains and reports the maximum.

    sc_dockq.py --selftest
    sc_dockq.py --designed <genie.pdb> --predicted <cofold.cif> [--design-binder A] [--pred-binder D]
"""
import itertools
import math
import sys

CONTACT_CA = 8.0        # A, CA-CA contact definition for Fnat
IFACE_CA = 10.0         # A, residues counted as interface for iRMS

# 🔴 THE CONTROL THAT MUST PASS BEFORE ANY DockQ IS REPORTED. Both structures hold the
# same target, so superposing the target chains on each other must give a SMALL RMSD -- TNF-alpha
# is a solved structure and ESMFold2 reproduces it to a couple of Angstroms. On 2026-10-08 this
# tool reported sc_DockQ ~0.03 with Fnat 0.00 and ligand RMSD of 22-50 A across all 48 designs,
# which read as "the predictor never puts the binder where it was designed". The target-only
# control then came back at 13.9-19.3 A, which is impossible for a solved fold, and the cause was
# a RESIDUE REGISTER SHIFT: the designed target is 152 residues numbered 1-152 (mature 6-157)
# while the predicted target is 157 numbered 1-157 (mature 1-157), so residue 1 was matched to
# residue 1 and every pair was five residues out. The numbers were a bug, not a finding.
# This tool now REFUSES rather than reporting when its own control fails.
TARGET_RMSD_MAX = 4.0   # A, target-only CA RMSD above which the correspondence is wrong


def read_ca(path, want=None):
    """{(chain, resnum): (x,y,z)} for CA atoms, from PDB or mmCIF."""
    out = {}
    if path.endswith(".cif"):
        import gemmi
        st = gemmi.read_structure(path); st.setup_entities()
        for ch in st[0]:
            for r in ch:
                a = r.find_atom("CA", "*")
                if a:
                    out[(ch.name, r.seqid.num)] = (a.pos.x, a.pos.y, a.pos.z)
    else:
        for l in open(path):
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                out[(l[21], int(l[22:26]))] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return out


def kabsch_rt(P, Q):
    """Rotation+translation taking P onto Q (both lists of 3-tuples, equal length)."""
    n = len(P)
    cp = [sum(p[i] for p in P)/n for i in range(3)]
    cq = [sum(q[i] for q in Q)/n for i in range(3)]
    A = [[sum((P[k][i]-cp[i])*(Q[k][j]-cq[j]) for k in range(n)) for j in range(3)] for i in range(3)]
    R = _svd_rot(A)
    return R, cp, cq


def _svd_rot(A):
    """Rotation from the 3x3 covariance A by Jacobi eigen-decomposition of A^T A. Pure stdlib."""
    AtA = [[sum(A[k][i]*A[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    V = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    M = [row[:] for row in AtA]
    for _ in range(64):
        p, q = max(((i, j) for i in range(3) for j in range(3) if i < j),
                   key=lambda ij: abs(M[ij[0]][ij[1]]))
        if abs(M[p][q]) < 1e-12:
            break
        th = 0.5 * math.atan2(2*M[p][q], M[p][p]-M[q][q])
        c, s = math.cos(th), math.sin(th)
        for k in range(3):
            mkp, mkq = M[k][p], M[k][q]
            M[k][p], M[k][q] = c*mkp + s*mkq, -s*mkp + c*mkq
        for k in range(3):
            mpk, mqk = M[p][k], M[q][k]
            M[p][k], M[q][k] = c*mpk + s*mqk, -s*mpk + c*mqk
        for k in range(3):
            vkp, vkq = V[k][p], V[k][q]
            V[k][p], V[k][q] = c*vkp + s*vkq, -s*vkp + c*vkq
    ev = [(M[i][i], i) for i in range(3)]
    ev.sort(reverse=True)
    v = [[V[r][i] for r, _ in [(0, 0)]] for i in range(3)]   # placeholder, rebuilt below
    cols = [[V[r][idx] for r in range(3)] for _, idx in ev]
    # U = A V S^-1
    U = []
    for k, (lam, _) in enumerate(ev):
        s = math.sqrt(max(lam, 1e-18))
        U.append([sum(A[i][j]*cols[k][j] for j in range(3))/s for i in range(3)])
    # R = U V^T with a determinant fix
    R = [[sum(U[k][i]*cols[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    if _det(R) < 0:
        U[2] = [-x for x in U[2]]
        R = [[sum(U[k][i]*cols[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    return R


def _det(R):
    return (R[0][0]*(R[1][1]*R[2][2]-R[1][2]*R[2][1])
            - R[0][1]*(R[1][0]*R[2][2]-R[1][2]*R[2][0])
            + R[0][2]*(R[1][0]*R[2][1]-R[1][1]*R[2][0]))


def apply_rt(p, R, cp, cq):
    d = [p[i]-cp[i] for i in range(3)]
    return tuple(sum(R[j][i]*d[j] for j in range(3)) + cq[i] for i in range(3))


def dockq_once(dca, pca, d_bind, p_bind, chain_map, offset=0):
    """DockQ for one target-chain correspondence. chain_map: predicted chain -> designed chain.

    `offset` is added to the DESIGNED residue number to reach the predicted numbering, for the
    case where the two constructs start at different mature positions.
    """
    dt = [k for k in dca if k[0] != d_bind]
    pt = [k for k in pca if k[0] != p_bind]
    pairs = [(k, (chain_map[k[0]], k[1] - offset)) for k in pt]
    pairs = [(p, d) for p, d in pairs if d in dca]
    if len(pairs) < 10:
        return None
    R, cp, cq = kabsch_rt([pca[p] for p, _ in pairs], [dca[d] for _, d in pairs])
    # the control: the shared target must superpose tightly or the correspondence is wrong
    tgt_rmsd = math.sqrt(sum(sum((apply_rt(pca[p], R, cp, cq)[i] - dca[d][i])**2
                                 for i in range(3)) for p, d in pairs) / len(pairs))
    if tgt_rmsd > TARGET_RMSD_MAX:
        return None
    db = sorted(k for k in dca if k[0] == d_bind)
    pb = sorted(k for k in pca if k[0] == p_bind)
    n = min(len(db), len(pb))
    if n < 5:
        return None
    moved = {k: apply_rt(pca[k], R, cp, cq) for k in pca}
    lrms = math.sqrt(sum(sum((moved[pb[i]][j]-dca[db[i]][j])**2 for j in range(3))
                         for i in range(n))/n)
    # designed contacts
    dcon = {(b, t) for b in db for t in dt if math.dist(dca[b], dca[t]) <= CONTACT_CA}
    if not dcon:
        return None
    inv = {v: k for k, v in chain_map.items()}
    pcon = set()
    for b in pb:
        for t in pt:
            if math.dist(pca[b], pca[t]) <= CONTACT_CA:
                pcon.add((db[pb.index(b)], (chain_map[t[0]], t[1] - offset)))
    fnat = len(dcon & pcon)/len(dcon)
    # interface residues in the DESIGNED complex, RMSD after the same superposition
    ifr = sorted({t for _, t in dcon} | {b for b, _ in dcon},
                 key=lambda k: (k[0], k[1]))
    tot, cnt = 0.0, 0
    for k in ifr:
        pk = ((inv.get(k[0], k[0]), k[1] + offset) if k[0] != d_bind
              else (pb[db.index(k)] if k in db else None))
        if pk is None or pk not in moved:
            continue
        tot += sum((moved[pk][j]-dca[k][j])**2 for j in range(3)); cnt += 1
    if cnt < 3:
        return None
    irms = math.sqrt(tot/cnt)
    return (fnat + 1/(1+(lrms/8.5)**2) + 1/(1+(irms/1.5)**2))/3, fnat, lrms, irms


def sc_dockq(designed, predicted, d_bind="A", p_bind="D", offsets=(0, 5, -5)):
    """Max DockQ over all target-chain permutations AND candidate numbering offsets.

    Returns None if no combination passes the target-superposition control, which means the two
    constructs could not be put in correspondence -- never a DockQ of zero.
    """
    dca, pca = read_ca(designed), read_ca(predicted)
    dt = sorted({k[0] for k in dca if k[0] != d_bind})
    pt = sorted({k[0] for k in pca if k[0] != p_bind})
    if len(dt) != len(pt):
        return None
    best = None
    for off in offsets:
        for perm in itertools.permutations(dt):
            r = dockq_once(dca, pca, d_bind, p_bind, dict(zip(pt, perm)), offset=off)
            if r and (best is None or r[0] > best[0]):
                best = r
    return best


def selftest():
    import tempfile, os, random
    def write(path, binder, target, bchain="A", tchains="BCD"):
        with open(path, "w") as fh:
            i = 1
            for n, p in enumerate(binder, start=1):
                fh.write(f"ATOM  {i:5d}  CA  GLY {bchain}{n:4d}    "
                         f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00\n"); i += 1
            for c, chain in zip(tchains, target):
                for n, p in enumerate(chain, start=1):
                    fh.write(f"ATOM  {i:5d}  CA  GLY {c}{n:4d}    "
                             f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00\n"); i += 1
    rng = random.Random(7)
    binder = [(rng.uniform(0, 12), rng.uniform(0, 12), rng.uniform(0, 12)) for _ in range(20)]
    target = [[(rng.uniform(10, 24), rng.uniform(0, 12), rng.uniform(0, 12)) for _ in range(25)]
              for _ in range(3)]
    d = tempfile.mkdtemp()
    a, b = os.path.join(d, "a.pdb"), os.path.join(d, "b.pdb")
    write(a, binder, target); write(b, binder, target, bchain="D", tchains="ABC")
    r = sc_dockq(a, b, "A", "D")
    assert r and r[0] > 0.99, f"a structure against its own copy must give DockQ 1.0, got {r}"
    print(f"  ok  identical complex -> DockQ {r[0]:.4f}, Fnat {r[1]:.3f}, "
          f"LRMS {r[2]:.3f}, iRMS {r[3]:.3f}")
    # MUTATION: move the binder far away; DockQ must collapse
    far = [(p[0]+40, p[1], p[2]) for p in binder]
    c = os.path.join(d, "c.pdb"); write(c, far, target, bchain="D", tchains="ABC")
    r2 = sc_dockq(a, c, "A", "D")
    assert r2 and r2[0] < 0.1, f"a displaced binder must collapse DockQ, got {r2}"
    print(f"  ok  MUTATION: binder displaced 40 A -> DockQ {r2[0]:.4f}, Fnat {r2[1]:.3f}")
    # MUTATION: a rigid-body move of the WHOLE complex must not change DockQ
    sh = [(p[0]+100, p[1]-50, p[2]+7) for p in binder]
    st = [[(p[0]+100, p[1]-50, p[2]+7) for p in ch] for ch in target]
    e = os.path.join(d, "e.pdb"); write(e, sh, st, bchain="D", tchains="ABC")
    r3 = sc_dockq(a, e, "A", "D")
    assert r3 and abs(r3[0]-r[0]) < 1e-3, f"DockQ must be invariant to rigid motion: {r3[0]} vs {r[0]}"
    print(f"  ok  MUTATION: whole complex moved 100 A -> DockQ unchanged at {r3[0]:.4f}")
    # MUTATION: a register-shifted copy must be REFUSED by the control, not scored as bad.
    # This is the 2026-10-08 bug: designed 1-152 against predicted 1-157, five residues out.
    f = os.path.join(d, "f.pdb")
    with open(f, "w") as fh:
        i = 1
        for n, p in enumerate(binder, start=1):
            fh.write(f"ATOM  {i:5d}  CA  GLY D{n:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00\n"); i += 1
        for c, chain in zip("ABC", target):
            for n, p in enumerate(chain, start=6):      # starts at 6, not 1
                fh.write(f"ATOM  {i:5d}  CA  GLY {c}{n:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00\n"); i += 1
    r4 = sc_dockq(a, f, "A", "D", offsets=(0,))          # offset 0 only: must REFUSE
    assert r4 is None, f"a 5-residue register shift must be refused, got {r4}"
    r5 = sc_dockq(a, f, "A", "D", offsets=(0, 5, -5))    # with the offset tried: must recover
    assert r5 and r5[0] > 0.99, f"with the right offset it must recover 1.0, got {r5}"
    print(f"  ok  MUTATION: a 5-residue register shift is REFUSED (not scored 0), and recovers")
    print(f"      to DockQ {r5[0]:.4f} once the offset is tried")
    print(f"  ok  the target-superposition control refuses above {TARGET_RMSD_MAX} A")
    print(f"  ok  CA-only, contact cut {CONTACT_CA} A -- not comparable to all-atom DockQ")
    print("\nself-tests passed: 6")
    return 0


if __name__ == "__main__":
    a = sys.argv
    if len(a) == 1 or "--selftest" in a:
        sys.exit(selftest())
    r = sc_dockq(a[a.index("--designed")+1], a[a.index("--predicted")+1],
                 a[a.index("--design-binder")+1] if "--design-binder" in a else "A",
                 a[a.index("--pred-binder")+1] if "--pred-binder" in a else "D")
    print(f"{r[0]:.4f}\t{r[1]:.4f}\t{r[2]:.4f}\t{r[3]:.4f}" if r else "REFUSED")
