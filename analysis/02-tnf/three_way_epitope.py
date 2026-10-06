#!/usr/bin/env python3
"""Objective 2, finally: can one binder satisfy human AND mouse -- and now also avoid LT-alpha?

WHY THIS WAS OVERDUE. Mouse cross-reactivity is the SECOND-RANKED objective and it has had one
day-1 sequence analysis (s2) all week; everything since went at the pH switch. And s14 added a
third protein to the problem after the fact: lymphotoxin-alpha binds the same groove, so the
design now has to hold TWO proteins and release a THIRD, and nothing has looked at all three
together.

WHAT THIS COMPUTES, at each of s1's 21 consensus epitope positions:
  - human / mouse / LT-alpha residue, by alignment (s2's method, extended to a third sequence)
  - sidechain VOLUME change, because a pocket shaped to human must still admit mouse (s2's
    H149Y argument, now applied to every position instead of the two that were noticed)
  - CHARGE change, which is what the s4 mechanism actually depends on
  - the structural check s2 never did: superpose the real mouse trimer (2TNF, verified clean in
    s12) on canonical human and measure how far each epitope sidechain actually moves

THE THREE-WAY LOGIC, stated before the numbers:
  conserved human->mouse   = SAFE to design against (objective 2 survives)
  divergent human->mouse   = HAZARD; a contact here costs cross-reactivity
  divergent human->LT-alpha = SELECTIVITY HANDLE (s14's axis)
  conserved in BOTH mouse and LT-alpha = safe for objective 2, useless for selectivity
The interesting set is the intersection: conserved in mouse AND divergent in LT-alpha. Those are
the only positions that buy cross-reactivity and selectivity at the same time.
"""
import os, sys, urllib.request
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from species_and_histidines import HUMAN, MOUSE, UNIPROT_START, align
import fetch

LTA_CACHE = os.path.join(HERE, "structures", "P01374.fasta")
CONSENSUS = [96, 97, 107, 108, 109, 149, 151, 153, 161, 162, 163, 166, 167, 173,
             189, 191, 219, 220, 221, 222, 225]
ANCHORS = {108: "R108 secondary", 166: "K166 PRIMARY"}
# sidechain volumes, A^3 (Zamyatnin 1972, standard table)
VOL = {'G':60,'A':89,'S':89,'C':109,'D':111,'P':113,'N':114,'T':116,'E':138,'V':140,
       'Q':144,'H':153,'M':163,'I':167,'L':167,'K':169,'R':174,'F':190,'Y':194,'W':228}
CHARGE = {'D':-1,'E':-1,'K':+1,'R':+1,'H':0}        # His neutral at 7.4


def lta():
    if not os.path.exists(LTA_CACHE):
        urllib.request.urlretrieve(f"https://rest.uniprot.org/uniprotkb/P01374.fasta", LTA_CACHE)
    return "".join(l.strip() for l in open(LTA_CACHE) if not l.startswith(">"))


def mapping(ref, other):
    A, B = align(ref, other)
    pos, m = UNIPROT_START - 1, {}
    for a, b in zip(A, B):
        if a != "-":
            pos += 1
            m[pos] = b
    return m


def epitope_rmsd():
    """Superpose the real mouse trimer on canonical human; per-epitope-residue CA displacement.

    s2 compared SEQUENCES. This asks whether the divergent residues sit where a binder would
    feel them, which a sequence cannot say."""
    import gemmi
    try:
        hs = gemmi.read_structure(os.path.join(HERE, "..", "..", "targets", "tnf",
                                               "tnf_canonical_trimer.pdb"))
        ms = gemmi.read_structure(fetch.pdb("2TNF"))
    except Exception as e:
        return None, f"structures unavailable: {e}"
    for st in (hs, ms):
        st.setup_entities(); st.remove_ligands_and_waters()
    hch = hs[0]["A"]; mch = ms[0]["A"]
    # map by sequence position: human uniprot = seqid+76; mouse observed[i] == MOUSE[i+8]
    hca = {r.seqid.num + 76: r for r in hch if r.find_atom("CA", "*")}
    mobs = [r for r in mch if r.find_atom("CA", "*")]
    mca = {}                      # mouse uniprot position (P06804) -> residue
    for i, r in enumerate(mobs):
        mca[i + 8 + 80] = r
    # align human->mouse positions through the pairwise alignment
    h2m = {}
    A, B = align(HUMAN, MOUSE)
    hp, mp = UNIPROT_START - 1, 80 - 1
    for a, b in zip(A, B):
        if a != "-": hp += 1
        if b != "-": mp += 1
        if a != "-" and b != "-": h2m[hp] = mp
    pairs = [(hca[h], mca[h2m[h]]) for h in CONSENSUS
             if h in hca and h in h2m and h2m[h] in mca]
    if len(pairs) < 8:
        return None, f"only {len(pairs)} epitope pairs mapped"
    P = np.array([[r.find_atom("CA", "*").pos.x, r.find_atom("CA", "*").pos.y,
                   r.find_atom("CA", "*").pos.z] for r, _ in pairs])
    Q = np.array([[r.find_atom("CA", "*").pos.x, r.find_atom("CA", "*").pos.y,
                   r.find_atom("CA", "*").pos.z] for _, r in pairs])
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    Qr = Qc @ R.T
    per = np.linalg.norm(Pc - Qr, axis=1)
    return (float(np.sqrt((per ** 2).mean())), dict(zip([h for h in CONSENSUS
            if h in hca and h in h2m and h2m[h] in mca], per))), None


def main():
    L = lta()
    m_mouse, m_lta = mapping(HUMAN, MOUSE), mapping(HUMAN, L)
    rms, err = epitope_rmsd()

    print(f"{'pos':>5} {'hu':^4}{'mo':^4}{'LT':^4} {'dVol(mo)':>9}{'dChg(mo)':>9}"
          f"{'dVol(LT)':>9}{'dChg(LT)':>9}  {'CA shift':>9}  verdict")
    print("-" * 104)
    safe_sel, hazard, safe_only = [], [], []
    for u in CONSENSUS:
        h = HUMAN[u - UNIPROT_START]; mo = m_mouse.get(u, "-"); lt = m_lta.get(u, "-")
        dv_m = VOL.get(mo, 0) - VOL.get(h, 0) if mo != "-" else 0
        dc_m = CHARGE.get(mo, 0) - CHARGE.get(h, 0) if mo != "-" else 0
        dv_l = VOL.get(lt, 0) - VOL.get(h, 0) if lt != "-" else 0
        dc_l = CHARGE.get(lt, 0) - CHARGE.get(h, 0) if lt != "-" else 0
        shift = rms[1].get(u) if rms else None
        mo_same, lt_same = (mo == h), (lt == h)
        if mo_same and not lt_same:
            v = "SAFE + SELECTIVE"; safe_sel.append(u)
        elif mo_same:
            v = "safe, no selectivity"; safe_only.append(u)
        else:
            v = "HAZARD (mouse differs)"; hazard.append(u)
        tag = "  <<< " + ANCHORS[u] if u in ANCHORS else ""
        print(f"{u:>5} {h:^4}{mo:^4}{lt:^4} {dv_m:>+9}{dc_m:>+9}{dv_l:>+9}{dc_l:>+9}"
              f"  {shift if shift is None else round(shift,2):>9}  {v}{tag}")

    print("-" * 104)
    if rms:
        print(f"epitope CA RMSD, real mouse trimer superposed on canonical human: {rms[0]:.2f} A "
              f"over {len(rms[1])} positions")
    else:
        print(f"structural superposition unavailable: {err}")
    print(f"\n  SAFE + SELECTIVE (conserved in mouse, divergent in LT-alpha): {safe_sel}")
    print(f"  safe but not selective (conserved in both):                   {safe_only}")
    print(f"  HAZARD (mouse differs):                                       {hazard}")
    print(f"\n  {len(safe_sel)} of {len(CONSENSUS)} consensus positions satisfy BOTH objective 2 "
          f"and selectivity.")
    for u in ANCHORS:
        where = ("SAFE + SELECTIVE" if u in safe_sel else
                 "safe, no selectivity" if u in safe_only else "HAZARD")
        print(f"    {ANCHORS[u]:<16} -> {where}")


if __name__ == "__main__":
    main()
