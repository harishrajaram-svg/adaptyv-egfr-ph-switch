#!/usr/bin/env python3
"""rAC1 vs the 4UIP crystal interface, by STRUCTURAL CONTACT RECOVERY.

Reviewer, 2026-10-05: "Compare predictions with the 4UIP interface using structural
contact recovery. Do not attach predicted PAE to crystallographic coordinates."

The second half is the important half. A PAE matrix is a property of a PREDICTION -- it
describes the predictor's uncertainty about its own output. Scoring crystallographic
coordinates with a predicted PAE produces a number with no defined meaning, because the
crystal has no PAE. The only honest comparison between a prediction and a crystal is
geometric: does the prediction place the same residues in contact?

METHOD. An interface contact is a residue pair (target_i, binder_j) whose closest
heavy atoms lie within CUT angstroms. Compute that set for the crystal and for each
predicted pose, then report:

    recall     fraction of CRYSTAL contacts the prediction reproduces
    precision  fraction of PREDICTED contacts that are in the crystal
    epitope recall / binder-paratope recall, at residue rather than pair level

Residue numbering is NOT assumed to agree. Both sides are mapped through a global
sequence alignment of each chain to its crystal counterpart, so a cropped construct or a
renumbered model still lines up. Pairs involving a residue absent from the alignment are
dropped and counted, rather than silently treated as non-contacts.

    rac1_contact_recovery.py [--cut 5.0]
"""
import glob, json, os, sys
import gemmi

XTAL = 'targets/egfr/controls/4UIP.pdb'
POSE_GLOB = 'runs/esmfold2/w1_rac1/**/*.cif'
CUT = 5.0


def chains_with_seq(path):
    st = gemmi.read_structure(str(path))
    st.setup_entities(); st.remove_ligands_and_waters()
    out = []
    for ch in st[0]:
        res = [r for r in ch if r.name not in ('HOH',)]
        if len(res) < 10:
            continue
        out.append((ch.name, gemmi.one_letter_code([r.name for r in res]).upper(), res))
    return out


def align_map(seq_a, seq_b):
    """index_in_a -> index_in_b via global alignment. gemmi's aligner, no gap penalties tuned."""
    res = gemmi.align_string_sequences(list(seq_a), list(seq_b), [])
    # gemmi exposes the alignment as a CIGAR string, e.g. "12M3I40M". M consumes both
    # sequences, I consumes seq_a only, D consumes seq_b only.
    import re
    m, ia, ib = {}, 0, 0
    for length, op in re.findall(r'(\d+)([MID])', res.cigar_str()):
        length = int(length)
        if op == 'M':
            for _ in range(length):
                m[ia] = ib; ia += 1; ib += 1
        elif op == 'I':
            ia += length
        elif op == 'D':
            ib += length
    return m


def contacts(res_t, res_b, cut):
    """{(i_t, i_b)} residue-index pairs within cut A, heavy atoms only."""
    import math
    # bucket binder atoms on a grid for a cheap neighbour search
    cell = cut
    grid = {}
    for jb, rb in enumerate(res_b):
        for a in rb:
            if a.element == gemmi.Element('H'):
                continue
            k = (int(a.pos.x // cell), int(a.pos.y // cell), int(a.pos.z // cell))
            grid.setdefault(k, []).append((jb, a.pos))
    out = set()
    for it, rt in enumerate(res_t):
        for a in rt:
            if a.element == gemmi.Element('H'):
                continue
            kx, ky, kz = int(a.pos.x // cell), int(a.pos.y // cell), int(a.pos.z // cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for jb, p in grid.get((kx+dx, ky+dy, kz+dz), ()):
                            if a.pos.dist(p) <= cut:
                                out.add((it, jb))
    return out


def main():
    cut = CUT
    global POSE_GLOB
    if '--cut' in sys.argv:
        cut = float(sys.argv[sys.argv.index('--cut') + 1])
    # --poses lets the SAME geometric test run on another predictor's structures.
    # Chai-1 emits no residue-level PAE (its npz carries aggregate_score, ptm, iptm,
    # per_chain_pair_iptm and clashes only), so ipSAE cannot be computed on its output
    # and a like-for-like SCORE comparison is impossible. The geometric test does not
    # need a score, which is why it is the comparison that survives.
    if '--poses' in sys.argv:
        POSE_GLOB = sys.argv[sys.argv.index('--poses') + 1]
    if '--label' in sys.argv:
        print(f"[{sys.argv[sys.argv.index('--label') + 1]}]")
    xc = chains_with_seq(XTAL)
    if len(xc) < 2:
        sys.exit(f"{XTAL}: need 2 chains, got {[c[0] for c in xc]}")
    xc.sort(key=lambda c: -len(c[2]))
    (tn, tseq, tres), (bn, bseq, bres) = xc[0], xc[1]
    print(f"crystal {XTAL}: target chain {tn} ({len(tres)} res), "
          f"binder chain {bn} ({len(bres)} res), contact cutoff {cut} A")
    xtal = contacts(tres, bres, cut)
    ep_x = {i for i, _ in xtal}
    pa_x = {j for _, j in xtal}
    print(f"crystal interface: {len(xtal)} residue-residue contacts, "
          f"{len(ep_x)} epitope residues, {len(pa_x)} paratope residues\n")

    poses = sorted(glob.glob(POSE_GLOB, recursive=True))
    if not poses:
        sys.exit(f"no poses under {POSE_GLOB}")
    rows = []
    print(f"{'pose':<46}{'recall':>8}{'prec':>7}{'epi_rec':>9}{'par_rec':>9}{'n_pred':>8}")
    for p in poses:
        pc = chains_with_seq(p)
        if len(pc) < 2:
            print(f"{os.path.basename(p)[:45]:<46}  <2 chains, skipped"); continue
        pc.sort(key=lambda c: -len(c[2]))
        (ptn, ptseq, ptres), (pbn, pbseq, pbres) = pc[0], pc[1]
        mt = align_map(ptseq, tseq)      # predicted target idx -> crystal target idx
        mb = align_map(pbseq, bseq)
        pred_raw = contacts(ptres, pbres, cut)
        pred, dropped = set(), 0
        for i, j in pred_raw:
            if i in mt and j in mb:
                pred.add((mt[i], mb[j]))
            else:
                dropped += 1
        inter = pred & xtal
        rec = len(inter) / len(xtal) if xtal else 0.0
        pre = len(inter) / len(pred) if pred else 0.0
        epi = len({i for i, _ in pred} & ep_x) / len(ep_x) if ep_x else 0.0
        par = len({j for _, j in pred} & pa_x) / len(pa_x) if pa_x else 0.0
        rows.append(dict(pose=os.path.basename(p), recall=round(rec, 4),
                         precision=round(pre, 4), epitope_recall=round(epi, 4),
                         paratope_recall=round(par, 4), n_pred=len(pred),
                         n_dropped_unaligned=dropped))
        print(f"{os.path.basename(p)[:45]:<46}{rec:>8.3f}{pre:>7.3f}{epi:>9.3f}"
              f"{par:>9.3f}{len(pred):>8}")
    if rows:
        import statistics as st
        print(f"\nmedian over {len(rows)} poses: recall "
              f"{st.median([r['recall'] for r in rows]):.3f}, precision "
              f"{st.median([r['precision'] for r in rows]):.3f}, epitope recall "
              f"{st.median([r['epitope_recall'] for r in rows]):.3f}")
        tot_drop = sum(r['n_dropped_unaligned'] for r in rows)
        if tot_drop:
            print(f"  ({tot_drop} predicted contacts dropped across all poses because a "
                  f"residue had no crystal counterpart in the alignment)")
    out = ('analysis/01-egfr/rac1_contact_recovery.json' if '--poses' not in sys.argv
           else sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv
           else 'analysis/01-egfr/rac1_contact_recovery_alt.json')
    json.dump(dict(crystal=XTAL, cut=cut, pose_glob=POSE_GLOB, n_xtal_contacts=len(xtal),
                   n_epitope=len(ep_x), n_paratope=len(pa_x), poses=rows),
              open(out, 'w'), indent=1)
    print(f"\nwrote {out}")
    print("NOTE: no PAE or ipSAE value is attached to the crystal anywhere in this analysis.")


if __name__ == '__main__':
    main()
