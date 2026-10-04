#!/usr/bin/env python3
"""Structural + sequence novelty gate for the Adaptyv 2026 de novo filter.

WHY THIS IS A HARD GATE, NOT A RANKING
--------------------------------------
Measured against Adaptyv's own pre-computed FoldSeek results for the 393 tested
round-2 EGFR designs (results/foldseek_results.tsv in adaptyvbio/egfr_competition_2):

    median qtmscore vs PDB          0.883
    median qtmscore vs AFDB         0.891
    TM < 0.5 vs BOTH databases      3.6%   (14/393)
    sequence identity <= 30%        ~21%
    FULL Level-4 (both halves)      2.8%   (11/393)

And of the 53 designs that actually BOUND, exactly ONE would have passed.
Binders passed at 2%, non-binders at 4% -- novelty is mildly ANTI-correlated with
binding on this target. So this gate removes most of what a binder-optimised
pipeline produces, and it has to run BEFORE ranking compute, not after.

THE RULE (adaptyvbio.com/blog/novelty)
--------------------------------------
Level 4 "de novo" = sequence identity <= 30% AND less than moderate structural
similarity, where moderate = >70% of the sequence covered by domains matching a
known structure at TM >= 0.5.

APPROXIMATION, STATED PLAINLY: the official pipeline segments the design into
domains (consensus of 3 predictors, TED-style) and weights TM by coverage. This
script uses whole-chain qtmscore from a TM-align search. For a 55-80 aa single-domain
minibinder the two should agree closely, but they are not the same computation.

COVERAGE (fixed 2026-10-02): the ">70% covered" half of the rule was never applied.
alnlen was parsed into every hit and then never read, so a design failed on a high
TM over any-length segment. union_coverage() now implements it, and both verdicts
are reported -- passes_struct (official: TM AND coverage) and passes_struct_tmonly
(the stricter rule used until now). On the 20-design submission set the coverage
clause flips nothing, which is why the 0/120 novelty result stands either way.
The calibration numbers above were produced the same way, so our designs and the
round-2 baseline are at least measured on the same ruler.

USAGE
    python3 bin/novelty_gate.py design.pdb [design2.pdb ...] --db PATH/pdb
    python3 bin/novelty_gate.py designs/*.pdb --db ~/fsdb/pdb --tsv novelty.tsv
"""
import argparse, os, subprocess, sys, tempfile, shutil, csv

# ---------------------------------------------------------------------------------------
# ADAPTYV'S ACTUAL RULE, implemented verbatim.  https://www.adaptyvbio.com/blog/novelty
# (linked in the competition Slack by both Tudor Cotet and Simon as the rule they run.)
#
# Structural:   TM >= 0.80 = HIGH     TM >= 0.50 = MODERATE
# Sequence:     70% and 30% identity
#
# GENERAL PROTEINS
#   L1 Essentially Known  seq > 70% AND at least moderate structural
#   L2 Familiar           seq > 70% alone, OR **HIGH STRUCTURAL ALONE**, OR (seq > 30% AND moderate)
#   L3 Partly Novel       seq > 30% OR moderate structural -- exactly ONE of them
#   L4 De Novo            seq <= 30% AND less than moderate structural
#
# ANTIBODIES (nanobody, scFv, Fab, VHH, VNAR)
#   L1 CDRH3 >= 95%, OR CDRH3 >= 70% AND global >= 95%
#   L2 CDRH3 70-95% AND global < 95%
#   L3 CDRH3 < 70% AND global >= 70%
#   L4 CDRH3 < 70% AND global < 70%
#
# WHAT THIS SCRIPT USED TO DO, AND WHY IT WAS WRONG (fixed 2026-10-04):
# it tested a single TM < 0.50 bar and called everything above it a failure, reporting
# "0/20 pass". The submission gate is **level >= 3**, and a design with moderate (not high)
# structural similarity and <30% sequence identity is Level 3 and CLEARS. The true state of
# the 20-design submission was 19 clear / 1 does not -- `cf_short120_r031` at TM 0.811, which
# crosses the HIGH line by 0.011 and is Level 2 despite 19.7% sequence identity. That clause
# -- high structural similarity ALONE, with no sequence condition -- is the one we missed.
#
# TWO LIMITS OF THIS IMPLEMENTATION, stated rather than hidden:
#  1. Adaptyv run MMseqs2 against SwissProt, PDB, patent sequences, the therapeutic-antibody
#     database and PLAbDab. We search PDB only, via FoldSeek. Our sequence identities are
#     therefore LOWER BOUNDS; a design clean here may still hit a patent or SwissProt entry.
#  2. The antibody branch needs CDRH3 identity, which requires ANARCI numbering we do not run.
#     Antibody-format designs are reported as ANTIBODY/UNSCORED rather than given a level.
TM_HIGH = 0.80    # "high" structural similarity
TM_MOD  = 0.50    # "moderate" structural similarity
SEQ_HI  = 0.70
SEQ_MID = 0.30
GATE_LEVEL = 3    # submissions must clear level >= 3

TM_BAR = TM_MOD   # retained: callers and the old TSV columns reference it
FID_BAR = SEQ_MID
COV_BAR = 0.70    # ">70% of the sequence covered" -- used for the coverage-qualified variant
FOLDSEEK = os.environ.get("FOLDSEEK", "foldseek")

FMT = "query,target,fident,alnlen,qtmscore,ttmscore,alntmscore,evalue,prob,qstart,qend"


def chains_in(pdb):
    out = {}
    for line in open(pdb):
        if line[:4] == "ATOM":
            out.setdefault(line[21], set()).add(int(line[22:26]))
    return {c: len(v) for c, v in out.items()}


def extract_chain(pdb, chain, dest):
    with open(dest, "w") as fh:
        for line in open(pdb):
            if line[:4] == "ATOM" and line[21] == chain:
                fh.write(line)
        fh.write("END\n")
    return dest


def search(query_pdb, db, workdir):
    res = os.path.join(workdir, "res.tsv")
    tmp = os.path.join(workdir, "tmp")
    os.makedirs(tmp, exist_ok=True)
    cmd = [FOLDSEEK, "easy-search", query_pdb, db, res, tmp,
           "--alignment-type", "1",          # TMalign
           "--format-output", FMT,
           "--max-seqs", "1000", "-v", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(res):
        raise RuntimeError(f"foldseek failed:\n{r.stderr[-800:]}")
    hits = []
    with open(res) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 11:
                continue
            hits.append(dict(target=row[1], fident=float(row[2]), alnlen=int(row[3]),
                             qtm=float(row[4]), ttm=float(row[5]), alntm=float(row[6]),
                             evalue=float(row[7]), prob=float(row[8]),
                             qstart=int(row[9]), qend=int(row[10])))
    return hits


def union_coverage(hits, nres, tm_bar=TM_BAR):
    """Fraction of query residues covered by ANY hit at TM >= tm_bar.

    The published rule gates on ">70% of the sequence covered by domains matching
    a known structure at TM >= 0.5". alnlen was parsed and discarded here for the
    whole project, so the coverage clause was never applied -- a design could fail
    on a high TM over a short segment. This takes the UNION of aligned spans
    rather than summing alnlen, which would double-count overlapping hits.
    """
    spans = sorted((h["qstart"], h["qend"]) for h in hits if h["qtm"] >= tm_bar)
    if not spans or nres <= 0:
        return 0.0
    merged = [list(spans[0])]
    for s, e in spans[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return min(1.0, sum(e - s + 1 for s, e in merged) / nres)


AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def seq_of(pdb, chain=None):
    """One-letter sequence of `chain` (or the only chain) straight from the PDB CA records."""
    out, seen = [], set()
    for line in open(pdb):
        if line[:4] != "ATOM" or line[12:16].strip() != "CA":
            continue
        if chain and line[21] != chain:
            continue
        key = (line[21], line[22:27])
        if key in seen:
            continue
        seen.add(key)
        out.append(AA3to1.get(line[17:20].strip().upper(), "X"))
    return "".join(out)


def novelty_level(tm, seq_id):
    """Adaptyv's general-protein level, verbatim from the blog post. Returns 1-4."""
    hi_struct  = tm >= TM_HIGH
    mod_struct = tm >= TM_MOD
    hi_seq     = seq_id > SEQ_HI
    mid_seq    = seq_id > SEQ_MID
    if hi_seq and mod_struct:
        return 1
    if hi_seq or hi_struct or (mid_seq and mod_struct):
        return 2
    if mid_seq or mod_struct:
        return 3
    return 4


def looks_like_antibody(seq):
    """Crude format flag. The real pipeline uses ANARCI; we only need to know when NOT to
    apply the general-protein rule, because the antibody rule is far more permissive."""
    import re as _re
    motifs = (_re.search(r'W[GRS][QRK]G[TA]', seq or ''), _re.search(r'[YF][YF]C[AGSV]', seq or ''),
              _re.search(r'C[AV][AV]S', seq or ''))
    return sum(1 for m in motifs if m) >= 2


def assess(pdb, db, chain=None, tm_bar=TM_BAR, fid_bar=FID_BAR, cov_bar=COV_BAR):
    ch = chains_in(pdb)
    if not ch:
        return dict(pdb=pdb, error="no ATOM records")
    if chain is None:
        # the binder is the SHORTER chain when a target is present
        chain = min(ch, key=lambda c: ch[c])
    with tempfile.TemporaryDirectory() as wd:
        q = extract_chain(pdb, chain, os.path.join(wd, "q.pdb"))
        hits = search(q, db, wd)
    if not hits:
        return dict(pdb=pdb, chain=chain, nres=ch[chain], hits=0,
                    best_qtm=0.0, best_fident=0.0, best_target="(none)",
                    fident_target="(none)",
                    cov=0.0, passes_struct=True, passes_struct_tmonly=True,
                    passes_seq=True, passes=True)
    top_tm = max(hits, key=lambda h: h["qtm"])
    top_id = max(hits, key=lambda h: h["fident"])
    nres = ch[chain]
    cov = union_coverage(hits, nres, tm_bar)
    # Official rule: structural similarity is disqualifying only when a hit clears
    # the TM bar AND covers more than cov_bar of the chain.
    ps = not (top_tm["qtm"] >= tm_bar and cov > cov_bar)
    # Kept alongside it: the stricter TM-only verdict this project ran all along.
    # Reporting both is what let us state that the coverage clause flips 0 of 20.
    ps_tmonly = top_tm["qtm"] < tm_bar
    pq = top_id["fident"] <= fid_bar
    return dict(pdb=pdb, chain=chain, nres=nres, hits=len(hits),
                best_qtm=top_tm["qtm"], best_target=top_tm["target"],
                best_fident=top_id["fident"], fident_target=top_id["target"],
                cov=cov, passes_struct=ps, passes_struct_tmonly=ps_tmonly,
                passes_seq=pq, passes=ps and pq,
                level=novelty_level(top_tm["qtm"], top_id["fident"]),
                antibody=looks_like_antibody(seq_of(pdb, chain)),
                clears_gate=novelty_level(top_tm["qtm"], top_id["fident"]) >= GATE_LEVEL)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbs", nargs="+")
    ap.add_argument("--db", required=True, help="foldseek database prefix, e.g. fsdb/pdb")
    ap.add_argument("--chain", default=None, help="binder chain (default: shortest)")
    ap.add_argument("--tm-bar", type=float, default=TM_BAR)
    ap.add_argument("--fid-bar", type=float, default=FID_BAR)
    ap.add_argument("--tsv")
    a = ap.parse_args()
    tm_bar, fid_bar = a.tm_bar, a.fid_bar

    rows = []
    print(f"{'design':<36} {'ch':>3} {'res':>4} {'best_TM':>8} {'best_id':>8} {'LVL':>4}  verdict")
    print("-" * 86)
    for p in a.pdbs:
        try:
            r = assess(p, a.db, a.chain, tm_bar, fid_bar)
        except Exception as e:
            print(f"{os.path.basename(p):<36} ERROR {e}")
            continue
        if "error" in r:
            print(f"{os.path.basename(p):<36} {r['error']}")
            continue
        rows.append(r)
        L = r["level"]
        if r.get("antibody"):
            v = "ANTIBODY/UNSCORED -- needs ANARCI CDRH3, antibody rule is more permissive"
        elif L >= GATE_LEVEL:
            v = f"clears gate (level {L})" + ("  [LEVEL 4 de novo]" if L == 4 else "")
        else:
            why = "high structural TM>=0.80" if r["best_qtm"] >= TM_HIGH else "sequence"
            v = f"*** LEVEL {L} -- REJECTED AT UPLOAD ({why})"
        print(f"{os.path.basename(p):<36} {r['chain']:>3} {r['nres']:>4} "
              f"{r['best_qtm']:>8.3f} {r['best_fident']:>8.3f} {r['level']:>4}  {v}  <- {r['best_target'][:22]}")
    if rows:
        import collections
        lv = collections.Counter(r["level"] for r in rows)
        clears = sum(1 for r in rows if r["clears_gate"])
        ab = sum(1 for r in rows if r.get("antibody"))
        print("-" * 96)
        for L in (4, 3, 2, 1):
            if lv[L]: print(f"  level {L}: {lv[L]:>4}")
        print(f"CLEARS THE SUBMISSION GATE (level >= {GATE_LEVEL}): {clears}/{len(rows)} "
              f"= {100*clears/len(rows):.1f}%")
        if lv[2] or lv[1]:
            print(f"  *** {lv[1]+lv[2]} design(s) would be REJECTED AT UPLOAD")
            for r in rows:
                if r["level"] < GATE_LEVEL:
                    print(f"      {os.path.basename(r.get('pdb','?')):<40} TM {r['best_qtm']:.3f} "
                          f"seq {r['best_fident']:.3f} -> level {r['level']}")
        if ab: print(f"  {ab} design(s) flagged ANTIBODY -- scored by the general-protein rule here, "
                     f"which is STRICTER than the antibody rule. Re-check with ANARCI before excluding.")
        print(f"Level 4 (de novo) for reference: {lv[4]}/{len(rows)}; "
              f"Adaptyv round-2 baseline was 2.8%")
    if a.tsv and rows:
        with open(a.tsv, "w") as fh:
            fh.write("design\tchain\tnres\tbest_qtm\tbest_target\tbest_fident\tcov\t"
                     "level\tclears_gate\tantibody\tpasses_struct\tpasses_struct_tmonly\tpasses_seq\tpasses\n")
            for r in rows:
                # The header declares 11 columns. Emit all 11: `cov` and
                # passes_struct_tmonly were missing, which silently SHIFTED every
                # column after best_fident and made the TSV read as its own reversal.
                fh.write(f"{r['pdb']}\t{r['chain']}\t{r['nres']}\t{r['best_qtm']:.4f}\t"
                         f"{r['best_target']}\t{r['best_fident']:.4f}\t{r['cov']:.4f}\t"
                         f"{r['level']}\t{r['clears_gate']}\t{r.get('antibody',False)}\t"
                         f"{r['passes_struct']}\t{r['passes_struct_tmonly']}\t"
                         f"{r['passes_seq']}\t{r['passes']}\n")
        print(f"wrote {a.tsv}")

if __name__ == "__main__":
    main()
