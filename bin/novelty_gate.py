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

TM_BAR = 0.50     # "less than moderate" structural similarity
FID_BAR = 0.30    # sequence identity ceiling
COV_BAR = 0.70    # ">70% of the sequence covered" -- the clause this script used to ignore
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
                passes_seq=pq, passes=ps and pq)


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
    print(f"{'design':<36} {'ch':>3} {'res':>4} {'best_TM':>8} {'best_id':>8}  verdict")
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
        v = ("PASS" if r["passes"] else
             "FAIL struct+seq" if not r["passes_struct"] and not r["passes_seq"] else
             "FAIL structure" if not r["passes_struct"] else "FAIL sequence")
        print(f"{os.path.basename(p):<36} {r['chain']:>3} {r['nres']:>4} "
              f"{r['best_qtm']:>8.3f} {r['best_fident']:>8.3f}  {v}  <- {r['best_target'][:26]}")
    if rows:
        n = sum(1 for r in rows if r["passes"])
        ns = sum(1 for r in rows if r["passes_struct"])
        print("-" * 86)
        print(f"structural (TM < {tm_bar}):      {ns}/{len(rows)} = {100*ns/len(rows):.1f}%")
        print(f"BOTH halves (Level 4):       {n}/{len(rows)} = {100*n/len(rows):.1f}%")
        print(f"round-2 baseline for comparison: 3.6% structural, 2.8% full Level 4")
    if a.tsv and rows:
        with open(a.tsv, "w") as fh:
            fh.write("design\tchain\tnres\tbest_qtm\tbest_target\tbest_fident\tcov\tpasses_struct\tpasses_struct_tmonly\tpasses_seq\tpasses\n")
            for r in rows:
                # The header declares 11 columns. Emit all 11: `cov` and
                # passes_struct_tmonly were missing, which silently SHIFTED every
                # column after best_fident and made the TSV read as its own reversal.
                fh.write(f"{r['pdb']}\t{r['chain']}\t{r['nres']}\t{r['best_qtm']:.4f}\t"
                         f"{r['best_target']}\t{r['best_fident']:.4f}\t{r['cov']:.4f}\t"
                         f"{r['passes_struct']}\t{r['passes_struct_tmonly']}\t"
                         f"{r['passes_seq']}\t{r['passes']}\n")
        print(f"wrote {a.tsv}")

if __name__ == "__main__":
    main()
