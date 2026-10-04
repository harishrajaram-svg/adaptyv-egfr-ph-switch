#!/usr/bin/env python3
"""Adaptyv's ANTIBODY novelty rule, applied to our VHH-format designs.

WHY THIS FILE EXISTS
--------------------
`bin/novelty_gate.py` scores every design with the GENERAL-PROTEIN rule. For an
antibody-format design that rule is the wrong one and it is far stricter: it
disqualifies on global sequence identity, which for a VHH is mostly framework --
shared with every other VHH on earth by construction, and not what Adaptyv
measure. Their antibody branch gates on CDRH3 instead:

    L1  CDRH3 >= 95%,  OR  CDRH3 >= 70% AND global >= 95%
    L2  CDRH3 70-95% AND global < 95%
    L3  CDRH3 < 70%  AND global >= 70%     <- clears the submission gate
    L4  CDRH3 < 70%  AND global < 70%

https://www.adaptyvbio.com/blog/novelty

This cost us a real candidate: `rimA02/d3_rimA_14` was dropped on 77.5% global
identity when under the antibody rule it is Level 3 and eligible -- and it carries
the strongest pH switch in the project (5.55x).

WHAT THIS IS NOT
----------------
Adaptyv run ANARCI numbering and MMseqs2 against the therapeutic-antibody database
and PLAbDab. We have neither. So:
  * CDRH3 is delimited by the conserved framework motifs -- the last Cys of FR3 and
    the W of the FR4 W[GRS][QRK]G motif -- not by IMGT numbering. On a well-formed
    VHH these coincide; on a malformed one this script says UNDELIMITED rather than
    guessing.
  * The comparison set is PDB via FoldSeek, not an antibody database. Our identities
    are LOWER BOUNDS: a design clean here can still hit a patent or a PLAbDab entry.
  * Identity is computed over the structural (TMalign) alignment FoldSeek returns,
    per region, from qaln/taln -- so the global and CDRH3 numbers come off the same
    alignment and are directly comparable to each other.

Both verdicts are printed. A design is reported eligible only if it clears under the
antibody rule AND is genuinely antibody-format; the general-protein level is shown
alongside so the cost of the format call is visible.

USAGE
    python3 bin/antibody_novelty.py design.pdb [...] --db tools/fsdb/pdb [--tsv out.tsv]
    python3 bin/antibody_novelty.py --selftest
"""
import argparse, csv, os, re, subprocess, sys, tempfile

FOLDSEEK = os.environ.get("FOLDSEEK", "tools/foldseek/bin/foldseek")
FMT = "query,target,fident,alnlen,qtmscore,qaln,taln,qstart,qend"
CDRH3_HI, CDRH3_LO, GLOBAL_HI, GLOBAL_LO = 0.95, 0.70, 0.95, 0.70
GATE_LEVEL = 3

AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
FR4 = re.compile(r'W[GRS][QRKE]G[TAS]')


def chains_in(pdb):
    out = {}
    for line in open(pdb):
        if line[:4] == "ATOM":
            out.setdefault(line[21], set()).add(line[22:27])
    return {c: len(v) for c, v in out.items()}


def seq_of(pdb, chain):
    out, seen = [], set()
    for line in open(pdb):
        if line[:4] != "ATOM" or line[12:16].strip() != "CA" or line[21] != chain:
            continue
        key = line[22:27]
        if key in seen:
            continue
        seen.add(key)
        out.append(AA3to1.get(line[17:20].strip().upper(), "X"))
    return "".join(out)


def extract_chain(pdb, chain, dest):
    with open(dest, "w") as fh:
        for line in open(pdb):
            if line[:4] == "ATOM" and line[21] == chain:
                fh.write(line)
        fh.write("END\n")
    return dest


def is_antibody(seq):
    """Two of three VHH framework motifs. Deliberately the same crude test
    novelty_gate.py uses, so the two scripts never disagree about FORMAT."""
    motifs = (FR4.search(seq or ''), re.search(r'[YF][YF]C[AGSVT]', seq or ''),
              re.search(r'C[AV][AV][SKTR]', seq or ''))
    return sum(1 for m in motifs if m) >= 2


def cdrh3_span(seq):
    """(start, end) 0-based half-open span of CDRH3, or None if undelimitable.

    CDRH3 runs from just after the conserved Cys that closes FR3 to just before
    the W of the FR4 WGxG motif. We take the LAST FR4 match (a VHH has one, but a
    scFv has two and the heavy chain is the later one) and the last Cys before it.
    A span outside 1..40 residues means the motifs are not really framework and we
    refuse to report a number rather than invent one.
    """
    if not seq:
        return None
    m = None
    for m in FR4.finditer(seq):
        pass
    if m is None:
        return None
    w = m.start()
    c = seq.rfind('C', 0, w)
    if c < 0:
        return None
    s, e = c + 1, w
    if not (1 <= e - s <= 40):
        return None
    return s, e


def search(query_pdb, db, workdir):
    res = os.path.join(workdir, "res.tsv")
    tmp = os.path.join(workdir, "tmp")
    os.makedirs(tmp, exist_ok=True)
    r = subprocess.run([FOLDSEEK, "easy-search", query_pdb, db, res, tmp,
                        "--alignment-type", "1", "--format-output", FMT,
                        "--max-seqs", "1000", "-v", "1"], capture_output=True, text=True)
    if not os.path.exists(res):
        raise RuntimeError(f"foldseek failed:\n{r.stderr[-800:]}")
    hits = []
    for row in csv.reader(open(res), delimiter="\t"):
        if len(row) < 9:
            continue
        hits.append(dict(target=row[1], fident=float(row[2]), alnlen=int(row[3]),
                         qtm=float(row[4]), qaln=row[5], taln=row[6],
                         qstart=int(row[7]), qend=int(row[8])))
    return hits


def region_identity(hit, span):
    """Identity over the whole alignment and over `span` (query coordinates).

    Walks qaln/taln together, tracking the query index so a CDRH3 that is partly
    unaligned is scored over the residues it actually has. Denominator for the
    region is the span length, NOT the aligned columns -- an unaligned CDRH3 is
    dissimilar, not unmeasured.
    """
    qi = hit["qstart"] - 1
    g_match = g_cols = 0
    r_match = 0
    for qc, tc in zip(hit["qaln"], hit["taln"]):
        if qc != '-' and tc != '-':
            g_cols += 1
            if qc == tc:
                g_match += 1
                if span and span[0] <= qi < span[1]:
                    r_match += 1
        if qc != '-':
            qi += 1
    g = g_match / g_cols if g_cols else 0.0
    r = (r_match / (span[1] - span[0])) if span else None
    return g, r


def antibody_level(cdrh3, glob):
    if cdrh3 >= CDRH3_HI or (cdrh3 >= CDRH3_LO and glob >= GLOBAL_HI):
        return 1
    if cdrh3 >= CDRH3_LO:
        return 2
    return 3 if glob >= GLOBAL_LO else 4


def assess(pdb, db, chain=None):
    ch = chains_in(pdb)
    if not ch:
        return dict(pdb=pdb, error="no ATOM records")
    if chain is None:
        chain = min(ch, key=lambda c: ch[c])      # binder = shorter chain
    seq = seq_of(pdb, chain)
    ab = is_antibody(seq)
    span = cdrh3_span(seq)
    with tempfile.TemporaryDirectory() as wd:
        hits = search(extract_chain(pdb, chain, os.path.join(wd, "q.pdb")), db, wd)
    if not hits:
        return dict(pdb=pdb, chain=chain, nres=ch[chain], seq=seq, antibody=ab,
                    span=span, hits=0, level=4, clears=True, note="no PDB hit")
    # The rule asks about the MOST SIMILAR antibody, so rank candidates by global
    # identity over the alignment -- not by TM, which is near-1.0 for every
    # immunoglobulin fold and cannot discriminate.
    scored = []
    for h in hits:
        g, r = region_identity(h, span)
        scored.append((g, r, h))
    scored.sort(key=lambda t: -t[0])
    g, r, h = scored[0]
    out = dict(pdb=pdb, chain=chain, nres=ch[chain], seq=seq, antibody=ab, span=span,
               hits=len(hits), best_target=h["target"], best_qtm=h["qtm"],
               global_id=g, fident_fs=h["fident"], cdrh3_id=r,
               cdrh3=seq[span[0]:span[1]] if span else None)
    if r is None:
        out.update(level=None, clears=None, note="CDRH3 UNDELIMITED -- motifs absent")
    else:
        L = antibody_level(r, g)
        out.update(level=L, clears=L >= GATE_LEVEL, note="")
    return out


SELFTEST = [
    # (seq, expected CDRH3)  -- a real VHH, and the shapes that must NOT be scored
    ("QVQLQESGGGLVQAGGSLRLSCAASGRTFSEYAMGWFRQAPGKEREFVATISWSGGSTYYADSVKGRFTISRDNAKNTVYL"
     "QMNSLKPEDTAVYYCAAGSRFSSYWGQGTQVTVSS", "AAGSRFSSY"),   # IMGT CDR3 = (C104, W118)
    #  ^ the span INCLUDES the two residues after the conserved Cys. My first expectation
    #    here was "GSRFSSY" (a Kabat-style start) and the code was right, not the test.
    ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", None),            # no framework at all
]


def selftest():
    ok = True
    for seq, want in SELFTEST:
        span = cdrh3_span(seq)
        got = seq[span[0]:span[1]] if span else None
        flag = "ok " if got == want else "FAIL"
        if got != want:
            ok = False
        print(f"  {flag} cdrh3_span -> {got!r} (want {want!r})")
    for cdrh3, glob, want in [(0.96, 0.10, 1), (0.80, 0.96, 1), (0.80, 0.50, 2),
                              (0.10, 0.75, 3), (0.10, 0.50, 4), (0.136, 0.705, 3)]:
        got = antibody_level(cdrh3, glob)
        flag = "ok " if got == want else "FAIL"
        if got != want:
            ok = False
        print(f"  {flag} level(cdrh3={cdrh3}, global={glob}) -> {got} (want {want})")
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbs", nargs="*")
    ap.add_argument("--db", default="tools/fsdb/pdb")
    ap.add_argument("--chain", default=None)
    ap.add_argument("--tsv")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.pdbs:
        ap.error("give at least one pdb, or --selftest")
    rows = []
    for p in a.pdbs:
        try:
            r = assess(p, a.db, a.chain)
        except Exception as e:
            print(f"{os.path.basename(p):<34} ERROR {e}")
            continue
        if "error" in r:
            print(f"{os.path.basename(p):<34} {r['error']}")
            continue
        rows.append(r)
        name = os.path.basename(p)[:-4]
        print(f"\n=== {name}")
        print(f"  chain {r['chain']}  {r['nres']} res   antibody-format: {r['antibody']}")
        print(f"  closest PDB chain by sequence : {r.get('best_target','-')}  (TM {r.get('best_qtm',0):.3f})")
        print(f"  global identity to it         : {r.get('global_id',0)*100:.1f}%")
        if r.get("cdrh3") is not None:
            print(f"  CDRH3 ({len(r['cdrh3'])} res)              : {r['cdrh3']}")
            print(f"  CDRH3 identity                : {r['cdrh3_id']*100:.1f}%")
        if r["level"] is None:
            print(f"  ANTIBODY LEVEL                : -- {r['note']}")
        else:
            v = "CLEARS the gate" if r["clears"] else "*** REJECTED AT UPLOAD"
            print(f"  ANTIBODY LEVEL                : {r['level']}  {v}")
        if not r["antibody"]:
            print("  NOTE: format flag is FALSE -- the antibody rule may not apply. "
                  "Score with the general-protein rule.")
    if a.tsv and rows:
        with open(a.tsv, "w") as fh:
            fh.write("design\tchain\tnres\tantibody\tbest_target\tbest_qtm\tglobal_id\t"
                     "cdrh3\tcdrh3_id\tlevel\tclears_gate\tnote\n")
            for r in rows:
                cid = r.get("cdrh3_id")
                cid = "" if cid is None else f"{cid:.4f}"
                fh.write(f"{r['pdb']}\t{r['chain']}\t{r['nres']}\t{r['antibody']}\t"
                         f"{r.get('best_target','')}\t{r.get('best_qtm',0):.4f}\t"
                         f"{r.get('global_id',0):.4f}\t{r.get('cdrh3') or ''}\t"
                         f"{cid}\t{r['level']}\t{r['clears']}\t{r['note']}\n")
        print(f"\nwrote {a.tsv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
