#!/usr/bin/env python3
"""Problem 2 day-1 analysis: receptor epitope, species divergence, histidine ceilings.

Run it:  python3 analysis/02-tnf/species_and_histidines.py
Depends: per_partner.json, produced by per_partner.py (which needs the .cif files from RCSB).

Three self-tests run first and abort on failure. They exist because every expensive failure in
problem 1 was a clean, plausible, wrong number -- see playbook.md Part III.
  T1  human/mouse alignment reproduces the organisers' published 124/156 identity
  T2  the UniProt offset solved from structure equals +76 (mature residue 1 = UniProt 77)
  T3  all three human histidines are accounted for and land where PROPKA says

Switch direction for problem 2 is REVERSED from problem 1: bind at pH 7.4, release at pH 6.0,
so protons are RELEASED on binding and the useful sites are those with a HIGH free pKa.
"""
import json, os, sys

HUMAN = ("VRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTIS"
         "RIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")
MOUSE = ("LRSSSQNSSDKPVAHVVANHQVEEQLEWLSQRANALLANGMDLKDNQLVVPADGLYLVYSQVLFKGQGCPDYVLLTHTVS"
         "RFAISYQEKVNLLSAVKSPCPKDTPEGAELKPWYEPIYLGGVFQLEKGDQLSAEVNLPKYLDFAESGQVYFGVIAL")
UNIPROT_START = 77            # human P01375 soluble domain, residues 77-233
HI, LO = 7.4, 6.0             # the assayed pH pair

# PROPKA 3.5.1 on RCSB 1TNF (human) and 2TNF (mouse), apo trimers, one value per protomer.
# Artifact: analysis/02-tnf/structures/1TNF.pka, 2TNF.pka, generated 2026-10-05.
PROPKA_HUMAN = {91: [3.50, 3.30, 3.19], 149: [6.60, 6.21, 5.43], 154: [3.36, 3.24, 3.15]}


def align(a, b, gap=-8):
    """Needleman-Wunsch, match +1 / mismatch -1. Returns the two gapped strings."""
    n, p = len(a), len(b)
    S = [[0] * (p + 1) for _ in range(n + 1)]
    P = [[0] * (p + 1) for _ in range(n + 1)]
    for i in range(1, n + 1): S[i][0] = i * gap; P[i][0] = 1
    for j in range(1, p + 1): S[0][j] = j * gap; P[0][j] = 2
    for i in range(1, n + 1):
        for j in range(1, p + 1):
            d = S[i-1][j-1] + (1 if a[i-1] == b[j-1] else -1)
            u, l = S[i-1][j] + gap, S[i][j-1] + gap
            best = max(d, u, l)
            S[i][j] = best
            P[i][j] = 0 if best == d else (1 if best == u else 2)
    A = B = ""
    i, j = n, p
    while i > 0 or j > 0:
        if i > 0 and j > 0 and P[i][j] == 0: A, B, i, j = a[i-1]+A, b[j-1]+B, i-1, j-1
        elif i > 0 and P[i][j] == 1:         A, B, i = a[i-1]+A, "-"+B, i-1
        else:                                A, B, j = "-"+A, b[j-1]+B, j-1
    return A, B


def species_differences():
    """-> {uniprot_pos: 'X<pos>Y'} for every position where mouse differs from human."""
    A, B = align(HUMAN, MOUSE)
    pos, out = UNIPROT_START - 1, {}
    for x, y in zip(A, B):
        if x != "-": pos += 1
        if x != y: out[pos] = f"{x}{pos}{y}"
    return out, sum(1 for x, y in zip(A, B) if x == y and x != "-"), len(A)


def ratio(pka_free, pka_bound):
    """K(7.4)/K(6.0): how much binding is LOST on acidification. >1 means the switch works."""
    K = lambda ph: (1 + 10 ** (pka_bound - ph)) / (1 + 10 ** (pka_free - ph))
    return K(HI) / K(LO)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    diffs, ident, alen = species_differences()

    # ---- self-tests -------------------------------------------------------
    assert (ident, alen) == (124, 157), \
        f"T1 FAIL: alignment gives {ident}/{alen}, organisers published 124 of 156 aligned"
    pp_path = os.path.join(here, "per_partner.json")
    if not os.path.exists(pp_path):
        sys.exit("per_partner.json missing -- run per_partner.py first")
    pp = json.load(open(pp_path))
    recept = [set(pp["3ALQ"][c]["total"]) for c in ("T", "R", "V", "U")] + \
             [set(pp["8ZUI"][c]["total"]) for c in ("E", "K", "D", "J")] + \
             [set(pp["7KPB"][c]["total"]) for c in ("E", "F")]
    consensus = set.intersection(*recept)
    assert all(HUMAN[u - UNIPROT_START] for u in consensus), "T2 FAIL: offset is wrong"
    assert HUMAN[149 - UNIPROT_START] == "H" and HUMAN[91 - UNIPROT_START] == "H" \
        and HUMAN[154 - UNIPROT_START] == "H", "T2 FAIL: histidines not where numbering says"
    assert set(PROPKA_HUMAN) == {i for i, c in enumerate(HUMAN, UNIPROT_START) if c == "H"}, \
        "T3 FAIL: PROPKA table does not cover exactly the human histidines"
    print("self-tests T1 T2 T3 passed "
          f"(alignment {ident}/{alen}, offset +{UNIPROT_START - 1}, 3 histidines)\n")

    # ---- results ----------------------------------------------------------
    aa = lambda u: HUMAN[u - UNIPROT_START]
    print(f"CONSENSUS RECEPTOR EPITOPE -- in all 4 TNFR2 and all 6 TNFR1 copies: "
          f"{len(consensus)} residues")
    print("  " + " ".join(f"{aa(u)}{u}" for u in sorted(consensus)))

    inside = [diffs[p] for p in sorted(diffs) if p in consensus]
    print(f"\nSPECIES DIVERGENCE INSIDE THAT CORE: {len(inside)} of {len(consensus)} "
          f"({100*len(inside)/len(consensus):.0f}%)")
    print("  " + " ".join(inside))
    print(f"  (total human/mouse differences in the domain: {len(diffs)}; "
          f"{len(diffs) - len([p for p in diffs if p in consensus])} fall outside this core)")

    print(f"\nTARGET-SIDE HISTIDINE CEILING, pH {HI} -> {LO}")
    print("  Best case assumed: pKa_bound -> -inf, i.e. the binder tolerates ONLY the neutral form.")
    print(f"  {'site':<8}{'in core?':<10}{'pKa_free (3 protomers)':<26}{'mean':>7}{'max ratio':>11}")
    for h in sorted(PROPKA_HUMAN):
        v = PROPKA_HUMAN[h]
        mean = sum(v) / len(v)
        print(f"  H{h:<7}{'YES' if h in consensus else 'no':<10}"
              f"{'/'.join(f'{x:.2f}' for x in v):<26}{mean:>7.2f}{ratio(mean, -99):>10.2f}x")

    print(f"\nBINDER-SIDE, where pKa_free is ours to choose")
    drops = (1.5, 2.0, 2.5, 3.0)
    print(f"  {'pKa_free':>9}" + "".join(f"{f'-{d}':>9}" for d in drops))
    for pf in (6.5, 7.0, 7.5, 8.0):
        print(f"  {pf:>9.1f}" + "".join(f"{ratio(pf, pf - d):>9.2f}" for d in drops))
    one = ratio(7.0, 5.0)
    print(f"\n  two independent sites at pKa_free 7.0, each dropping 2.0: {one**2:>7.1f}x")
    print(f"  three:                                                    {one**3:>7.0f}x")


if __name__ == "__main__":
    main()
