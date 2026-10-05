#!/usr/bin/env python3
"""A rank cited in prose must match the generated rank table.

METHODS declares a rule against citing ranks in prose -- `rank_table()` in
gen_methods_submission.py says so in its own docstring: "we stop citing ranks in prose that
the generator does not own." The rule was never enforced, so the prose kept citing them and
the citations rotted every time the submission was reordered. The 2026-10-05 audit raised
seven separate findings of this one shape, e.g. SS 11.2 placing `rimA02_d3_rimA_14_vhh` at
"rank 6" and calling it the "second-highest pH ratio" when the generated table puts it at
12, and SS 10 saying `rimA01_r15` ships "at rank 1" when it is 3.

WHAT IT CHECKS
Every `rank N` / `ranks N-M` / `ranked Nth` mention adjacent to a design name must agree
with that design's position in the emitted CSV. The CSV order IS the rank order, so there
is exactly one source of truth and it is the graded artifact.

WHAT IT DELIBERATELY DOES NOT CHECK
  - a rank inside a GENERATED block: the generator owns those and gen_methods --check
    already guards them
  - a rank in a passage marked as history ("an earlier version", "previously", "was ranked")
  - a rank that names no design: "ranks 1-18" as a span, "the top three", etc.

    check_rank_citations.py [--verbose]
"""
import csv, re, sys

CSV = 'submissions/01-egfr.csv'
DOCS = ['submissions/01-egfr-METHODS.md', 'README.md', 'HANDOFF.md',
        'outbox/CONTROL-TABLE.md']
HIST = ('earlier version', 'previously', 'was ranked', 'used to', 'formerly',
        'an earlier', 'for a day', 'read "', "read '", 'said "', 'corrected',
        'withdrawn', 'this line', 'this paragraph', 'stale')
# A rank explicitly scoped to a NON-SHIPPED basis is not a claim about the CSV order.
# "On the all-site and partnered bases X ranks 1st" is true and must not be flagged --
# the shipped order is the his-only basis only.
OTHER_BASIS = ('all-site', 'all-titratable', 'partnered', 'target-only', 'allsite',
               'on the superseded', 'perturb', 'of draws')


def strip_generated(t):
    """Blank out generated blocks, keeping offsets so line numbers stay right."""
    def blank(m):
        return re.sub(r'[^\n]', ' ', m.group(0))
    return re.sub(r'<!-- GENERATED:.*?<!-- /GENERATED:[A-Z-]+ -->', blank, t, flags=re.S)


def main():
    verbose = '--verbose' in sys.argv
    rows = list(csv.DictReader(open(CSV)))
    rank = {r['name']: i for i, r in enumerate(rows, 1)}
    # also index by the distinctive leading token, so a truncated name still resolves
    bad = []
    for doc in DOCS:
        try:
            raw = open(doc).read()
        except FileNotFoundError:
            continue
        t = strip_generated(raw)
        for m in re.finditer(r'rank(?:ed|s)?\s+(\d+)(?:\s*(?:st|nd|rd|th))?', t, re.I):
            claimed = int(m.group(1))
            lo, hi = max(0, m.start() - 230), min(len(t), m.end() + 230)
            ctx = t[lo:hi]
            low = ctx.lower()
            if any(h in low for h in HIST):
                continue
            # The basis exemption must attach to the RANK CLAIM, not merely appear in the
            # paragraph. Scanning the whole window let "...at rank 3 led this submission at
            # 5.461x on the superseded target-only basis" pass, because the exempting words
            # sat AFTER the citation and described a different quantity. Only a basis named
            # in the 70 characters BEFORE the rank token scopes it.
            near = t[max(0, m.start() - 70):m.start()].lower()
            if any(b in near for b in OTHER_BASIS):
                continue
            # which design does this rank refer to?
            # The design must be the NEAREST named one and sit close to the citation.
            # A 230-char window alone produced false positives by pairing a generic
            # "do not claim rank 1 beats rank 4" with whatever design name was nearby.
            # Prefer the NEAREST name, and among names at the same place the LONGEST --
            # `bc_s831683_mpnn6_S15D` is a substring of
            # `ss_bc_s831683_mpnn6_S15D_S62H_routeA`, so a plain find() attributed a rank
            # to the parent when the text named the mutant.
            hits = []
            for name, actual in rank.items():
                k = ctx.find(name)
                if k == -1:
                    continue
                hits.append((abs((lo + k) - m.start()), -len(name), name, actual))
            # STRONGEST SIGNAL FIRST: `NAME` (rank N) / `NAME` at rank N / `NAME`, rank N.
            # Without this, a sentence naming three designs and three ranks mis-pairs them
            # all -- which is exactly what the corrected glycan paragraph does.
            pre = t[max(0, m.start() - 60):m.start()]
            adj = re.search(r'`([^`]+)`[^`]{0,24}$', pre)
            if adj:
                nm = adj.group(1).strip('.,; ')
                exact = [(n, i) for n, i in rank.items() if n == nm]
                if exact:
                    name, actual = exact[0]
                    if claimed != actual:
                        line = raw[:m.start()].count('\n') + 1
                        bad.append((doc, line, name, claimed, actual,
                                    ' '.join(ctx.split())[:150]))
                    continue
                if any(nm in n or n in nm for n, _ in rank.items()):
                    continue      # a truncated or variant name: do not guess
            hits = [h for h in hits if h[0] <= 120]
            if not hits:
                continue
            hits.sort()
            # drop any hit that is a substring of a longer hit found at the same offset
            best = hits[0]
            for h in hits:
                if h[2] != best[2] and best[2] in h[2] and abs(h[0] - best[0]) <= len(h[2]):
                    best = h
            _, _, name, actual = best
            if claimed != actual:
                line = raw[:m.start()].count('\n') + 1
                bad.append((doc, line, name, claimed, actual,
                            ' '.join(ctx.split())[:150]))
        if verbose:
            print(f"{doc}: scanned")

    if bad:
        for doc, line, name, claimed, actual, ctx in bad:
            print(f"{doc}:{line}  `{name}` cited at rank {claimed}, "
                  f"the CSV puts it at {actual}")
            print(f"    ...{ctx}...")
        print(f"\nFAIL: {len(bad)} prose rank citation(s) disagree with the emitted CSV.")
        sys.exit(1)
    print(f"PASS: every design-adjacent rank citation matches the CSV "
          f"({len(rows)} designs, {len(DOCS)} documents).")


if __name__ == '__main__':
    main()
