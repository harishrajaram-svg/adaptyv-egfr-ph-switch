#!/usr/bin/env python3
"""Every number inside a CSV `assessment` string must agree with that row's own columns.

WHY THIS GATE EXISTS
The `assessment` column is prose, and prose was the one part of the graded upload with no
mechanical check. It was hand-written when the submission held 10-12 designs and scored at
n = 5, carried through emit_submission_csv.py verbatim ever since, and never re-derived
when the pose sets grew to n = 20/26 and the ranking basis changed from target-only to
his-only. The seven-lens audit of 2026-10-05 found 12 findings in this one column, 9 of
them HIGH, including:

  - a row quoting "3.15x -> 5.46x" as its switch while its own graded column reads 1.023,
    because 5.46 is the SUPERSEDED target-only basis and the string never says so
  - a row calling 3.738x "the all-titratable-site basis" when 3.738 is its HIS-ONLY column
    and its all-site column reads 88.593 -- the exact mislabel SS 11.7 exists to correct
  - the rank-14 row shipping "human 0.594 -> 0.616", "pose spread 1.31x the median" and
    "L133D reads 0.723x", all three of which SS 11.6 formally WITHDRAWS by name after the
    15-seed triad landed, and all three contradicting columns in the row itself
  - two rows both claiming "<0.47 identity to anything else submitted" while being 86.7%
    identical to each other

THREE CHECKS
  BASIS     every `<float>x` token must match one of the row's four pH columns within
            rounding. If it matches ONLY a non-headline basis, the string must name that
            basis nearby, or a reader takes a superseded number for the result.
  COLUMNS   pose count, pose spread and ipSAE values quoted in prose must equal the row's
            own columns.
  WITHDRAWN a literal that SS 11.6 or any "what did not survive" retraction block names as
            withdrawn may not appear in any assessment string.

    check_assessment_strings.py            # report + exit non-zero
    check_assessment_strings.py --verbose  # show every token checked
"""
import csv, json, os, re, sys

CSV = 'submissions/01-egfr.csv'
METHODS = 'submissions/01-egfr-METHODS.md'
TOL = 0.011          # rounding slack for a 2-3 d.p. quote

# A count may be written as a digit or as a word; both must mean the same thing.
WORDS = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six',
         7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten', 11: 'eleven', 12: 'twelve',
         13: 'thirteen', 14: 'fourteen', 15: 'fifteen', 16: 'sixteen', 17: 'seventeen',
         18: 'eighteen', 19: 'nineteen', 20: 'twenty'}


def _eq(got, expected):
    g = str(got).strip().lower()
    try:
        if abs(float(g) - float(expected)) < 1e-9:
            return True
    except ValueError:
        pass
    return g == WORDS.get(int(expected))

BASIS_COLS = {
    'ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE': 'his-only (headline)',
    'ph_ratio_allsite_SENSITIVITY': 'all-site',
    'ph_ratio_partnered_SENSITIVITY': 'partnered',
    'ph_ratio_target_only_SUPERSEDED': 'target-only (SUPERSEDED)',
}
# Tokens that, appearing near a ratio, tell the reader which basis it is on.
BASIS_WORDS = {
    'all-site': ('all-site', 'all-titratable', 'allsite', 'multi-site'),
    'target-only (SUPERSEDED)': ('target-only', 'target only', 'superseded', 'target-side',
                                 'matched wild-type', 'same gate', 'same batch'),
    'partnered': ('partnered', 'counter-charge'),
}


def floats(s):
    """(value, span) for every `<float>x` ratio token."""
    return [(float(m.group(1)), m.span())
            for m in re.finditer(r'(\d+\.\d+)\s*[x×]', s)]


def withdrawn_literals():
    """Numbers a retraction block names as withdrawn -> the claim that was withdrawn.

    Keyed on the BLOCK, not on a quote pattern. The first version matched
    `"<claim>" ... At n = N` and found nothing, because the real text wraps the withdrawn
    claim in markdown emphasis (`*"affinity held: human went **up**, 0.594 -> 0.616."*`)
    and the correction begins on the next line. So: locate the block, split it into its
    numbered items, and in each item take the numbers that appear BEFORE the correction
    marker -- those are the withdrawn ones. Numbers after the marker are the new, correct
    values and must NOT be blocklisted.
    """
    if not os.path.exists(METHODS):
        return {}
    t = open(METHODS).read()
    out = {}
    for bm in re.finditer(r'\*\*What did not survive[^\n]*\*\*(.*?)(?=\n\*\*[A-Z]|\n#{2,})',
                          t, re.S):
        block = bm.group(1)
        for item in re.split(r'\n\s*\d+\.\s', block)[1:]:
            # the correction starts at the first of these; everything before it is the
            # claim being withdrawn
            cut = len(item)
            for marker in ('At n =', 'It is now', 'it is now', 'is **flat'):
                k = item.find(marker)
                if k != -1:
                    cut = min(cut, k)
            claim = item[:cut]
            for mm in re.finditer(r'(\d+\.\d{2,})(?!\s*(?:\u00c5|A\b))', claim):
                v = float(mm.group(1))
                out.setdefault(v, []).append(' '.join(claim.split())[:100])
    # the "Why it ranks" paragraph withdraws the five-pose spread the same way
    for mm in re.finditer(r'pose spread was reported as \*\*([\d.]+)\s*[x\u00d7]\*\*', t):
        out.setdefault(float(mm.group(1)), []).append('pose spread reported on five poses')
    return out


def main():
    verbose = '--verbose' in sys.argv
    rows = list(csv.DictReader(open(CSV)))
    wd = withdrawn_literals()
    bad = []

    # ---- pairwise identity, for the "<0.47 to anything else" style claim --------------
    def ident(a, b):
        if len(a) != len(b):
            import difflib
            sm = difflib.SequenceMatcher(None, a, b)
            return sm.ratio()
        return sum(1 for x, y in zip(a, b) if x == y) / len(a)

    maxid = {}
    for i, r in enumerate(rows):
        best = 0.0
        for j, q in enumerate(rows):
            if i != j:
                best = max(best, ident(r['sequence'], q['sequence']))
        maxid[r['name']] = best

    for i, r in enumerate(rows, 1):
        a = r.get('assessment', '')
        if not a:
            continue
        cols = {}
        for c, label in BASIS_COLS.items():
            try:
                cols[label] = float(r[c])
            except (KeyError, TypeError, ValueError):
                pass
        head = cols.get('his-only (headline)')

        # ---- BASIS ----------------------------------------------------------------
        for v, (s0, s1) in floats(a):
            hits = [lab for lab, cv in cols.items() if abs(cv - v) <= TOL]
            if not hits:
                continue                      # not one of this row's pH numbers
            if 'his-only (headline)' in hits:
                continue                      # it IS the graded number
            ctx = a[max(0, s0 - 160):min(len(a), s1 + 160)].lower()
            named = any(w in ctx for lab in hits for w in BASIS_WORDS.get(lab, ()))
            if not named:
                bad.append((i, r['name'], 'BASIS',
                            f"{v}x matches only {', '.join(hits)} "
                            f"(headline is {head}) and the string does not name the basis"))

        # ---- COLUMNS ---------------------------------------------------------------
        for m in re.finditer(r'(?:over|at n\s*=\s*|across)\s*(\d+)\s*(?:refold\s*)?poses?', a, re.I):
            n = int(m.group(1))
            if r.get('ph_poses_n') and n != int(r['ph_poses_n']):
                bad.append((i, r['name'], 'COLUMNS',
                            f"prose says {n} poses, column ph_poses_n = {r['ph_poses_n']}"))
        for m in re.finditer(r'(\d+\.\d+)x the median', a):
            v = float(m.group(1))
            col = r.get('ph_pose_spread_over_median')
            if col and abs(float(col) - v) > TOL:
                bad.append((i, r['name'], 'COLUMNS',
                            f"prose says spread {v}x the median, column reads {col}"))

        # ---- IDENTITY CLAIM --------------------------------------------------------
        for m in re.finditer(r'identity to any other submitted design:?\s*<?\s*(\d+\.\d+)', a, re.I):
            claim = float(m.group(1))
            real = maxid[r['name']]
            if real > claim + 0.005:
                bad.append((i, r['name'], 'IDENTITY',
                            f"claims <{claim} identity to anything else submitted; "
                            f"measured max is {real:.3f}"))

        # ---- SCOPE ------------------------------------------------------------------
        # The first version of this gate only checked numbers that MATCHED a column, so
        # stale SCOPE sailed through: rank 17 shipped "Six of the twelve submitted rows"
        # long after the submission reached eighteen. A denominator describing the
        # submission must be the submission's size.
        for m in re.finditer(r'(?:of the|the)\s+(' + '|'.join(
                [r'\d+'] + sorted(set(WORDS.values()), key=len, reverse=True)) +
                r')\s+(?:submitted|shipped)\s+(?:rows?|designs?)', a, re.I):
            pre = a[max(0, m.start() - 170):m.start()].lower()
            if any(w in pre for w in ('read ', 'previously', 'earlier', 'withdraw')):
                continue
            if not _eq(m.group(1), len(rows)):
                bad.append((i, r['name'], 'SCOPE',
                            f"says {m.group(1)!r} submitted rows/designs; "
                            f"the submission has {len(rows)}"))

        # ---- BASIS MISLABEL ---------------------------------------------------------
        # A sentence may not name a non-headline basis while quoting the headline value.
        # rank 4 shipped 'On the all-titratable-site basis ... it reads 3.738x' where
        # 3.738 IS its his-only column and its all-site column reads 88.593 -- the
        # headline number presented a second time as an independent multi-site check.
        for m in re.finditer(r'(all-titratable[- ]site|all-site|multi-site)\s*basis'
                             r'[^.]{0,140}?(\d+\.\d+)\s*[x\u00d7]', a, re.I):
            v = float(m.group(2))
            ctx = a[max(0, m.start() - 170):m.start()].lower()
            if any(w in ctx for w in ('read ', 'previously', 'mislabel', 'is not')):
                continue
            if head is not None and abs(v - head) <= TOL:
                allsite = cols.get('all-site')
                bad.append((i, r['name'], 'BASIS-LABEL',
                            f"names the all-site basis and quotes {v}x, which is this "
                            f"row's HIS-ONLY column"
                            + (f" (its all-site column reads {allsite})" if allsite else "")))

        # ---- WITHDRAWN -------------------------------------------------------------
        for v, claims in wd.items():
            # A withdrawn literal that still equals one of THIS row's live column values
            # is not a stale quote -- e.g. 0.594 appears in the withdrawn sentence
            # "human went up, 0.594 -> 0.616" but 0.594 is the parent's correct value.
            # Only 0.616 is withdrawn. Suppressing these keeps the blocklist honest.
            live = any(abs(float(r[c]) - v) <= TOL for c in
                       ('ipsae_min_human', 'ipsae_min_mouse', 'ph_poses_n',
                        'ph_pose_spread_over_median', *BASIS_COLS)
                       if r.get(c) not in (None, ''))
            if live:
                continue
            for m in re.finditer(r'(?<![\d.])' + re.escape(f"{v:g}") + r'(?![\d])', a):
                # A string may NAME a withdrawn figure in order to withdraw it -- that is
                # the correction, not the error. Exempt only when a withdrawal marker sits
                # in the SAME sentence, so "the earlier 0.616 ... is withdrawn" passes
                # while a bare "0.616" does not. Sentence-scoped on purpose: a wider
                # window would let one withdrawal note launder every stale number after it.
                pre = a[:m.start()]
                k = max(pre.rfind('. '), pre.rfind('; '), pre.rfind(': '))
                sent = a[k + 1: a.find('.', m.end()) + 1 or len(a)].lower()
                if any(w in sent for w in ('withdraw', 'earlier', 'superseded',
                                           'was five-pose', 'retract', 'too few poses',
                                           'no longer')):
                    continue
                bad.append((i, r['name'], 'WITHDRAWN',
                            f"{v:g} appears here but METHODS retracts it: \"{claims[0]}\""))
                break

        if verbose:
            print(f"rank {i:>2} {r['name'][:40]:<42} "
                  f"{len(floats(a))} ratio token(s), headline {head}")

    # ---- SUPERLATIVE UNIQUENESS -----------------------------------------------------
    # Two rows both shipped "the closest agreement of any submitted design". A superlative
    # is by definition claimable by one row; two rows asserting it is a contradiction a
    # reviewer can see without leaving the CSV.
    SUPER = ('the closest agreement of any submitted design',
             'the largest in the submission', 'the largest causal swing',
             'the most reproducible', 'the tightest seed reproducibility',
             'the least reproducible row in this submission',
             'the only design', 'the highest pH ratio in the submission')
    for phrase in SUPER:
        holders = []
        for i, r in enumerate(rows, 1):
            a = r.get('assessment', '')
            k = a.lower().find(phrase)
            if k == -1:
                continue
            pre = a[max(0, k - 170):k].lower()
            if any(w in pre for w in ('read ', 'previously', 'earlier', 'both claimed')):
                continue
            holders.append((i, r['name']))
        if len(holders) > 1:
            bad.append((holders[0][0], holders[0][1], 'SUPERLATIVE',
                        f"{len(holders)} rows all claim {phrase!r}: "
                        + ", ".join(f"rank {i} {n}" for i, n in holders)))

    if bad:
        print(f"{'rank':>4}  {'check':<10} design / problem")
        print('-' * 96)
        for i, n, kind, msg in bad:
            print(f"{i:>4}  {kind:<10} {n[:38]}")
            print(f"        {msg}")
        kinds = {}
        for _, _, k, _ in bad:
            kinds[k] = kinds.get(k, 0) + 1
        print(f"\nFAIL: {len(bad)} problem(s) in the graded assessment column "
              f"({', '.join(f'{k} {v}' for k, v in sorted(kinds.items()))}).")
        sys.exit(1)
    print(f"PASS: assessment strings on all {len(rows)} rows agree with their own columns.")


if __name__ == '__main__':
    main()
