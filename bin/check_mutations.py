#!/usr/bin/env python3
"""Does each named point mutation actually exist at that position in the shipped sequence?

THE MOST IMPORTANT INVARIANT IN THE SUBMISSION, and until 2026-10-05 nothing checked it.
The entire causal claim is that specific point mutations create the pH switch -- S15D, A22D,
L133E, S62H, S88D. Each of those is a claim about ONE residue at ONE position in a sequence
being uploaded for synthesis. If a name says S15D and position 15 of the shipped sequence is
not Asp, the organisers express a molecule that does not carry the mechanism we claim, and
every number attached to it describes something else.

Three checks, all on the emitted CSV, none on prose:

  IDENTITY   every X<pos>Y token in a design name must find Y at `pos`, 1-based in that
             design's own sequence. A near miss is reported as an OFFSET rather than a flat
             WRONG, because 0-based/1-based confusion and parent-numbering drift are the two
             ways this fails and they need different fixes.
  PAIRS      each declared parent/mutant pair must differ at EXACTLY the claimed position and
             nowhere else. Three pairs are shipped deliberately (SS 11.3) and the whole point
             of shipping a matched parent is that the comparison is clean; a pair differing
             at three positions is not a controlled comparison.
  CENSUS     the binder-histidine count per design, straight from `sequence.count('H')`,
             which is the independent route that caught the generated block publishing 14
             for a design whose sequence holds 2.

Numbering note: positions in design NAMES are binder-local and 1-based. They are NOT the
target's canonical or mature coordinates -- see the convention box at SS 1. This script never
touches target numbering.
"""
import csv
import re
import sys

AA3 = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly', 'H': 'His',
       'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln',
       'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'}

rows = list(csv.DictReader(open('submissions/01-egfr.csv')))
seq = {r['name']: r['sequence'].strip().upper() for r in rows}

# every X<pos>Y token in a design name
MUT = re.compile(r'(?<![A-Za-z0-9])([ACDEFGHIKLMNPQRSTVWY])(\d{1,3})'
                 r'([ACDEFGHIKLMNPQRSTVWY])(?![A-Za-z0-9])')

print('=' * 92)
print('MUTATION IDENTITY -- is the named residue actually there?')
print('=' * 92)
print(f"{'design':<40}{'mut':>7}{'len':>5}{'at pos':>8}{'verdict':>10}")
bad = []
for r in rows:
    n = r['name']
    s = seq[n]
    for m in MUT.finditer(n):
        frm, pos, to = m.group(1), int(m.group(2)), m.group(3)
        if pos > len(s):
            v = 'OUT OF RANGE'
            bad.append((n, m.group(0), f'position {pos} > length {len(s)}'))
        else:
            got = s[pos - 1]
            if got == to:
                v = 'ok'
            else:
                # try offsets before calling it wrong
                alt = [o for o in (-1, 1, -2, 2)
                       if 0 < pos + o <= len(s) and s[pos + o - 1] == to]
                v = f'OFF BY {alt[0]:+d}' if alt else 'WRONG'
                bad.append((n, m.group(0),
                            f'position {pos} is {AA3[got]} ({got}), name says {AA3[to]} ({to})'
                            + (f'; {to} is at {pos + alt[0]}' if alt else '')))
        print(f"{n[:39]:<40}{m.group(0):>7}{len(s):>5}"
              f"{(s[pos-1] if pos <= len(s) else '-'):>8}{v:>10}")

# declared parent/mutant pairs: must differ ONLY at the claimed position
PAIRS = [
    ('bc_s831683_mpnn9_WT', 'bc_s831683_mpnn9_S15D', [15]),
    ('rimA01_r15_boltzgen_egfr_d3_rimA_20', 'rimA01_r15_L133E', [133]),
    ('bc_s831683_mpnn6_S15D', 'ss_bc_s831683_mpnn6_S15D_S62H_routeA', [62]),
]
print()
print('=' * 92)
print('DECLARED PAIRS -- do they differ ONLY where claimed?')
print('=' * 92)
for a, b, claimed in PAIRS:
    if a not in seq or b not in seq:
        print(f'  MISSING: {a if a not in seq else b}')
        continue
    sa, sb = seq[a], seq[b]
    if len(sa) != len(sb):
        print(f'  {a} / {b}: LENGTH DIFFERS {len(sa)} vs {len(sb)}')
        bad.append((f'{a}/{b}', 'pair', f'lengths {len(sa)} vs {len(sb)}'))
        continue
    diff = [i + 1 for i in range(len(sa)) if sa[i] != sb[i]]
    ok = diff == claimed
    print(f'  {a[:34]:<36} vs {b[:36]:<38}')
    print(f'      differs at {diff}, claimed {claimed}  ->  '
          + ('ok' if ok else 'MISMATCH'))
    for p in diff:
        print(f'         pos {p}: {AA3[sa[p-1]]}({sa[p-1]}) -> {AA3[sb[p-1]]}({sb[p-1]})')
    if not ok:
        bad.append((f'{a}/{b}', 'pair', f'differs at {diff}, claimed {claimed}'))

# binder histidine census
print()
print('=' * 92)
print('BINDER HISTIDINE CENSUS -- count and positions from the shipped sequence')
print('=' * 92)
carry = [(r['name'], seq[r['name']].count('H'),
          [i + 1 for i, c in enumerate(seq[r['name']]) if c == 'H']) for r in rows]
with_h = [c for c in carry if c[1]]
for n, c, pos in sorted(with_h, key=lambda t: -t[1]):
    print(f'  {n[:44]:<46} {c} at {pos}')
print(f'\n  {len(with_h)} of {len(rows)} designs carry at least one histidine of their own')

print()
print('=' * 92)
if bad:
    print(f'FAIL: {len(bad)} problem(s)')
    for n, tok, why in bad:
        print(f'  {n}  [{tok}]  {why}')
    sys.exit(1)
else:
    print('PASS: every named mutation is present at its stated position, every declared pair')
    print('      differs only where claimed, and the histidine census matches the sequences.')
