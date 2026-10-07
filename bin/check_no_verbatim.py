#!/usr/bin/env python3
"""No passage of the external reviewer's emails may be reproduced verbatim in the public tree.

The reviewer asked, in writing on 2026-10-05, that the feedback be carried anonymously and in
our own words, with no name or initials and no attribution. Their corrections are load-bearing
throughout this project, so the substance stays and is still credited to an anonymous reviewer --
but the wording is theirs, and a future edit must not quietly reintroduce it.

Matching is against the emails themselves rather than against quotation marks, because
unquoted reproductions are the ones that slip through: any run of >= 8 consecutive words
appearing in both is flagged.

SCOPE IS THE WHOLE PUBLISHED TREE, from `git ls-files`, not a hand-kept list. The hand-kept
list of five files omitted `submissions/01-egfr.csv` -- the one artifact that is actually
uploaded -- and that file carried a quoted passage through every green run of this gate. A
fix is not applied until every consumer of the broken artifact is re-pointed at the fixed one.

This gate ALSO refuses the reviewer's initials anywhere in the published tree, which is the
other half of what was asked and was not checked at all.

ALLOWED, and excluded by name below: the competition's own phrasing, quotations of OUR OWN
withdrawn claims (which must appear as written for the correction to mean anything), and
published third-party data they relayed, such as the G532 SPR table, which is from the paper
and not their wording.

The inbox it reads is outside this repository and is not published; when it is absent the
check cannot run and says so rather than passing.
"""
import glob
import os
import re
import subprocess

PUB = '/Users/harish/code/adaptyv-2026/'
INBOX = ('/Users/harish/code/context-directory/projects/'
         'anthropic-adaptyv-2026/inbox/')
# Every published text file, so a new document cannot be born outside the gate's scope.
SKIP_EXT = {'.pdb', '.cif', '.npz', '.pt', '.png', '.jpg', '.gz', '.zip', '.a3m'}
FILES = sorted(
    f for f in subprocess.run(['git', 'ls-files'], cwd=PUB, capture_output=True,
                              text=True, check=True).stdout.split()
    if os.path.splitext(f)[1].lower() not in SKIP_EXT
)
# The reviewer's initials, which were asked to be kept out of the published tree entirely.
# Word-boundary matched so pKa, PROPKA and similar do not trip it.
# Built from character codes so that this gate does not itself publish the initials it
# forbids. _I1/_I2 are the two letters; spelling them out here would fail the check below.
_I1, _I2 = chr(80), chr(75)
INITIALS = re.compile(
    r"\b%s%s\b|\b%s\.\s?%s\." % (_I1, _I2, _I1, _I2))
MIN_RUN = 8


def norm(t):
    t = t.replace('’', "'").replace('“', '"').replace('”', '"')
    t = t.replace('—', ' ').replace('–', ' ').replace('§', 'section')
    return re.sub(r'[^a-z0-9 ]+', ' ', t.lower()).split()


his = []
for f in sorted(glob.glob(INBOX + '*.md')):
    his += norm(open(f).read())
hisgrams = set()
for i in range(len(his) - MIN_RUN + 1):
    hisgrams.add(' '.join(his[i:i + MIN_RUN]))
print(f'reviewer corpus: {len(his):,} words, {len(hisgrams):,} distinct {MIN_RUN}-grams\n')

ALLOWED = (
    'at ph 6 5 than at 7 4',                              # the competition's own wording
    'the 95th percentile of a matched calibration null',   # our own withdrawn claim
    'target kd at ph 6 5 kd at ph 7 4',                    # the published G532 SPR table
    'human egfr 294 nm 3 900 nm',
)

if not os.path.isdir(INBOX):
    print(f'SKIP: reviewer corpus not found at {INBOX} -- cannot verify.')
    raise SystemExit(0)

named = []
for rel in FILES:
    try:
        body = open(PUB + rel, encoding='utf-8').read()
    except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError, OSError):
        continue
    for m in INITIALS.finditer(body):
        line = body.count('\n', 0, m.start()) + 1
        named.append((rel, line, body[max(0, m.start() - 60):m.end() + 60].replace('\n', ' ')))
if named:
    print(f'{len(named)} use(s) of the reviewer\'s initials in the published tree:')
    for rel, line, ctx in named:
        print(f'    {rel}:{line}  ...{ctx.strip()}...')
    print()
else:
    print('initials: none in the published tree.\n')

total = 0
for rel in FILES:
    try:
        s = open(PUB + rel, encoding='utf-8').read()
    except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError, OSError):
        continue
    w = norm(s)
    hits, i = [], 0
    while i < len(w) - MIN_RUN + 1:
        if ' '.join(w[i:i + MIN_RUN]) in hisgrams:
            j = i
            while j < len(w) - MIN_RUN + 1 and ' '.join(w[j:j + MIN_RUN]) in hisgrams:
                j += 1
            run = ' '.join(w[i:j + MIN_RUN - 1])
            if not any(a in run for a in ALLOWED):
                hits.append(run)
            i = j
        else:
            i += 1
    total += len(hits)
    if hits:
        print(f'{rel}  —  {len(hits)} verbatim run(s), '
              f'{sum(len(h.split()) for h in hits)} words')
        for h in hits:
            print(f'    [{len(h.split()):>3}w] {h[:150]}')
        print()

if total or named:
    if total:
        print(f"FAIL: {total} verbatim run(s) of the reviewer's words in the public tree.")
        print('  Paraphrase them, or add a phrase that is genuinely not theirs')
        print('  to ALLOWED in bin/check_no_verbatim.py with the reason.')
    if named:
        print(f'FAIL: {len(named)} use(s) of the reviewer\'s initials. Anonymise them.')
    raise SystemExit(1)
print(f"PASS: {len(FILES)} published files carry neither the reviewer's words "
      'nor their initials.')
