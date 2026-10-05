#!/usr/bin/env python3
"""No passage of the external reviewer's emails may be reproduced verbatim in the public tree.

He asked, in writing on 2026-10-05: "Please keep it anonymous and paraphrase the feedback
without my name or initials. No attribution is required." His corrections are load-bearing
throughout this project, so the substance stays and is still credited to an anonymous
reviewer -- but his words are his, and a future edit must not quietly reintroduce them.

Matching is against the emails themselves rather than against quotation marks, because
unquoted reproductions are the ones that slip through: any run of >= 8 consecutive words
appearing in both is flagged.

ALLOWED, and excluded by name below: the competition's own phrasing, quotations of OUR OWN
withdrawn claims (which must appear as written for the correction to mean anything), and
published third-party data he relayed, such as the G532 SPR table, which is from the paper
and not his wording.

The inbox it reads is outside this repository and is not published; when it is absent the
check cannot run and says so rather than passing.
"""
import glob
import os
import re

PUB = '/Users/harish/code/adaptyv-2026/'
INBOX = ('/Users/harish/code/context-directory/projects/'
         'anthropic-adaptyv-2026/inbox/')
FILES = ['submissions/01-egfr-METHODS.md', 'README.md', 'HANDOFF.md',
         'outbox/CONTROL-TABLE.md', 'outbox/PREREGISTRATION.md']
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

total = 0
for rel in FILES:
    s = open(PUB + rel).read()
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
    print(f'{rel}  —  {len(hits)} verbatim run(s), {sum(len(h.split()) for h in hits)} words')
    for h in hits:
        print(f'    [{len(h.split()):>3}w] {h[:150]}')
    print()

if total:
    print(f"\nFAIL: {total} verbatim run(s) of the reviewer's words in the public tree.")
    print('  He asked for paraphrase. Rewrite them, or add a genuinely-not-his phrase')
    print('  to ALLOWED in bin/check_no_verbatim.py with the reason.')
    raise SystemExit(1)
print("PASS: no verbatim reproduction of the reviewer's words in the public tree.")
