#!/usr/bin/env python3
"""Every live surface must state the SAME limitation count as section 14 actually has.

WHY THIS EXISTS. The lesson "a number fixed in one document is not fixed" has now been
learned three times in this project: a +20% effect size corrected in METHODS section 11
and not in the public methodology box (2026-10-07 morning), then corrected in one of the
box's two prose fields and not the other (found that night), and a limitation count that
went 48 -> 53 while two other documents still said 48. Each time it was written down as a
lesson. Writing it down did not work. This is the check.

WHAT IS A LIVE CLAIM, AND WHAT IS HISTORY. "48 limitations" in a sentence recording what a
past change did ("43 -> 48 limitations") is CORRECT as history and must not be rewritten --
this project does not edit its own record. So a match is exempt when a transition arrow
appears before it on the same line.

Cross-repo by design: the published box and the roadmap live in the notes vault, not here.
A missing vault is SKIPPED and named, never silently passed, so a clone stays green.

    check_published_counts.py              the gate
    check_published_counts.py --selftest   self-test, including two mutation tests
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAULT = os.path.join(ROOT, "..", "context-directory", "projects", "anthropic-adaptyv-2026")

METHODS = os.path.join(ROOT, "submissions", "02-tnf-METHODS.md")
SURFACES = [
    METHODS,
    os.path.join(VAULT, "outbox", "02-proteinbase-methodology-box.md"),
    os.path.join(VAULT, "ROADMAP.md"),
]

# "53 limitations", "limitations ... 53", "**53 as of"
CLAIM = re.compile(r"(\d{2,3})\s+limitations|limitations[^.\n]{0,40}?\*\*(\d{2,3})\b")
ARROW = re.compile(r"(->|→)")


def count_limitations(text):
    """Numbered items under '## 14. Limitations', up to the next '## ' heading."""
    i = text.find("## 14. Limitations")
    if i < 0:
        return None
    rest = text[i + 1:]
    j = rest.find("\n## ")
    sec = rest[:j] if j >= 0 else rest
    nums = [int(m.group(1)) for m in re.finditer(r"^(\d{1,3})\.\s+\*\*", sec, re.M)]
    return max(nums) if nums else None


def live_claims(text):
    """[(lineno, claimed_int)] for live claims only; history (a transition arrow) is exempt."""
    out = []
    for n, line in enumerate(text.splitlines(), 1):
        for m in CLAIM.finditer(line):
            if ARROW.search(line[:m.start()]):
                continue                      # "43 -> 48 limitations": a record, not a claim
            v = m.group(1) or m.group(2)
            out.append((n, int(v)))
    return out


def check():
    if not os.path.exists(METHODS):
        print(f"SKIP  {METHODS} absent"); return 0
    truth = count_limitations(open(METHODS).read())
    if truth is None:
        print("FAIL  could not count section 14's limitations"); return 1
    print(f"section 14 has {truth} numbered limitations")
    bad = 0
    for f in SURFACES:
        if not os.path.exists(f):
            print(f"  SKIP   {os.path.relpath(f, ROOT)} -- not present (clone without the vault)")
            continue
        claims = live_claims(open(f).read())
        if not claims:
            print(f"  ok     {os.path.relpath(f, ROOT)} -- states no count")
            continue
        for ln, v in claims:
            if v == truth:
                print(f"  ok     {os.path.relpath(f, ROOT)}:{ln} says {v}")
            else:
                print(f"  FAIL   {os.path.relpath(f, ROOT)}:{ln} says {v}, section 14 has {truth}")
                bad += 1
    print("PASS: every live surface agrees" if not bad else f"FAIL: {bad} stale count(s)")
    return 1 if bad else 0


def selftest():
    doc = ("## 14. Limitations\n\n1. **a** x\n2. **b** y\n3. **c** z\n\n## 15. Next\n")
    assert count_limitations(doc) == 3, count_limitations(doc)
    # a heading-less document must refuse rather than report zero
    assert count_limitations("no section here") is None
    # the next '## ' heading bounds the section, so later numbered lists are not counted
    assert count_limitations(doc + "\n9. **not a limitation** q\n") == 3
    # live claim found
    assert live_claims("the attachment carries 53 limitations\n") == [(1, 53)]
    # the box's actual phrasing, not a paraphrase -- a '.' between the two halves would
    # break the pattern, so the test uses the real line
    assert live_claims(
        "- [ ] number of limitations in the attached document \u2014 **53 as of 2026-10-07**\n"
    ) == [(1, 53)]
    # MUTATION 1: history must stay exempt, or this gate would force rewriting the record
    assert live_claims("- §14 -- 43 -> 48 limitations, after C+\n") == []
    assert live_claims("- §14 — 43 → **48 limitations**\n") == []
    # MUTATION 2: a stale live claim must be CAUGHT, not shrugged at
    assert live_claims("four retired instruments, 48 limitations, a fresh clone\n") == [(1, 48)]
    # and a one-digit number is not a count (avoids matching "§4 limitations")
    assert live_claims("see 4 limitations\n") == []
    print("check_published_counts.py --selftest PASS")
    return 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else check())
