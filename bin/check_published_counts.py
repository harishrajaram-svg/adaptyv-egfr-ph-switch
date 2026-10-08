#!/usr/bin/env python3
"""Every live surface must state the SAME counts the code actually has.

Two counts are checked: the number of limitations in METHODS section 14, and the number of
gates in bin/gate_sweep.py's GATES list. Both have gone stale in a document before.

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
import ast
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAULT = os.path.join(ROOT, "..", "context-directory", "projects", "anthropic-adaptyv-2026")

METHODS = os.path.join(ROOT, "submissions", "02-tnf-METHODS.md")
GATE_SWEEP = os.path.join(ROOT, "bin", "gate_sweep.py")
SURFACES = [
    METHODS,
    os.path.join(ROOT, "HANDOFF.md"),
    os.path.join(VAULT, "outbox", "02-proteinbase-methodology-box.md"),
    os.path.join(VAULT, "ROADMAP.md"),
]

# "53 limitations", "limitations ... 53", "**53 as of"
CLAIM = re.compile(r"(\d{2,3})\s+limitations|limitations[^.\n]{0,40}?\*\*(\d{2,3})\b")
# "23 gates", "23 green gates", "23 of 23 gates". The "N of N" form is deliberately the only
# accepted way to write the number twice on one line, so both halves get checked.
# The negative lookahead is because "gates" is also a verb: "item 10 gates every writing item" is
# not a count, and without it that sentence makes this check red.
GATE_CLAIM = re.compile(r"(\d{1,3})\s+of\s+(\d{1,3})\s+gates\b"
                        r"|(\d{1,3})\s+(?:green\s+)?gates\b(?!\s+(?:every|all|each|the|this|that|it|them|us))")
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


def count_gates(path=GATE_SWEEP):
    """How many gates gate_sweep.py actually RUNS -- by parsing, never by importing it.

    That is len(GATES) PLUS every gate appended under a literal name in main(). `references`
    is appended there rather than listed, so counting the list alone undercounts by one --
    which is how "22" was written into two documents while 23 gates were running.
    """
    if not os.path.exists(path):
        return None
    tree = ast.parse(open(path).read())
    listed = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
                getattr(t, "id", None) == "GATES" for t in node.targets):
            if isinstance(node.value, (ast.List, ast.Tuple)):
                listed = len(node.value.elts)
    if listed is None:
        return None
    appended = 0
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute) and node.func.attr == "append"
                and getattr(node.func.value, "id", None) == "results"
                and node.args and isinstance(node.args[0], ast.Tuple)
                and node.args[0].elts
                and isinstance(node.args[0].elts[0], ast.Constant)
                and isinstance(node.args[0].elts[0].value, str)):
            appended += 1
    return listed + appended


# Everything below one of these markers is a closed record of a past state. The guard's own
# doctrine is that this project does not rewrite its own log, and an explicit archive marker
# is a far less fragile signal of "history" than a transition arrow on the same line.
ARCHIVE_MARKERS = ("PROBLEM 1 \u2014 CLOSED", "ARCHIVE BELOW", "CLOSED, ARCHIVE")


def live_text(text):
    """The part of a document that still makes claims: everything above the archive marker."""
    cut = len(text)
    for marker in ARCHIVE_MARKERS:
        i = text.find(marker)
        if 0 <= i < cut:
            cut = i
    return text[:cut]


def limitation_order(text):
    """[] if section 14's numbered items run 1..N in order, else the offending jumps.

    WHY. Twice on 2026-10-08 a new limitation was inserted BEFORE the one it should follow,
    leaving the list reading 52, 54, 53. count_limitations takes the MAX, so the count stayed
    right and the gate stayed green while the document was visibly out of order. A reader who
    cites "limitation 54" would land on the wrong item.
    """
    i = text.find("## 14. Limitations")
    if i < 0:
        return ["section 14 not found"]
    rest = text[i + 1:]
    j = rest.find("\n## ")
    sec = rest[:j] if j >= 0 else rest
    nums = [int(m.group(1)) for m in re.finditer(r"^(\d{1,3})\.\s+\*\*", sec, re.M)]
    bad = [f"{a} then {b}" for a, b in zip(nums, nums[1:]) if b != a + 1]
    return bad


def live_claims(text, pattern=CLAIM):
    """[(lineno, claimed_int)] for live claims only; history (a transition arrow) is exempt."""
    out = []
    for n, line in enumerate(text.splitlines(), 1):
        for m in pattern.finditer(line):
            if ARROW.search(line[:m.start()]):
                continue                      # "43 -> 48 limitations": a record, not a claim
            for v in m.groups():
                if v is not None:
                    out.append((n, int(v)))
    return out


def _audit(label, truth, pattern):
    """Compare every live claim of one kind, on every surface, against the code's own count."""
    print(f"{label}: {truth}")
    bad = 0
    for f in SURFACES:
        rel = os.path.relpath(f, ROOT)
        if not os.path.exists(f):
            print(f"  SKIP   {rel} -- not present (clone without the vault)")
            continue
        claims = live_claims(live_text(open(f).read()), pattern)
        if not claims:
            print(f"  ok     {rel} -- states no count")
            continue
        for ln, v in claims:
            if v == truth:
                print(f"  ok     {rel}:{ln} says {v}")
            else:
                print(f"  FAIL   {rel}:{ln} says {v}, the code has {truth}")
                bad += 1
    return bad


def check():
    if not os.path.exists(METHODS):
        print(f"SKIP  {METHODS} absent"); return 0
    lim = count_limitations(open(METHODS).read())
    if lim is None:
        print("FAIL  could not count section 14's limitations"); return 1
    gates = count_gates()
    if gates is None:
        print("FAIL  could not count gate_sweep.py's GATES list"); return 1
    order = limitation_order(open(METHODS).read())
    if order:
        print(f"FAIL  section 14's limitations are out of order: {'; '.join(order)}")
        return 1
    print("section 14's limitations run 1..N in order")
    bad = _audit("section 14 has N numbered limitations, N =", lim, CLAIM)
    print()
    bad += _audit("gate_sweep.py runs N gates, N =", gates, GATE_CLAIM)
    print("\nPASS: every live surface agrees" if not bad else f"\nFAIL: {bad} stale count(s)")
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
    # --- the gate count, added 2026-10-08 after HANDOFF.md and ROADMAP.md both still said 22
    #     while the list had grown to 23. Same failure mode, a different number.
    src = ("GATES = [\n    ('a', [], None),\n    ('b', [], None),\n]\n"
           "def main():\n    results = []\n    results.append(('references', True, ''))\n"
           "    results.append((name, None, ''))\n")
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write(src); tmp = fh.name
    # 2 listed + 1 appended under a literal name; the `name`-variable append is NOT a gate
    assert count_gates(tmp) == 3, count_gates(tmp)
    os.unlink(tmp)
    # a file with no GATES assignment must refuse, not report zero
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
        fh.write("X = [1, 2, 3]\n"); tmp = fh.name
    assert count_gates(tmp) is None
    os.unlink(tmp)
    assert count_gates("/nonexistent/gate_sweep.py") is None
    # MUTATION 7: an out-of-order list must be caught even though the MAX is still correct.
    # This is the error made twice on 2026-10-08.
    good = "## 14. Limitations\n\n1. **a** x\n2. **b** y\n3. **c** z\n\n## 15. Next\n"
    assert limitation_order(good) == [], limitation_order(good)
    swapped = "## 14. Limitations\n\n1. **a** x\n3. **c** z\n2. **b** y\n\n## 15. Next\n"
    assert limitation_order(swapped), "a 1,3,2 list must be refused"
    assert count_limitations(swapped) == 3, "and the MAX alone would not have caught it"
    # the three accepted phrasings
    assert live_claims("all 23 gates\n", GATE_CLAIM) == [(1, 23)]
    assert live_claims("53 limitations, 23 green gates, and a fresh clone\n", GATE_CLAIM) == [(1, 23)]
    # MUTATION 3: "N of N" must check BOTH halves, or the box could say "23 of 22 gates"
    assert live_claims("the box says 23 of 23 gates\n", GATE_CLAIM) == [(1, 23), (1, 23)]
    assert live_claims("the box says 23 of 22 gates\n", GATE_CLAIM) == [(1, 23), (1, 22)]
    # MUTATION 4: a stale gate count must be CAUGHT
    assert live_claims("all 22 gates\n", GATE_CLAIM) == [(1, 22)]
    # history stays exempt here too
    assert live_claims("- gate count 22 -> 23 gates\n", GATE_CLAIM) == []
    # MUTATION 5: the archive is history. A past status report below the marker carries no
    # arrow, so only the marker saves it -- HANDOFF.md's problem-1 archive says "14 gates".
    archived = "live says 23 gates\n# PROBLEM 1 \u2014 CLOSED, ARCHIVE BELOW\nit said 14 gates\n"
    assert live_claims(live_text(archived), GATE_CLAIM) == [(1, 23)]
    assert live_claims(archived, GATE_CLAIM) == [(1, 23), (3, 14)]
    # and a limitation count must not be read as a gate count
    assert live_claims("53 limitations\n", GATE_CLAIM) == []
    # MUTATION 6: "gates" as a VERB is not a count. Dropping the lookahead makes this red.
    assert live_claims("item 10 gates every writing item\n", GATE_CLAIM) == []
    assert live_claims("item 10 gates the upload\n", GATE_CLAIM) == []
    # ...but the noun reading still has to be caught
    assert live_claims("10 gates, all green\n", GATE_CLAIM) == [(1, 10)]
    print("  ok  MUTATION: a 1,3,2 limitation list is refused though its MAX is still 3")
    print("check_published_counts.py --selftest PASS")
    return 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else check())
