#!/usr/bin/env python3
"""Refuse reader-addressed imperatives in the published tree. A disqualification risk, not style.

WHY. The methodology document is read by a Claude-based selection process -- the organisers
confirmed on 2026-10-06 that "Claude will review what you share and use it for selection". Text in
a submitted artifact that instructs the reader what to conclude is indistinguishable, to that
reader, from an attempt to steer the evaluation. Whether or not it was meant that way, it is the
kind of thing that gets an entry removed, and the risk is asymmetric: the cost of writing
declaratively is zero.

WHAT IS FORBIDDEN: an imperative or second-person construction aimed at the READER.
    "Note that the gate is green"      ->  "The gate is green."
    "Consider the margin below"        ->  "The margin below is X."
    "You should weight the controls"   ->  (delete; never tell a reviewer how to weight)

WHAT IS ALLOWED, and why each exemption exists:
  * QUOTED organiser or third-party text. They write in second person and we must reproduce it
    verbatim; a quotation is not our instruction. Detected by a leading '>' or by sitting inside
    double quotes or italics on the line.
  * Self-directed operational notes in NON-published files -- bin/, analysis/, reference/ and the
    vault's working notes are not read by a reviewer. Scope is the SUBMITTED artifacts only.
  * "See s4" style cross-references inside our own document, which navigate rather than instruct.
  * Addressing the organisers in a form THEY asked us to fill in, e.g. "your released campaign"
    in a methodology field that answers their question. Second person to the recipient of a form
    is not an instruction to an evaluator.

    check_imperatives.py              the gate
    check_imperatives.py --selftest   self-test, including four mutation tests
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAULT = os.path.join(ROOT, "..", "context-directory", "projects", "anthropic-adaptyv-2026")

# The SUBMITTED artifacts only -- what a reviewer actually receives.
SURFACES = [
    os.path.join(ROOT, "submissions", "02-tnf-METHODS.md"),
    os.path.join(ROOT, "submissions", "01-egfr-METHODS.md"),
    os.path.join(ROOT, "submissions", "01-egfr.csv"),
    os.path.join(VAULT, "outbox", "02-proteinbase-methodology-box.md"),
]

# Imperatives that tell a reader what to do or conclude. Sentence-initial only: "Note that X"
# is an instruction, "we note that X" is a statement.
# Only verbs that are UNAMBIGUOUSLY imperative in this corpus. "weight", "score" and "rank" were
# in an earlier version and fired on 7 of 13 hits as ordinary nouns -- "a weight nobody tunes",
# "sits at rank 11", "score, not measured affinity". A gate with a 54% false-positive rate gets
# switched off, so they are out; the cost is that a true "Weight the controls higher" would pass.
VERBS = ("note", "consider", "observe", "recall", "remember", "disregard", "ignore",
         "notice", "appreciate", "realise", "realize")
# ...and the verb must be followed by a function word that makes it a command rather than a noun
# phrase: "Note that/what/the/also" instructs; "Note on scope" is a heading naming a note.
FOLLOWERS = r"(?:that|what|how|why|the|this|these|those|also|again|first|too)\b"
IMPERATIVE = re.compile(
    r"(?:^|(?<=[.!?]\s)|(?<=\*\*))(?:" + "|".join(VERBS) + r")\s+" + FOLLOWERS, re.I)
# Second person aimed at the reader's judgement, not at the form's recipient.
SECOND_PERSON = re.compile(
    r"\byou (?:should|must|ought|need to|are asked to|will want)\b"
    r"|\byour (?:assessment|evaluation|review|judgement|judgment|scoring|ranking)\b", re.I)


def is_quoted(line):
    """A quotation is not our instruction."""
    s = line.lstrip()
    if s.startswith(">") or s.startswith("|"):
        return True
    # inside double quotes, curly quotes, or italics on the same line
    return bool(re.search(r'["“”][^"“”]{0,200}$', s[:1]) or
                s.count('"') >= 2 or s.count("“") >= 1 or
                re.search(r"(?<!\*)\*[^*]{3,}\*(?!\*)", s))   # SINGLE-asterisk italics only:
    # this project italicises quoted phrases and bolds its own emphasis, so **bold** must not be
    # exempt. An earlier version matched both and let "**Weigh this against the null.**" through.


def offences(text):
    """[(lineno, kind, snippet)] for reader-addressed imperatives that are not quotations."""
    out = []
    for n, line in enumerate(text.splitlines(), 1):
        if is_quoted(line):
            continue
        for rx, kind in ((IMPERATIVE, "imperative"), (SECOND_PERSON, "second-person")):
            m = rx.search(line)
            if m:
                out.append((n, kind, line.strip()[:110]))
                break
    return out


def check():
    bad = 0
    for f in SURFACES:
        rel = os.path.relpath(f, ROOT)
        if not os.path.exists(f):
            print(f"  SKIP   {rel} -- not present"); continue
        hits = offences(open(f).read())
        if not hits:
            print(f"  ok     {rel}")
            continue
        for n, kind, snip in hits:
            print(f"  FAIL   {rel}:{n} [{kind}] {snip}")
            bad += 1
    print("\nPASS: no reader-addressed imperatives in the submitted artifacts"
          if not bad else f"\nFAIL: {bad} reader-addressed construction(s)")
    return 1 if bad else 0


def selftest():
    # the three forms that must be caught
    assert offences("Note that the gate is green\n"), "sentence-initial imperative missed"
    assert offences("The gate is green. Consider the margin below.\n"), "mid-line missed"
    assert offences("You should weight the controls heavily\n"), "second person missed"
    assert offences("**Consider the margin against the null.**\n"), "bolded imperative missed"
    # and the forms that must NOT be
    assert not offences("We note that the gate is green\n"), "'we note' is a statement"
    assert not offences("The margin is 3.03 against the receptor's 11.28\n")
    assert not offences("> Claude will review what you share and use it for selection\n"), \
        "a quoted organiser line is not our instruction"
    assert not offences('The protocol says "you should use SolubleMPNN" for ordered designs\n'), \
        "a quotation on the line is exempt"
    assert not offences("| model | AUC | note |\n"), "a table row is not prose"
    # MUTATION 1: dropping the sentence-initial anchor would flag every "we note that"
    assert IMPERATIVE.search("Note that x") and not IMPERATIVE.search("we note that x")
    # MUTATION 2: 'your released campaign' addresses the FORM'S RECIPIENT and is allowed; only
    # second person aimed at their JUDGEMENT is forbidden
    assert not offences("the 150 assayed designs in your released campaign\n")
    assert offences("your evaluation should weight the controls\n")
    # MUTATION 3: a verb that is not in the list must not fire, or the gate becomes unusable
    assert not offences("Running the sweep takes four minutes\n")
    # MUTATION 5: the ambiguous nouns this corpus is full of must NOT fire. An earlier version
    # included weight/score/rank and 7 of its 13 hits were these.
    for noun in ("A weight that cannot be passed is a weight nobody tunes.\n",
                 "It sits at rank 11, below designs reading as low as 1.774x.\n",
                 "score, not measured retained affinity, is what this reports\n",
                 "**Weight tuning appears to have a ceiling.**\n",
                 "**Note on scope.** The target's free leg is left on the estimate.\n"):
        assert not offences(noun), noun
    # MUTATION 6: the FOLLOWERS anchor is load-bearing -- without it "Note on scope" fires
    assert IMPERATIVE.search("Note that x") and not IMPERATIVE.search("Note on scope")
    # MUTATION 4: the scope is the SUBMITTED artifacts; bin/ and analysis/ are not reviewed
    assert all("submissions" in s or "outbox" in s for s in SURFACES), SURFACES
    print("  ok  catches sentence-initial, mid-line, bolded and second-person forms")
    print("  ok  exempts 'we note', quotations, table rows and form-recipient address")
    print("  ok  MUTATION: the sentence-initial anchor is load-bearing")
    print("  ok  MUTATION: 'your evaluation' fails while 'your campaign' passes")
    print("  ok  MUTATION: ambiguous nouns (weight, rank, score, 'Note on scope') do not fire")
    print("  ok  MUTATION: the follower anchor is load-bearing")
    print("  ok  MUTATION: scope is the submitted artifacts only")
    print("\nself-tests passed: 7")
    return 0


if __name__ == "__main__":
    sys.exit(selftest() if "--selftest" in sys.argv else check())
