#!/usr/bin/env python3
"""Every pre-submission gate, one command, non-zero exit if any fails.

The gates existed but were run by hand and from memory, which meant the set actually run
varied between sessions -- and the drift guard in the ipSAE bundle only fired at all
because someone happened to run that bundle. A sweep that is one command is a sweep that
gets run.

    gate_sweep.py            # everything
    gate_sweep.py --fast     # skip the structure-scoring gates (fixtures, fail-closed)

Gates, in order of what they protect:

  emit --selftest        the CSV emitter's own ranking fixtures
  gen_methods --check    the four GENERATED blocks still match the CSV
  check_claims           every numeric claim traces to its artifact (rules, citations,
                         inventory) -- bin/check_claims.py
  check_discards         no discarded design outranks a shipped one unmeasured
  instrument_v2 --self-test   seed grouping and max/median divergence
  control_family_balance      the control panel recomputes
  run_fixtures --check   12 ipSAE cases through all three production parsers
  test_failclosed        10 regressions on the two fail-open faults
  references             no § pointer in METHODS lacks a heading

The reference check is here because it catches a class the numeric gate cannot: §11.8 was
cited by §11.7 before it existed, and §5b pointed at a heading that lives in a different
document. Neither is a number, so check_claims would never have seen them.
"""
import os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VENV = os.path.join(ROOT, '.venv', 'bin', 'python')
PY = VENV if os.path.exists(VENV) else sys.executable

FAST_SKIP = {'run_fixtures --check', 'test_failclosed'}

GATES = [
    ('emit --selftest',          [PY, 'bin/emit_submission_csv.py', '--selftest'], None),
    ('gen_methods --check',      [PY, 'bin/gen_methods_submission.py', '--check'], None),
    ('check_claims',             [PY, 'bin/check_claims.py'], None),
    ('check_discards',           [PY, 'bin/check_discards.py'], None),
    ('instrument_v2 --self-test',[PY, 'bin/instrument_v2.py', '--self-test'], None),
    ('control_family_balance',   [PY, 'bin/control_family_balance.py'], None),
    ('run_fixtures --check',     [PY, 'run_fixtures.py', '--check'], 'outbox/ipsae-fixtures'),
    ('test_failclosed',          [PY, 'test_failclosed.py'], 'outbox/ipsae-fixtures'),
]


def dangling_refs():
    """§ pointers in METHODS with no matching heading. Cross-document refs are exempt
    when they name the other document explicitly."""
    p = os.path.join(ROOT, 'submissions/01-egfr-METHODS.md')
    t = open(p).read()
    have = {m.group(1) for m in re.finditer(r'^#{2,4}\s+(\d+[a-c]?(?:\.\d+)?)', t, re.M)}
    bad = []
    for m in re.finditer(r'§\s?(\d+[a-c]?(?:\.\d+)?)', t):
        pre = t[max(0, m.start() - 24):m.start()]
        if re.search(r'(CONTROL-TABLE|PREREGISTRATION|HANDOFF|README)\s*$', pre):
            continue
        if m.group(1) not in have:
            line = t[:m.start()].count('\n') + 1
            bad.append(f"METHODS:{line} points at §{m.group(1)}, no such heading")
    return bad


def main():
    fast = '--fast' in sys.argv
    results = []
    for name, cmd, cwd in GATES:
        if fast and name in FAST_SKIP:
            results.append((name, None, 'skipped (--fast)')); continue
        wd = os.path.join(ROOT, cwd) if cwd else ROOT
        if not os.path.exists(os.path.join(wd, cmd[1])):
            results.append((name, None, f'missing: {cmd[1]}')); continue
        r = subprocess.run(cmd, cwd=wd, capture_output=True, text=True)
        tail = [l for l in (r.stdout or '').strip().splitlines() if l.strip()]
        results.append((name, r.returncode == 0, tail[-1][:88] if tail else
                        (r.stderr or '').strip().splitlines()[-1][:88] if r.stderr else ''))

    bad = dangling_refs()
    results.append(('references', not bad,
                    'no dangling § pointers' if not bad else f'{len(bad)} dangling'))

    print(f"{'gate':<28}{'result':<8}note")
    print('-' * 96)
    failed = []
    for name, ok, note in results:
        tag = 'SKIP' if ok is None else ('PASS' if ok else 'FAIL')
        if ok is False:
            failed.append(name)
        print(f"{name:<28}{tag:<8}{note}")
    for b in bad:
        print(f"    {b}")
    print('-' * 96)
    if failed:
        print(f"FAIL: {len(failed)} gate(s) failed -> {', '.join(failed)}")
        sys.exit(1)
    print(f"PASS: all gates green"
          f"{' (fast mode: 2 skipped)' if fast else ''}")


if __name__ == '__main__':
    main()
