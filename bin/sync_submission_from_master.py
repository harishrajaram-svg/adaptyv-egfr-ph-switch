#!/usr/bin/env python3
"""Re-point every submission row's pH ratio at the sequence-pooled value.

`submission_final.json` carried ratios written by several different scripts over three
days, keyed on run names. `analysis/01-egfr/master_rank.json` is keyed on the binder
SEQUENCE and pools every human-leg refold pose of that exact molecule, so it is the
only column where one number means one thing. This script copies ratio / ratio_n /
site across by sequence and PRINTS EVERY CHANGE -- a silent re-point would be worse
than the inconsistency it fixes.

Rows whose sequence is absent from master_rank.json are left untouched and listed; that
means no ESMFold2 human pose of that exact sequence exists, which is itself a finding.

    python3 bin/sync_submission_from_master.py [--write]
"""
import json, sys

SUB = 'analysis/01-egfr/submission_final.json'
MASTER = 'analysis/01-egfr/master_rank.json'

rows = json.load(open(SUB))
master = {r['seq']: r for r in json.load(open(MASTER))}

changed, missing = [], []
for x in rows:
    m = master.get(x['seq'])
    if not m:
        missing.append(x['name']); continue
    old, oldn = float(x.get('ratio', 0)), int(x.get('ratio_n', 0))
    new, newn = m['ratio'], m['ratio_n']
    if new is None: missing.append(x['name']); continue
    if abs(new - old) > 1e-3 or newn != oldn:
        changed.append((x['name'], old, oldn, new, newn))
    x['ratio'], x['ratio_n'], x['ratio_site'] = new, newn, m['site']

print(f"{len(rows)} rows; {len(changed)} re-pointed, {len(missing)} not in master\n")
print(f"{'was':>7}{'n':>4}  ->{'now':>7}{'n':>4}   design")
for n, o, on, w, wn in sorted(changed, key=lambda t: -abs(t[3] - t[1])):
    print(f"{o:7.3f}{on:4d}  ->{w:7.3f}{wn:4d}   {n[:52]}")
if missing:
    print(f"\nno pooled pH for (sequence never refolded against human):")
    for n in missing: print(f"   {n}")

if '--write' in sys.argv:
    json.dump(rows, open(SUB, 'w'), indent=1)
    print(f"\nwrote {SUB}")
else:
    print("\n(dry run -- pass --write to apply)")
