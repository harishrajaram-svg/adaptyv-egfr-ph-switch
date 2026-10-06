# VOID — LT-α constructs carrying an uncleaved signal peptide

These four files are **quarantined, not deleted**. Do not launch anything from them.

Every LT-α chain here is the pre-§27 **alignment-derived** construct, and it is wrong in two
independent ways:

| | broken (here) | correct (`targets/e2/`) |
|---|---|---|
| N-terminus | `ERFLPRTHLLLLGLLLVLLPGAQ…` — 37 res of signal peptide | `KPAAHLIGDPSK…` |
| mature domain | **120 res** | **144 res** (~24 missing internally) |
| provenance | sequence alignment | observed **1TNR chain A**, byte-identical |

The trimer never assembles. Every inter-chain ipSAE reads **exactly 0.0000**, including the
LT-α protomers' own mutual packing (`A:B`, `A:C`, `B:C`), where human TNF reads **0.71–0.74**.
pLDDT 0.27–0.39.

**Why quarantine instead of leaving them in place.** `bin/score-esmfold2.sh` accepts a
DIRECTORY and globs `*.faa` inside it. These files sat in `targets/validation-tnf/`, the
validation gate's own output directory, written at 1:32 PM — *after* the gate wrote its real
arms at 11:05 AM — so they looked like gate arms while nothing generated or consumed them.
Pointing the scorer at that directory would have folded them silently.

**Nothing current depends on them.** `runs/gate-dve/` is the void first attempt §27 already
supersedes; the live `validation_gate_tnf.py` builds its arms from 3ALQ and 3WD5 and never
touches LT-α; §27's 0.87 result came from `targets/e2/`, the corrected construct. The one live
consumer was `g-fab`'s G4 leg, now replaced by `targets/g4b/` (§30).

Correct construct: **`targets/e2/E2_LTa_TNFR1_POSITIVE.faa`** chains A/B/C, or rebuild from
`analysis/02-tnf/structures/1TNR.cif` chain A. Audit with
`python3 analysis/02-tnf/construct_audit.py`.

## 2026-10-06, later: the original G4 moved here too

`g-fab__G4_adalimumab_LTalpha.faa` was the one file with a LIVE consumer — it produced the
voided G4 leg of `g-fab` (§30). It is quarantined for the same reason as the other four:
`targets/g-fab/` is a directory, `score-esmfold2.sh` globs directories, and that is exactly
how the void run was launched. Re-running `score-esmfold2.sh targets/g-fab` now folds the four
G1 affinity-ladder complexes and nothing else, which is correct — the broken G4 should not be
refoldable from there. Its replacement is `targets/g4b/G4b_adalimumab_LTalpha_1tnr.faa`.
