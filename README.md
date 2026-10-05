# pH-switchable EGFR binder design

Anthropic × Adaptyv Bio protein design competition, 2026 — **Challenge 1, Track 3**.

Design a de novo binder to human EGFR that binds **more tightly at pH 6.5 than at pH 7.4**,
cross-reactive with mouse EGFR.

## Start here

| | |
|---|---|
| **[submissions/01-egfr-METHODS.md](submissions/01-egfr-METHODS.md)** | The methods document. Read §4.4 first. |
| **[submissions/01-egfr.csv](submissions/01-egfr.csv)** | The submission: 10 designs, ranked. |

## What this submission claims, in four lines

1. **One mutation makes the switch, and it replicated four times.** A single Ser→Asp on the
   best-binding BindCraft backbone took four independent ProteinMPNN sequences from 3.15–3.85×
   to 5.40–5.46× at no cost in predicted affinity. The matched wild-type is submitted alongside
   so the comparison gets made in the laboratory, not in our gate.
2. **A molecule with no reported KD reads a 5.27× switch.** Of 11 molecules Adaptyv measured on this
   platform, the highest-scoring one on our own ranking metric is a design already measured
   **not to bind** — and it ranks 8th of 2,009 on the competition's primary objective.
3. **Requiring both species is what catches it.** Both molecules that outrank the measured
   binder on human score exactly 0.0000 on mouse — no interface at all, not a narrow miss.
   Cross-reactivity was in the brief; it turns out to be the only specificity filter here that
   measured data supports. The panel is 10 negatives and one positive, so this is a statement
   about ranks and a mechanism, not a validated error rate.
4. **The single-site ceiling is 5.55× and the route past it is closed.** H433's free pKa is
   6.22, so no single-site design can beat 5.55× over a 0.9 pH-unit window. The only histidine
   pair close enough to bridge (H433+H370, 8.5 Å) is unreachable: 3 of 1,944 designs hit both,
   all by accident, and 9 of 10 that switched there did not bind.

We submitted **10 designs of the 20 allowed**. The other ten did not stand on a measurement —
eight read *below* 1.0× and sat on the 0.702× value that a large share of designs return.
(That value is **not** a "no-switch floor" — no linkage reads 1.0. 0.702× is the analytic
*acid-weakening* extreme, and a pile-up there indicates protonation-model saturation, so those
eight are uninterpretable on this gate rather than measured non-switchers. Corrected
2026-10-05.) [§11](submissions/01-egfr-METHODS.md) explains the cut.

## Reproducing it

`bin/` holds every analysis step as a standalone script; several carry `--selftest`.

```
bin/ph_gate_refolds.py        pH gate on independent ESMFold2 refolds
bin/ph_pool_by_sequence.py    one ratio per MOLECULE, pooled by binder sequence
bin/master_rank.py            pH + affinity joined BY SEQUENCE -- the single ranking table
bin/instrument_v2.py          ipSAE_min, median over 5 seeds   (--self-test)
bin/novelty_gate.py           Adaptyv's general-protein novelty levels
bin/antibody_novelty.py       their antibody branch, CDRH3     (--selftest)
bin/express_qc.py             cell-free expression liabilities
bin/emit_submission_csv.py    builds the CSV                   (--selftest)
bin/check_discards.py         refuses to ship while an unmeasured design outranks a
                              submitted one                    (--selftest)
```

Joins are on the **binder sequence**, never the run name: the same molecule appears in this
project under as many as three names, and name-keyed joins broke five separate analyses before
`master_rank.py` made it impossible.

## The error history is the point

§7 of the methods document lists thirteen errors we found in our own instrument, each one a
commit in this repository rather than a quiet rewrite. Several reversed a published-in-log
conclusion. Two are worth knowing before you read any number here:

- Our negative control turned out to be **81% mature human EGF** — the agonist. Every bar in
  the project had been calibrated against it. The thresholds are retired, not replaced.
- The novelty rule was implemented with one clause missing, which rejected **114 designs that
  actually clear**, including the entire pool that supplies ranks 1–4, 8 and 9.

Environment setup from the start of the project is in
[docs/setup-notes-2026-09-18.md](docs/setup-notes-2026-09-18.md) and is largely out of date.
