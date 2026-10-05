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

## Corrections — 2026-10-05

An external reviewer went through the public tree, the CSV and the code on 2026-10-05. What
follows is the full list of what changed as a result, dated, so that anyone reading an earlier
copy of this repository can tell what it got wrong. **This review was a critical read, not a
certification**: the reviewer found specific errors and asked for specific repairs, and nothing
below should be read as an endorsement of the submission or of the claims that survive.

**Scoring code — three fail-open paths, none of which changed a submitted number.**
`bin/ipsae_min.py` could return a score from a stale output file after the scorer crashed (no
return-code check, no invalidation), and a file holding only one of an interface's two
alignment directions was accepted as a minimum over two. The same single-direction fault, plus
a complete absence of chain-pair grouping, was also in `bin/master_rank.py` — the parser that
writes the file this submission is built from — so on a complex with more than two chains
`min()` could report an intra-binder interface as the binder's binding score. All three now
fail closed. Verified: all 6,840 cached scoring outputs are single-pair and two-direction, so
none of the three paths was ever taken; the canonical file regenerates with zero field changes
across 2,049 sequences and the CSV is byte-identical. The faults were latent, not harmless —
a scorer that can invent a number cannot be a reproducible path — and `outbox/ipsae-fixtures/`
now scores every case through all three production parsers, fails on any disagreement, and
carries ten regressions that were checked to fail against the pre-fix code.

**The vendored reference was missing from the public tree.** `.gitignore` carried `ipsae/`,
which matches a directory of that name at *any* depth, silently excluding
`outbox/ipsae-fixtures/vendor/ipsae/`. Anyone who cloned this repo could not run the fixtures.
Scoped to `/ipsae/`.

**The pH ranking basis was not what this repository said it was.** The gate was documented here
and in the preregistration as composing "every titratable site on both partners". It computed
pKa values for histidines, aspartates and glutamates and then composed **only the
histidines** — and the designed intervention in most submitted families *is* an acid (A22D,
S88D, S15D, L133E, T65D), so the gate was blind to the residue each design was built around.
The shipped column is renamed to say so. Three bases are now reported side by side
(histidine-only, all-site, and sites with a counter-charge within 6 Å), computed from the same
poses through one code path. **Kendall τ between the shipped basis and the partnered basis is
+0.000** — the orderings are uncorrelated, designs move by up to 8 ranks — so **every tier in
the CSV is now marked `provisional`**. The shipped order is retained because the
histidine-only value is the minimum of the three for all twelve designs, making it the
conservative envelope under one uniform rule; neither wider basis can rank, because the matched
negative control reads 5.3× on the partnered basis. See METHODS §11.7.

**"No linkage" is a ratio of 1.0, not 0.699×.** This repository described 0.702× as a
"no-switch floor" that any design touching a histidine returns. That was backwards: 0.699× is
the opposite extreme, the limit of maximal *acid-weakening* linkage, and designs piling up
within 0.003 of an analytic bound indicate protonation-model saturation. Those designs are
uninterpretable on the gate rather than measured flat. Also withdrawn: the argument that a
one-sided pKa shift must be real "because noise would scatter both ways" — a consistent
one-sided shift is what systematic model bias produces.

**The measured controls are right-censored, not non-binders.** The ten Adaptyv control
molecules were called "measured non-binders" throughout. They are molecules for which **no KD
was reported**: affinity known only to be weaker than the assay's quantifiable limit, a class
that also contains expression and QC failures. A missing KD bounds affinity from above, not
below. Every claim resting on them is restated as separating one quantified binder from ten
censored observations, which is weaker. Relatedly, the claim that G532 is "the only one that
binds" is withdrawn — G532Ctrl is also a measured EGFR binder.

**The pre-submission gate was checking the wrong submission.** `bin/check_discards.py`
compared against a 31-design candidate pool rather than the 12 designs actually shipped, so its
threshold came from a design we did not submit, and the 19 unshipped candidates were skipped by
a name-stub test and never checked against the shipped set at all. It now reads the graded CSV
and matches by sequence.

**The exclusion ledger is rebuilt and sequence-keyed.** METHODS §10 held two unreconciled
lists — 45 warn names said to be "38 distinct molecules", and a separate set of high-ranking
molecules that were never candidates, called 28 and then 59. `bin/exclusion_ledger.py` resolves
the first to **15 distinct sequences** (8 alias collapses), establishes that the two lists are
**disjoint**, and classifies the resulting 75-molecule union into eligibility failures,
inadequate assessment, and pH-gate-only exclusions — the last of which are reopened and
re-scored on the same gate as the finalists, because a gate-dependent rejection is not a
finding.

**Pose discovery was both undercounting and double-counting.** Every sequence-to-pose join went
through a `.faa` index that finds a pose directory only when a matching file happens to exist,
and that hits the same files twice when a sequence appears in two index directories. Measured
on the rank-1 design: 24 hits over 12 distinct files, every pose counted twice; 76 sequences
had an inflated pose count. Separately, three submitted designs had more poses on disk than the
submission was counting (n=5 → 20, 5 → 20, 11 → 26). `bin/pose_seq_index.py` now indexes poses
by the binder sequence read out of the structure. The medians barely moved (5.659 → 5.656);
the **spreads** moved a lot (1.31 → 4.38), because n=5 could not see them. No tier-1 gating
decision changed.

**Three claims about the L133E design did not survive deeper sampling** and are withdrawn by
name in METHODS §11.6: its human-leg affinity did *not* improve over its parent (0.598 vs
0.594 at n=20 — the apparent gain was five-pose noise); L133D's ratio is 2.221× not 0.723×; and
the glutamate-reach argument was geometric inference, now measured per pose. What strengthened
is the mechanism's specificity: a neutral isostere at the same position leaves the gate
unchanged and an acid one methylene shorter cannot reach its counter-charge and does nothing.

**Still open.** The partner-deletion free leg is a fixed-conformation diagnostic, not a
measurement of the apo state. Full-ECD, glycan and receptor-state checks have not been applied
to the finalist footprints. The rAC1 comparison needs structural contact recovery against
4UIP rather than predicted confidence attached to crystallographic coordinates. These are
recorded in METHODS §13 rather than resolved.

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
