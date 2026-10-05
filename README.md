# pH-switchable EGFR binder design

Anthropic × Adaptyv Bio protein design competition, 2026 — **Challenge 1, Track 3**.

Design a de novo binder to human EGFR that binds **more tightly at pH 6.5 than at pH 7.4**,
cross-reactive with mouse EGFR.

## Start here

| | |
|---|---|
| **[submissions/01-egfr-METHODS.md](submissions/01-egfr-METHODS.md)** | The methods document. Read §4.4 first. |
| **[submissions/01-egfr.csv](submissions/01-egfr.csv)** | The submission: 18 designs, ranked on the two-partner histidine-only pH product. |

## What this submission claims, in four lines

1. **One mutation makes the switch, and it replicated four times.** A single Ser→Asp on the
   best-binding BindCraft backbone took four independent ProteinMPNN sequences from 3.15–3.85×
   to 5.40–5.46× at no cost in predicted affinity. The matched wild-type is submitted alongside
   so the comparison gets made in the laboratory, not in our gate.
2. **A molecule with NO REPORTED KD reads a 5.27× switch.** Of 11 molecules Adaptyv ran on this
   platform, the highest-scoring one on our own ranking metric is a molecule with **no KD
   reported** — which is right-censored, *not* a measurement that it does not bind (PK,
   2026-10-04; METHODS §4.1). Its rank is **8th of 132 rankable molecules**; an earlier version
   of this file quoted "8th of 2,009", which METHODS §7 records as a ~15× overstatement because
   1,877 of that denominator were never rankable.
3. **Requiring both species is what catches it.** Both molecules that outrank the measured
   binder on human score exactly 0.0000 on mouse — no interface at all, not a narrow miss.
   Cross-reactivity was in the brief; it turns out to be the only specificity filter here that
   measured data supports. The panel is **10 right-censored molecules and one positive**, so
   this is a statement about ranks and a mechanism, not a validated error rate — and with
   n = 1 positive it cannot be an error rate at all.
4. **The single-site ceiling is 5.55× and we failed to build the route past it — which is not
   the same as the route being closed.** H433's free pKa is 6.22, so no single-site design can
   beat 5.55× over a 0.9 pH-unit window. We reached the pair that bridges it (H433+H370, 8.5 Å)
   only 3 times in 1,944 designs, all incidental. METHODS §8.1 and §3.4 both retract the
   stronger claim by name: **the two-site route is demonstrated, not closed** — G532 is a
   published molecule that does it. This section has been wrong twice in that direction.

We submitted **18 designs of the 20 allowed**. The ones we left out did not stand on a measurement —
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

**Second block, same day — the review's remaining items.**

- **The exclusion ledger was rebuilt and it did not favour the submission.** 45 warn names
  reconcile to 15 distinct sequences (not 38); that list and the separate "never a candidate"
  list are **disjoint**, so the union is 75 molecules; 63 were excluded on the pH gate alone and
  were reopened. Re-scored on the finalists' own gate, 35 of our own designs outranked the
  weakest shipped tier-1 design, and novelty cleared 25 of them at Level 3. **Five were added to
  the submission (12 → 17 of the 20 permitted)** by a rule fixed before the result was examined.
- **The pre-submission gate was checking a 31-design candidate pool, not the 12 shipped**, so
  its threshold came from a design we did not submit and 19 unshipped candidates were never
  checked at all.
- **Pose discovery was undercounting and double-counting simultaneously.** Three submitted
  designs had four times more poses on disk than the submission counted; 76 sequences had an
  inflated pose count because the same files were globbed twice. Medians barely moved; spreads
  moved a lot (1.31 → 4.38 on one design). No gating decision was wrong.
- **On its own co-crystal, the instrument scores a correct interface zero.** rAC1 against 4UIP
  by structural contact recovery: one pose in ten reproduces 72% of the crystal contacts and 93%
  of the epitope, and **all ten poses score ipSAE_min 0.0000**. A zero on this instrument does
  not mean "no interface".
- **The submission has nine backbone families and one epitope.** All 17 designs contact the same
  patch of domain III, 20 residues shared by ≥80% of them. The top-ranked design contacts an
  N-glycosylation sequon (Asn420). Human/mouse identity *at the contacted positions* is 0.86,
  not the whole-protein figure.
- **Control recovery is reported family-balanced**, because nine of the ten no-KD molecules come
  from one submitter group. Raw 8/10 = 0.800; family-balanced 0.889; leave-one-family-out spans
  0.778–1.000. No interval is reported and no effective-n is substituted into Clopper–Pearson.
- **The 0.2218 `affinity_above_null` column is removed from the CSV**, not relabelled — it was a
  percentile of a null that proved to be a point mass at zero.
- **Two corrections to this day's own work**, both caught before release: 6ARU is the
  cetuximab-Fab complex in the *tethered* conformation, which matches the assay, so an earlier
  claim here that receptor state was unmatched and unassessed was wrong; and the count of
  designs carrying binder histidines is **nine of eighteen** — a figure that was wrong five times
  by hand and is now a generated block (METHODS §11.1), never typed.
- **Near-duplicate pairs are declared, not smoothed over.** Three pairs exceed 90% identity
  (`rimA01_r15_d3_rimA_20`/`L133E` at 0.993, `ss_bc_s831683_mpnn6_S15D_S62H_routeA`/
  `bc_s831683_mpnn6_S15D` at 0.985, `bc_s831683_mpnn9_S15D`/`mpnn9_WT` at 0.985), so **6 of 18**
  designs sit in one. A fourth pair sits at 0.867
  (`bcr_d3acid3_l60_s647537_mpnn3`/`mpnn11`) — both are MPNN redesigns of one backbone, and both
  CSV rows previously claimed "<0.47 identity to anything else submitted". Corrected 2026-10-05.

**Since closed.** All four items listed here as open on 2026-10-04 were completed on 10-05 and
this paragraph described them as open for a day afterwards:

- **The free leg now has three rungs, not one.** Partner-deletion (shipped), a backbone-restrained
  relaxed leg and a separately-folded apo leg, with the bound leg identical in all three
  (METHODS §11.7, §11.8). The relaxed leg moves every testable design *down* — median 0.842×,
  so the shipped basis is optimistic — while leaving the order intact (Kendall τ = +0.956).
- **Full-ECD, glycan and receptor-state checks were applied to the finalist footprints.**
  `analysis/01-egfr/finalist_footprints.json` carries `n_outside_d3_crop`, `glycan_sequon_hits`
  and `hu_mo_identity_at_epitope` for all 18 designs (METHODS §10b). Three designs contact the
  Asn420 sequon, not one.
- **The rAC1 comparison runs on structural contact recovery**, 5.0 Å heavy-atom contact sets
  with numbering mapped by alignment (`bin/rac1_contact_recovery.py`), on both predictors
  (METHODS §4.4b) — not predicted confidence attached to crystallographic coordinates.

What remains genuinely open is in METHODS §13, now 24 numbered limitations, and in the one red
gate: **novelty eligibility is unresolved for four shipped designs** (`bin/check_novelty_coverage.py`).

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
