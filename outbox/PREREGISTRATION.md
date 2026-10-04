# Preregistration — Adaptyv challenge 1 (EGFR pH-switch), December analysis

Written **2026-10-04**, before any experimental outcome is known, in answer to PK:
*"Freeze the candidate universe, sequence-family assignments, model/code versions, scoring
settings, selection history, final ranks and exclusions."*

**Status of the work it covers: EXPLORATORY.** Everything below the "frozen" section was
produced by repeated, outcome-driven screening of one pool. This document cannot undo that.
It fixes what happens *next*, and it documents the selection history so the bias is visible
rather than hidden.

---

# PART 1 — FROZEN STATE

## 1.1 Candidate universe

| | |
|---|---|
| designs generated | **1,948** across 46 arms (BoltzGen 1,877; BindCraft 71) |
| pH-gated on generator poses | 1,944 |
| re-gated on independent ESMFold2 refolds | 1,659 |
| scored by Proton-PottsMPNN | 1,693 |
| distinct binder sequences scored on the instrument | **2,009** (of which 27 are controls) |
| complexes folded and ipSAE-cached | **6,011** |
| **submitted** | **10** |

The universe is the 1,948. Designs outside it (e.g. BindCraft trajectories never sequenced)
are **not** eligible and are not nonbinders; they were never candidates.

**Amended 2026-10-04, 18:00 EDT, and the amendment is itself part of the selection history.**
The first draft of this document (13:18 EDT) froze a **20-design** submission. At 15:32 the
submission was cut to **10** — ranks 11–20 of that build did not stand on a measurement: eight
read *below* 1.0×, i.e. on the 0.702× steric floor that any design touching a histidine with no
carboxylate nearby returns. That cut is a post-freeze selection change and is recorded here
rather than silently absorbed, on your own standard: *"A correction does not invalidate
December's analysis if every candidate can be rescored consistently before outcomes are
examined; silent selective changes would."* Every candidate was rescored consistently; no
outcome has been examined.

Two counts in the first draft were also stale and are corrected above: 1,948 designs *generated*
is not the same as 2,009 distinct binder sequences *scored* (the gap is controls and
resequencings), and the cache now holds 6,011 scored outputs, not 5,368.

## 1.2 Versions, pinned

| component | version |
|---|---|
| this repo | `a5bd3d6`, 2026-10-04 — **tree clean, nothing uncommitted.** Public at <https://github.com/harishrajaram-svg/adaptyv-egfr-ph-switch> (the first draft pinned `76c5b0c` and warned of an uncommitted working tree; that no longer applies) |
| ipSAE reference | `ipsae/ipsae.py` v4, Dunbrack, commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3` |
| structure model | ESMFold2 via Modal, loops=10, steps=68 |
| pKa | PROPKA **3.5.1** |
| protonation model | Proton-PottsMPNN `potts_v6_afdb_edge_his0.3_acid0.06/epoch-0125.ckpt`, 21,494,143 bytes |
| H-bonds | HBPLUS v3.06 (McDonald), compiled locally |
| cutoffs | ipSAE PAE 10 / dist 10; pH 6.5 vs 7.4; RATIO_BAR 1.20; MIN_N 5 poses |

## 1.3 Sequence families

**5 families among the 10 submitted designs.** The first draft of this section described a
submission that was never shipped — it named `d2d_101_l120_s869126_mpnn1`/`_mpnn5` as the one
multi-member family, and neither design is in the file. The real partition:

| family | n | members (submission rank) | note |
|---|---|---|---|
| **`d3acid_l65_s831683`** | **6** | 1, 2, 3, 4, 8, 9 | all 65 aa on one BindCraft backbone: four S15D point mutants, one unmutated MPNN sequence, and the matched wild-type of rank 3 |
| `rimA02_d3_rimA_14` | 1 | 5 | VHH, 129 aa |
| `rimA01_r15_d3_rimA_20` | 1 | 6 | 150 aa BoltzGen |
| `d2c_101_l147_s144898` | 1 | 7 | 147 aa, the S88D causal result |
| `h370_020` | 1 | 10 | VHH, 98 aa |

**Consequence, and it is the whole point of your ask: the effective n of this submission is 5
clusters, not 10 designs — and 6 of the 10 sit in one cluster.** Any December confidence
interval, hit rate or test must be computed on families, not on designs. Ranks 1, 2, 3 and 4
differ from one another by nothing but the ProteinMPNN sequence on an identical backbone, and
ranks 3 and 9 differ by **one residue** (Ser15 vs Asp15). Treating those six as six independent
tests would inflate n by a factor of six on the arm carrying our only causal claim.

This concentration is deliberate and it is a trade we made with eyes open: diversity given up to
buy wet-lab replication of the one causal result in the project. It is also the thing you warned
against twice, so it should be scored as a known cost rather than discovered in December.

Control families (binder chain): EGF-derived (2 runs, 1 molecule), cetuximab scFv (4 runs,
1 molecule), VHH (2 runs, 1 molecule). **Three molecules, eight runs.** Plus, added 2026-10-04:
10 experimentally measured non-binders + human EGF (one group supplies 8 of the 10 — see
CONTROL-TABLE §6, where that family structure is an open problem, not a solved one), 4
framework-preserving CDR decoys, and 12 nested shuffle nulls.

Control families (binder chain): EGF-derived (2 runs, 1 molecule), cetuximab scFv (4 runs,
1 molecule), VHH (2 runs, 1 molecule). **Three molecules, eight runs.**

## 1.4 Selection history — every narrowing applied to the pool, in order

1. generator-native `iptm` (BoltzGen's own filter, top-30 of each arm retained)
2. H433 carboxylate contact < 4.0 Å
3. the 2.4–3.2 Å "predictive band"
4. pH ratio ≥ 3.0, then ≥ 1.5, then ≥ 1.20
5. ipSAE_min ≥ 0.1493 (ECD) / 0.2218 (d3) — **both now retired**
6. n ≥ 5 poses for a ratio to rank
7. novelty, as the organisers actually define it — **not** the "sequence identity < 30%" rule
   the first draft of this document froze. That rule omits the clause whose absence produced our
   bogus "0 of 120 pass novelty" result. The implemented rule (`bin/novelty_gate.py`,
   `bin/antibody_novelty.py`, both self-tested):
   - **general proteins:** TM ≥ 0.80 **alone** is Level 2 regardless of sequence identity;
     TM ≥ 0.50 is "moderate"; identity bars at 70% and 30%. Level 1 = identity >70% AND
     moderate structure. Level 2 = identity >70%, **or** high structure alone, **or**
     (identity >30% AND moderate). Level 3 = exactly one of (identity >30%, moderate structure).
     Level 4 = neither. The upload gate is **Level ≥ 3**.
   - **antibodies** (VHH/scFv/Fab) are levelled on a separate CDRH3 branch, which is far more
     permissive, and classification depends on the organisers' own ANARCI call.
   - our search is **FoldSeek against PDB only**; the organisers additionally search SwissProt,
     patents, a therapeutic-antibody database and PLAbDab, so our identities are **lower bounds**

Steps 1–5 were applied, revised and re-applied to the **same pool**, with thresholds chosen
after seeing the data. Three designs were promoted and later withdrawn on re-measurement.
Surviving designs are therefore selected on the same measurements December would analyse.

## 1.5 Exclusions, with reasons

**Corrected 2026-10-04.** The first draft of this table listed `rimA02/d3_rimA_14` and
`h370_only/h370_020` as excluded. **Both are submitted, at ranks 5 and 10.** They were excluded
under the general-protein rule and then readmitted on the antibody branch once
`bin/antibody_novelty.py` existed; the exclusion table was never updated. A preregistration that
lists shipped designs as excluded is worse than none, so both the error and the correction are
recorded.

**Still excluded**

| design | reason |
|---|---|
| `h370_2site/2site_009` | 74.5% identity to `4tqe_L`. Antibody format flag false and CDRH3 undelimitable, so the antibody branch does not apply; top hits are antibody *light* chains. Not recoverable. |

**Readmitted on the antibody branch — in the submission**

| design | general-protein rule | antibody rule | shipped as |
|---|---|---|---|
| `rimA02_d3_rimA_14_vhh` | TM 0.8552 ⇒ **Level 2**, would be auto-rejected | CDRH3 **13.6%**, global **84.8%** ⇒ **Level 3, clears** | rank 5, `nanobody` |
| `h370_020_vhh` | TM 0.914 ⇒ **Level 2** | CDRH3 **27.3%**, global **52.6%** ⇒ **Level 4, clears** — the only Level 4 in the project | rank 10, `nanobody` |

**Both therefore depend on the organisers classifying them as antibodies.** If their ANARCI does
not, both rows fail at upload. This is stated as a property of the submission, not a prediction
about it. (The global identity for `rimA02` was reported as 70.5% in an earlier build of the
submission CSV; 70.5% exists nowhere in the pipeline except a selftest fixture, and the artifact
value is 84.8%. Fixed before upload.)

**VHH format, affinity not assessed — in the submission**

| | |
|---|---|
| 2 of the 10 shipped designs (ranks 5, 10) | **affinity inadequately assessed**, not rejected on score. The first draft said "3 VHH-format designs" were set aside; two shipped. §6 of CONTROL-TABLE now shows the blind spot is total: four framework-preserving CDR decoys and the real 5 nM VHH all read 0.0000. |

## 1.6 Final ranks, frozen

You asked to freeze *"final ranks and exclusions"*. The first draft froze the exclusions and
then contained no ranks and no design list anywhere, which made it undiffable against the file
you receive. Here is the submission as uploaded, in order. It should reconcile row-for-row with
`submissions/01-egfr.csv`.

| rank | design | class | family | aa | pH ratio | poses | ipSAE human | ipSAE mouse |
|---|---|---|---|---|---|---|---|---|
| 1 | `bc_s831683_mpnn8_S15D` | protein | d3acid_l65_s831683 | 65 | **5.461** | 5 | 0.7760 | 0.7637 |
| 2 | `bc_s831683_mpnn19_S15D` | protein | d3acid_l65_s831683 | 65 | **5.435** | 5 | 0.8077 | 0.7859 |
| 3 | `bc_s831683_mpnn9_S15D` | protein | d3acid_l65_s831683 | 65 | **5.428** | 5 | 0.8025 | 0.8026 |
| 4 | `bc_s831683_mpnn6_S15D` | protein | d3acid_l65_s831683 | 65 | **5.397** | 5 | 0.7803 | 0.7507 |
| 5 | `rimA02_d3_rimA_14_vhh` | nanobody | rimA02 (VHH) | 129 | **5.186** | 6 | 0.2186 | 0.4468 |
| 6 | `rimA01_r15_boltzgen_egfr_d3_rimA_20` | protein | rimA01_r15 | 150 | **4.582** | 6 | 0.5938 | 0.5668 |
| 7 | `d2c_mpnn13_S88D_serasp` | protein | d2c_101_l147 | 147 | **4.572** | 5 | 0.6031 | 0.5283 |
| 8 | `bc_d3acid_l65_s831683_mpnn11` | protein | d3acid_l65_s831683 | 65 | **4.010** | 6 | 0.7963 | 0.7838 |
| 9 | `bc_s831683_mpnn9_WT` | protein | d3acid_l65_s831683 | 65 | **3.522** | 11 | 0.7833 | 0.7829 |
| 10 | `h370_020_vhh` | nanobody | h370_020 (VHH) | 98 | **2.289** | 11 | 0.4171 | 0.7093 |

Reading conventions, so the columns are not over-read:

* **pH ratio** = predicted KD(7.4)/KD(6.5), the median over every ESMFold2 refold pose of that
  exact binder sequence, pooled across runs by **sequence** (`bin/ph_pool_by_sequence.py`). Not
  a measurement.
* **ipSAE human / mouse** = ipSAE_min on ESMFold2, median of 5 seeds, pooled by sequence from
  `analysis/01-egfr/master_rank.json`. **pH-agnostic** — it is structural compatibility, not
  affinity at either pH. For ranks 5 and 10 it is **not interpretable at all** (VHH format).
* Ranks 1–4 span **0.22 pKa units** and are not distinguishable from one another. The order
  among them is not a claim.
* Rank 9 is a **control**, not a candidate: the matched wild-type of rank 3, one residue apart.

---

# PART 2 — PRESPECIFIED ANALYSIS

## 2.1 Three separate outcomes. Not combined into one score.

1. **Expression** — binary, per design, from the platform's own readout.
2. **Binding** — KD against human EGFR, and separately against mouse EGFR. Direction: lower
   KD = stronger. Designs below the assay's quantifiable ceiling are **right-censored**, not
   zero, and are analysed as censored observations.
3. **pH selectivity** — the ratio KD(7.4)/KD(6.5). **Direction: > 1 means stronger binding at
   pH 6.5**, which is the competition's objective. A design with no measurable binding at 7.4
   has an **undefined** ratio, not an infinite one, and is reported as such.

Expression is a precondition for (2); binding at pH 6.5 is a precondition for (3). A design
failing an earlier outcome is not evidence about a later one.

## 2.2 Assay limits, stated in advance

Detection limits come from the platform, not from us, and will be recorded before unblinding.
Any design whose KD is reported as "None"/not-determined counts as **no detectable binding at
that pH**, not as a measured non-binder, and not as a KD of infinity.

## 2.3 Human and mouse assessed separately, then jointly

Human and mouse are separate endpoints. Joint success (both species) is a third, prespecified
endpoint. **Mouse may be measured at pH 7.4 only** (per Anthropic, Oct 3, unconfirmed). If so,
a design that genuinely switches off at 7.4 has an undefined mouse/human ratio, and mouse
cross-reactivity is **not assessable** for that design. This is recorded now so it cannot be
reinterpreted later.

## 2.4 Primary endpoint, comparator, minimum useful improvement

* **Primary endpoint:** the proportion of submitted designs with any detectable binding to
  human EGFR at pH 6.5.
* **Comparator, corrected 2026-10-04.** The first draft said *"10 of 11 round-1 designs tested
  showed no binding, i.e. a ~9% hit rate"*. That is wrong and it **tripled the bar**: the
  collection holds **10 de novo designs, all with no binding detected**, plus **human EGF**,
  which is the collection's positive control and the native agonist — not a design. So the
  comparator is **0 of 10 expressed de novo designs with binding detected**, **exact
  (Clopper–Pearson) two-sided 95% interval [0, 0.308]**, from `proteinbase.com/collections/egfr-round1-second-submission` (local copy
  `data/proteinbase/egfr_round1_second.csv`). CONTROL-TABLE §2 described the same file correctly;
  this document did not.
* **Minimum useful improvement:** prespecified as **a hit rate whose exact 95% interval excludes
  the comparator's upper bound of 0.308**, computed on **families (n=5)** and not on designs
  (n=10) — see §1.3. With 5 clusters, of which one holds 6 designs, that bar is demanding and we
  are stating so in advance rather than discovering it in December. Anything less is reported as
  **inconclusive**.
* Every proportion reported with an **exact (Clopper–Pearson) two-sided 95% interval**,
  **clustered by sequence family**, and an explicit **inconclusive** category. *Correction: an
  earlier draft specified Wilson intervals and attributed the choice to PK. He specified "exact
  two-sided 95% intervals"; Wilson was our word, not his.*
* **Discrimination, for the December outcome analysis only:** **AUROC and average precision,
  the latter reported relative to target prevalence**, computed **within target** and reported
  target-specific **before** any equal-target summary — as asked. An earlier draft removed these
  entirely while correcting a separate error (an ROC-AUC quoted off a control panel with a single
  positive, §2.5). That over-corrected: the control panel cannot support an AUROC at n=1
  positive, and the December outcome analysis both can and should carry one. **The two are
  different objects and the distinction is now explicit here so it is not collapsed again.**
* **Falsification rule and minimum useful effect are committed above, before outcomes**, and an
  imprecise result is reported as **inconclusive** rather than as a negative.

## 2.5 Three assessments kept separate

PK: *"Keep software correctness, control recovery and prospective predictive performance as
separate assessments."*

1. **Software correctness** — does `ipsae_min.py` reproduce the pinned reference on the
   fixture set? Independent of any biology.
2. **Control recovery** — does the instrument rank the 10 measured non-binders below the
   **one** measured binder? This is answerable *now*, before any new experiment. (The first
   draft wrote "the measured binders", plural, on a panel with a single positive — the exact
   weakness you raised as `arms-backlog §2a` weeks earlier.)

   **ANSWERED 2026-10-04, 17:00 EDT, before you read this. Recorded here rather than left for
   you to find, because §2.7 names it as a falsification criterion and on one leg it fires.**

   Panel: **1 measured binder** (human EGF, KD 55 nM, n=15 runs) and **10 measured non-binders**.
   With one positive no discrimination estimate is available, so we report ranks and deliberately
   report no AUROC.

   * **Human leg: 8 of 10 non-binders rank below the binder. 2 rank ABOVE it.** The top scorer in
     the whole measured panel is a measured non-binder (`EXPNEG_gitter-yolo10`, 0.5893 vs EGF's
     0.3549), and it also reads a 5.27× pH ratio. **As stated in §2.7, this fires the criterion.**
   * **Both species required: the binder outranks all ten.** Both molecules that beat it on human
     are at **exactly 0.0000 on mouse, 5 of 5 dead seeds** — no interface, not a near miss. So the
     separation rests on a mechanism rather than on a margin a larger panel could erode.
   * 8 of the 11 molecules are one group's series and we have counted them as 8 independent
     negatives. The effective n is nearer 3 families. **Not fixed** — open question in the cover.

   **Our reading, offered rather than asserted:** the instrument fails on a single-species read
   and survives on the dual-species read that the submission actually uses. Whether that counts
   as passing §2.7 is a judgement we would rather you made than we did, which is why the result
   is here with its own caveats attached. Full table: `outbox/CONTROL-TABLE.md` §6; artifact
   `analysis/01-egfr/control_recovery.json`.
3. **Prospective predictive performance** — do our scores predict the December outcomes?
   Only this one requires the experiment, and only this one speaks to the thesis.

Failing (1) or (2) does not invalidate (3) and vice versa; they are reported separately.

## 2.6 Selection bias — what we will and will not claim

Preregistration does not undo repeated selection from one pool. Therefore:

* **Inference is restricted to the designs actually tested.** Untested designs are not
  nonbinders. No weighting scheme can recover groups that had zero probability of selection.
* We do **not** claim our ranking generalises to the pool. We claim only what the tested
  designs show.
* **If testing allocation permits**, we request that a small prespecified sample be drawn
  beyond the top ranks — random or score-stratified — with selection probabilities recorded.
  Without it, the rank-vs-outcome relationship is unidentifiable above the cut.

## 2.7 What would falsify our thesis

Stated before outcomes, so it cannot be moved afterwards:

* **The pH thesis fails** if no submitted design shows a KD(7.4)/KD(6.5) ratio > 1 with a CI
  excluding 1.
* **The instrument fails** if ipSAE_min does not separate the 10 measured non-binders from the
  **one** measured binder, or if the measured hit rate among our high-ipSAE designs does not
  exceed that among our low-ipSAE designs. **Status: fired on the human leg, survives on the
  dual-species leg — see §2.5, answered 2026-10-04 before outcomes.**
* **The mechanism fails** if designs engaging H433 do not show higher pH ratios than designs
  that do not.

---

# PART 3 — WHAT THIS DOCUMENT DOES NOT COVER

It does not retrospectively legitimise the exploratory phase. The honest summary of that phase:
one pool, screened repeatedly, with thresholds set after seeing the data, against a negative
control that turned out to be 81% of the target's native agonist. The December analysis should
treat every pre-October-4 number as hypothesis-generating.
