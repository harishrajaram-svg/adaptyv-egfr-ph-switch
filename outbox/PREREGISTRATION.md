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
| complexes folded and ipSAE-cached | **5,368** |
| submitted | 20 |

The universe is the 1,948. Designs outside it (e.g. BindCraft trajectories never sequenced)
are **not** eligible and are not nonbinders; they were never candidates.

## 1.2 Versions, pinned

| component | version |
|---|---|
| this repo | `76c5b0c0920fde7be2f4b8324ca9be1db8b0f304` (+ uncommitted working tree; see fixtures README) |
| ipSAE reference | `ipsae/ipsae.py` v4, Dunbrack, commit `6174cf9e71cb1bd660cc805856a18c4871a6dec3` |
| structure model | ESMFold2 via Modal, loops=10, steps=68 |
| pKa | PROPKA **3.5.1** |
| protonation model | Proton-PottsMPNN `potts_v6_afdb_edge_his0.3_acid0.06/epoch-0125.ckpt`, 21,494,143 bytes |
| H-bonds | HBPLUS v3.06 (McDonald), compiled locally |
| cutoffs | ipSAE PAE 10 / dist 10; pH 6.5 vs 7.4; RATIO_BAR 1.20; MIN_N 5 poses |

## 1.3 Sequence families

**19 families among the 20 submitted designs.** One multi-member family:
`d2d_101_l120_s869126_mpnn1` and `_mpnn5` (same backbone, different MPNN sequence).
All December analysis must treat that pair as **one family**, not two independent tests.

Control families (binder chain): EGF-derived (2 runs, 1 molecule), cetuximab scFv (4 runs,
1 molecule), VHH (2 runs, 1 molecule). **Three molecules, eight runs.**

## 1.4 Selection history — every narrowing applied to the pool, in order

1. generator-native `iptm` (BoltzGen's own filter, top-30 of each arm retained)
2. H433 carboxylate contact < 4.0 Å
3. the 2.4–3.2 Å "predictive band"
4. pH ratio ≥ 3.0, then ≥ 1.5, then ≥ 1.20
5. ipSAE_min ≥ 0.1493 (ECD) / 0.2218 (d3) — **both now retired**
6. n ≥ 5 poses for a ratio to rank
7. novelty: sequence identity < 30% to any PDB entry (organiser hard gate at upload)

Steps 1–5 were applied, revised and re-applied to the **same pool**, with thresholds chosen
after seeing the data. Three designs were promoted and later withdrawn on re-measurement.
Surviving designs are therefore selected on the same measurements December would analyse.

## 1.5 Exclusions, with reasons

| design | reason |
|---|---|
| `rimA02/d3_rimA_14` | 77.5% identity to `7na9_D` — organiser novelty gate |
| `h370_2site/2site_009` | 74.5% identity to `4tqe_L` — organiser novelty gate |
| `h370_only/h370_020` | 41.3% identity to `7a4t_A` — organiser novelty gate |
| 3 VHH-format designs | **inadequately assessed** by this pipeline (0/2 on VHH positives), not rejected on score |

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
* **Comparator:** the published Adaptyv de novo hit rate on EGFR from prior rounds. The
  reference point we hold is **10 of 11 round-1 designs tested showed no binding**, i.e. a
  ~9% hit rate on that collection (`proteinbase.com/collections/egfr-round1-second-submission`).
* **Minimum useful improvement:** prespecified as **beating that comparator's point estimate
  with a 95% CI excluding it.** Anything less is reported as inconclusive.
* Every proportion reported with a 95% confidence interval (Wilson), clustered by sequence
  family, and an explicit **inconclusive** category.

## 2.5 Three assessments kept separate

PK: *"Keep software correctness, control recovery and prospective predictive performance as
separate assessments."*

1. **Software correctness** — does `ipsae_min.py` reproduce the pinned reference on the
   fixture set? Independent of any biology.
2. **Control recovery** — does the instrument rank the 10 measured non-binders below the
   measured binders? This is answerable *now*, before any new experiment.
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
  measured binders, or if the measured hit rate among our high-ipSAE designs does not exceed
  that among our low-ipSAE designs.
* **The mechanism fails** if designs engaging H433 do not show higher pH ratios than designs
  that do not.

---

# PART 3 — WHAT THIS DOCUMENT DOES NOT COVER

It does not retrospectively legitimise the exploratory phase. The honest summary of that phase:
one pool, screened repeatedly, with thresholds set after seeing the data, against a negative
control that turned out to be 81% of the target's native agonist. The December analysis should
treat every pre-October-4 number as hypothesis-generating.
