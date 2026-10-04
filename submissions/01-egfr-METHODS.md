# Methods — EGFR pH-switch binder design (Challenge 1, Track 3)

**Objective.** A de novo binder to human EGFR whose affinity is higher at pH 6.5 than at
pH 7.4, cross-reactive with mouse EGFR.

**Scale.** 1,948 designs generated across 46 arms; 2,009 distinct binder sequences scored on
the ranking instrument; 1,944 pH-gated on generator poses; 2,757 poses across 131 run
directories re-gated on independent refolds; 1,693 scored by a third, protonation-aware model;
238 levelled against PDB for novelty; **27 control molecules folded, of which 11 have published
experimental outcomes.** All compute on Modal (ESMFold2, BoltzGen, BindCraft) and locally
(PROPKA, Proton-PottsMPNN, FoldSeek).

This document leads with what did not work and with the errors we found in our own instrument,
because those are the results we are most confident in. If you read one section, read **§4.4**:
it is the only place where our ranking metric is graded against molecules whose binding was
measured in a laboratory, and the highest-scoring molecule in that panel is a measured
non-binder.

**Repository.** <https://github.com/harishrajaram-svg/adaptyv-egfr-ph-switch>. Each correction
in §7 is a commit rather than a rewrite, and the 20-design build this submission cut in half
(§11) is in the history.

**What is and is not reproducible from it, stated precisely, because an earlier version of this
line claimed "everything here is reproducible" and that is false.** The analysis code, the
submission, every artifact under `analysis/`, and the full ranking table are published. The
**6,011 cached pose outputs under `runs/` are not** — they are gitignored and run to tens of
gigabytes. Consequence, and we would rather you read it here than discover it: running the
documented command `python3 bin/emit_submission_csv.py` on a fresh clone emits a complete,
plausible CSV in which **every affinity column reads 0.0000 and every `affinity_above_null`
reads `no`**, with exit code 0 and no warning, because the pose index it needs is absent. That
is this project's own signature failure — a `0.0000` meaning *absent* being read as *measured* —
reproduced inside its reproducibility claim. The emitter now raises instead of emitting a
silently zeroed file when a submitted sequence is missing from the ranking table, but it cannot
conjure poses it does not have. **To reproduce the affinity columns you need the pose cache; ask
and we will supply it. Everything else reproduces from the clone.**

---

## 1. Mechanism

We used **thermodynamic linkage**. A titratable site with pKa_free in the unbound state and
pKa_bound in the complex changes the observed affinity with pH:

    ratio = K(6.5) / K(7.4),    K(pH) = (1 + 10^(pKa_bound − pH)) / (1 + 10^(pKa_free − pH))

Two mechanisms are available:

- **Mechanism B** — a carboxylate on the *binder* raises the pKa of a histidine on the
  *target*, so the salt bridge only forms once that histidine protonates. We used this.

**How this is assessed, from the organisers' own specification** (Adaptyv, 2026-10-04): human
EGFR is measured at **both pH 6.5 and 7.4**; mouse EGFR at **pH 6.5 only**. So the switch is
graded on the human leg and cross-reactivity on the acidic leg, where a working switch should
bind. Our columns map onto that without adjustment: the pH ratio is pooled over **human-leg**
refold poses only, and both affinity columns are the binding-competent state. An earlier
project note recorded mouse as being measured at 7.4 only and concluded that cross-reactivity
and pH selectivity were therefore in tension; the organisers retracted that premise and the
tension does not exist.

**The switching residue is the target's own H433 — a native domain III residue, conserved in
mouse, and not an affinity tag.** This is worth stating because the assayed construct carries a
C-terminal His tag, and the organisers have said explicitly that a binder engaging the tag
"might look pH-selective but it would bind to anything with a His tag", and that in-silico
evidence will be weighted more heavily to catch it. Across our whole pool of 50 designs that
switch at n ≥ 5, the engaged site is **H433 in 48, H370 in 1 and H383 in 1** — zero tag. All ten
submitted designs switch on H433. We designed against a tag-free crystal structure, so the tag
was not available to optimise against even accidentally.
- **Mechanism A** — the inverse: histidine on the binder, carboxylate on the target.

**The ceiling, and why it set the strategy.** The general one-proton upper bound over this pH
pair is **7.94×** per site; with H433's measured free pKa of 6.22 the single-site model spans
**0.699× to 5.55×**. Both figures were independently confirmed by our collaborator's review.
Worth noting that our empirically measured non-switch floor, **0.702×** over n=570 designs
(§6), lands on the analytic lower bound of 0.699× — the floor is the model's own minimum, not
a property of those designs. For EGFR's H433, whose free pKa we measure at 6.22,
the attainable ceiling is **5.55×**. Over a 0.9 pH-unit window two sites give 63× and three
give 500×; a single site cannot exceed 7.94× by any design.

A full PROPKA census of all 17 histidines in the EGFR ectodomain (6ARU, apo):

| site | pKa_free | own ceiling | ring-N distance to H433 | histidine in mouse |
|---|---|---|---|---|
| H433 | 6.22 | 5.55× | — | yes |
| **H370** | **4.93** | **7.76×** | **8.5 Å** | yes |
| H418 | 4.89 | 7.78× | 25.2 Å | yes |
| H383 | 6.38 | 4.95× | 26.1 Å | **no (Arg)** |
| H358 | 6.26 | 5.41× | 26.2 Å | yes |

H433 + H370 is the only pair close enough for one binder to bridge: a two-site ceiling of
**43.1×**. Everything in §3 is an attempt to reach it.

---

## 2. Pipeline

    backbone + sequence (BoltzGen, BindCraft)
      → free geometry screen: nearest binder carboxylate to a target histidine (gemmi, no GPU)
      → complex prediction (ESMFold2, 5 seeds) → ipSAE_min, min over alignment directions
      → pH gate: PROPKA on the complex AND on the same file with the binder deleted in place
      → knockout control: the causal carboxylate truncated to Ala, same coordinates
      → express QC, novelty (FoldSeek vs PDB)

Two choices are worth stating.

**The free leg is computed per design, from the design's own coordinates.** Our first pass
took pKa_free from a single apo structure on the reasoning that the target is rigid. It is
not: the generator repacks the target per design, and H370's own suppressor ASP344 sits
2.4–3.8 Å from its ring. A rotamer shift in ASP344 then appears as a "switch" in designs
whose binder is 27 Å away. That artifact turned 10 apparent switches into 1. Deleting the
binder from the *same file* holds every other coordinate fixed, so the difference is the
binder and nothing else.

**The pH gate runs on an independent refold, not the generator's pose.** See §5.

---

## 3. What we measured

### 3.1 Cropping the target fixes epitope selection. Confirmed four times.
Against the full 609-residue ectodomain, 4/30 designs touched the declared patch. Against a
170-residue domain III crop, **27/30 and 29/30** did. The generator was not ignoring the
binding-site declaration; it was taking a better offer, because much of the time another
domain was more designable. Remove the alternative and it complies.

### 3.2 H370 is reachable by accident and not by instruction.
Across all 1,944 designs, a binder carboxylate lands within 4.0 Å of H370 in **23 (1.5%)**
without being asked. We then ran two arms that declared the H370 patch explicitly:

| arm | H370 engaged | H433 engaged |
|---|---|---|
| H370 patch alone, n=120 | **0 (0.0%)** | 16 (13.3%) |
| H370 + H433 patches, n=120 | 3 (2.5%) | 16 (13.3%) |
| incidental baseline, n=1,944 | 23 (1.5%) | 5.4% |

Declaring the patch produced **fewer** engagements than not declaring it, and both arms gave
the *identical* H433 rate — the declaration changed nothing about where designs land. The
best site for 89/111 and 90/113 designs in those two arms was H433, not the H370 they were
told to hit.

### 3.3 H358 is unreachable, and we spent 60 designs learning it.
0/60 engaged, closest approach 13.6 Å. It is 26.2 Å from H433. We selected it before
measuring the pairwise distances; the distance matrix would have ruled it out for free.

### 3.4 We could not build the two-site route — but it exists. See §8.1.
**This section previously read "the two-site route is closed." It is not; we failed at it.**
G532, a published antibody, achieves **13.26×** on human EGFR using carboxylate contacts with
**H433 and H370** — this exact pair (§8.1). Everything below is an accurate account of our own
failure and should be read that way.

3 of 1,944 designs engage H433 and H370 simultaneously, all by accident. The best physically
composed multi-site product in the entire pool is **6.58×** — below the 7.94× single-proton
bound and 15% of the 43.1× two-site ceiling. Of 23 designs that reached H370 incidentally,
**9 of 10 that switched did not bind.** H370 has 153 heavy atoms within 10 Å of its ring
against H433's 57; reaching into that cleft appears to cost the interface area binding needs.

### 3.5 RETRACTED — "binding and switching trade off" was our central result and it is sign-reversed.

An earlier version of this section read: *"This is the central result. Designs with the tightest
carboxylate–histidine contacts switch most and bind least. It reproduces on the generator pose,
on the refold, and on a trained protonation-aware model that shares no machinery with either."*
It was four lines of prose with no table, no n and no coefficient. **Computed, it comes out with
the opposite sign.**

Across the 115 designs with n ≥ 5 pooled refold poses and an affinity reading, controls excluded,
pH ratio against human ipSAE_min:

| subset | n | Pearson | Spearman |
|---|---|---|---|
| all designs with n ≥ 5 | 115 | **+0.342** | **+0.273** |
| those that switch (ratio ≥ 1.20×) | 49 | +0.249 | +0.247 |
| those that bind human ≥ 0.50 | 34 | +0.199 | +0.170 |

A trade-off requires a negative coefficient. All three slices are positive: in this pool, designs
that switch more tend to bind slightly **better**, not worse. **The claim is withdrawn.** We have
no pool-wide trade-off, and the sentence asserting one as the project's central result should
never have been written without the table underneath it.

**What survives, and it is narrower and site-specific.** The trade-off is real at **H370** and
only there. Of the 23 designs that reached H370 incidentally, **9 of the 10 that switched did not
bind** (§3.4); the one design built deliberately to put a carboxylate on H370, `S60D`, returned
1.11× — inside PROPKA's noise — while dropping human ipSAE from 0.654 to **0.012** (§3.7). H370
carries 153 heavy atoms within 10 Å of its ring against H433's 57, so reaching it costs the
interface area binding needs. That is a statement about one cleft, not about the pool, and it is
what §3.4's conclusion actually rests on.

The three-predictor agreement claimed in the retracted text was never computed on affinity at
all; §5 reports what the three pH predictors agree on, which is a different question.

### 3.6 Contact distance predicts switching — on the generator pose only.
Measured on the **1,584** gated poses with a resolvable carboxylate–histidine distance at
H433 ∪ H370 (of 1,944 gated in total; the remainder have no carboxylate within range to measure).
An earlier version of this header said "all 1,944", which the table's own rows contradict —
they sum to 1,584. Generator pose:

| carboxylate–ring N | n | switch ≥1.20× |
|---|---|---|
| ≤2.8 Å | 40 | **65.0%** |
| 2.8–3.2 Å | 49 | **61.2%** |
| 3.2–4.0 Å | 136 | 19.9% |
| 4.0–6.0 Å | 193 | 0.5% |
| >6.0 Å | 1,166 | 0.3% |

On the refold the same relationship collapses to 19.0% vs 9.7%. **The generator-pose pH gate
is largely a readout of the generator's own geometry** — self-consistent, and weak as
validation.

---

### 3.7 One mutation converts a binder into a switch, and it replicates four times.

This is the only causal result in the project, and the only one where we changed one thing and
measured the consequence.

`d2c_mpnn13` binds and does not switch (0.70×). A single Ser→Asp at position 88 gives
**4.57×** over 5 poses. Truncating that Asp back to Ala in the same coordinates returns
**0.72×**. Knockout effect **+5.20 pKa units**, the largest in the project.

The replication came off a different backbone. On `d3acid_l65_s831683` — the best-binding
BindCraft backbone we have — four independent ProteinMPNN sequences each took the same single
Ser→Asp at position 15:

| sequence | WT ratio | S15D ratio | human ipSAE WT → S15D | mouse ipSAE WT → S15D |
|---|---|---|---|---|
| mpnn6  | 3.85× | **5.40×** | 0.789 → 0.780 | 0.769 → 0.751 |
| mpnn8  | 3.15× | **5.46×** | 0.821 → 0.776 | 0.807 → 0.764 |
| mpnn9  | 3.52× | **5.43×** | 0.783 → 0.803 | 0.782 → 0.803 |
| mpnn19 | 3.68× | **5.43×** | 0.828 → 0.808 | 0.815 → 0.786 |

Four for four in the same direction, at a small and consistent cost in predicted affinity
(−0.009 to −0.045 human, −0.018 to −0.043 mouse; the mpnn9 pair is the exception and rises).

**A correction to how this was previously stated.** An earlier version of this section said
these were "simultaneously the four highest-affinity and the four highest-switching designs in
the pool." **That is false on both halves**, and it was the sentence used to justify spending six
of ten submission slots on one backbone, so it matters. Three designs with n ≥ 5 out-switch them
(5.819×, 5.630×, 5.527×) and three out-bind them on the human leg (0.860, 0.828, 0.821) — two of
the latter being these designs' own unmutated siblings. What is true, and sufficient:

> Among the 12 designs that bind **both** species at ≥ 0.75 with n ≥ 5 poses, these four are the
> four highest-switching, and the gap to the fifth is **5.40× against 4.01×**.

That is a joint claim with its conditions attached, and it is the honest version. Taken
unconditionally on either axis alone, they are not first.
`bc_s831683_mpnn9_WT` is submitted alongside its S15D so the comparison gets made in the wet
lab instead of inferred from our gate. If the S15D designs switch and the WT does not, the
mechanism stands independently of every threshold in this document.

**Two things this result does not support, stated because they are easy to read into it.**

*The four are not distinguishable from each other.* Inverting the linkage equation of §1 at
pKa_free = 6.22, those ratios imply pKa_bound of 8.89, 9.11, 8.98 and 8.98 — the whole 5.40–5.46
spread is **0.22 pKa units wide**, far inside PROPKA's own error. They are four measurements of
one quantity, and the order among them is noise. The WT→S15D step, by contrast, is +1.2 to
+1.7 pKa units and sits well outside it.

*The implied shift is large.* 97–98% of the 5.55× single-site ceiling requires raising H433's
pKa by **+2.7 to +2.9 units**. Near the ceiling the ratio saturates, so the number is
insensitive to how wrong that shift is — which protects the result from small pKa errors and
equally stops it being evidence that the shift is right.

Finally, four MPNN sequences on one backbone are four sequences, not four independent designs,
and the ratios come from the same PROPKA-on-refold gate as everything else in this document.

---

## 4. The instrument: what it gets right, and what it gets wrong

Our ranking metric is **ipSAE_min** on ESMFold2, the minimum over both alignment directions
of one interface. Two findings about its calibration matter more than any design result.

### 4.1 The negative control is EGF-derived, and its activity is unknown.
`NEG_nonbinder`, the molecule that set our "does it bind" floor, is **81% mature human EGF**:
43 of 53 residues verbatim, all six cysteines in correct linear order, two insertions in the
loops. Human EGF's high-affinity contacts are on domain III — the surface every design
targets — and on this assay platform EGF's own measured KD is **55 nM** (median of 15 runs,
range 27–795 nM; §4.4). An earlier version of this document said ~2 nM, which was quoted from
memory rather than from the assay and is corrected here. The conclusion does not change: EGF
is a real binder, so an 81%-EGF molecule cannot serve as a negative control. There is no
provenance for this sequence anywhere in the project.

**What this does and does not establish.** It is sufficient to stop treating the molecule as
an established nonbinder. It does **not** establish that the altered sequence retains EGF's
affinity or agonist activity, and the higher domain III score (0.2218 crop vs 0.1493 ECD) is
not proof of binding either — changing the target crop changes prediction and scoring
behaviour. Both EGF-derived controls are therefore relabelled **activity-unknown**. Every
claim that rested on their supposed negative status is withdrawn; the original numbers are
preserved for the audit trail.

**Consequence for this submission.** We retired the 0.1493 / 0.2218 thresholds rather than
replacing them with another universal number. No design is excluded on affinity. Affinity is
reported as a continuous score with an explicit "above the computational null" flag, which
means exactly that and not "binds EGFR".

### 4.2 The instrument is blind to nanobodies.
The two "independent" positive controls are **95.9% identical** — one cetuximab-derived scFv
measured twice; their agreement had been cited as instrument validation. The one genuinely
different positive, a real 5 nM VHH, scores **0.0000 with 5 of 5 seeds dead**. Measured
sensitivity is **2/4 on true positives and 0/2 on nanobody format.** This panel is small and
partly redundant, so it does not estimate general sensitivity reliably — but it is enough to
say that any VHH design is **inadequately assessed by this pipeline**. Our three VHH-framework
candidates are marked as such rather than rejected for low scores. **§4.5 closes this: against
framework-preserving CDR decoys the blind spot is total, not partial.** (Two are separately
ineligible on the organisers' novelty rule at 77.5% and 74.5% sequence identity to solved
structures — an upload-time hard gate, distinct from our scoring.)

### 4.3 Control-derived bars are construct-specific.
Bars must come from controls folded against the *same* target construct:

| construct | nonbinder | weakest real binder |
|---|---|---|
| full ECD (621 aa) | 0.1493 | 0.6224 |
| domain III crop (170 aa) | **0.2218** | **0.4005** |

The strict ECD bar of 0.6224 is unreachable on the crop, where a construct-matched 1 nM
binder only reaches 0.4183.

---

### 4.4 Control recovery: the instrument graded against measured outcomes.

Adaptyv published 11 molecules expressed and measured against EGFR on the assay platform this
challenge uses — human EGF, and **10 designs with no binding detected**
(`proteinbase.com/collections/egfr-round1-second-submission`). We folded all 11 against both
target constructs, 5 seeds each, and scored them with the same ipSAE_min that ranks our
designs. **This is the only test in the project where the instrument is graded against an
experimental outcome rather than against itself.**

| molecule | measured | ipSAE_min human | ipSAE_min mouse |
|---|---|---|---|
| Human EGF | **KD 55 nM** (n=15) | 0.3549 | 0.3700 |
| gitter-yolo10 | no binding | **0.5893** | 0.0000 |
| gitter-yolo9 | no binding | **0.4326** | 0.0000 |
| deepsatflow-design7 | no binding | 0.2168 | 0.1823 |
| gitter-yolo5 | no binding | 0.0264 | 0.0108 |
| gitter-yolo4 | no binding | 0.0109 | 0.0000 |
| gitter-yolo3 | no binding | 0.0000 | 0.0121 |
| gitter-yolo2, 6, 7, 8 | no binding | 0.0000 | 0.0000 |

**What this panel can and cannot support, before the numbers.** It is **10 negatives and one
positive.** With a single positive there is no discrimination estimate to be had: every summary
statistic reduces to *where that one molecule ranks*, and an ROC-AUC computed from it would be
that same single comparison restated ten times. We report ranks and we do not report an AUC.

**On the human leg, 8 of 10 measured non-binders rank below the measured binder and 2 rank above
it.** The highest-scoring molecule in the whole measured panel is a measured non-binder:
gitter-yolo10 at 0.5893 against EGF's 0.3549.

**The cross-reactivity requirement is what rescues it.** Both molecules that outrank the binder
score exactly **0.0000 on the mouse leg**. Scored the way the submission is scored — requiring
both species — **the one measured binder outranks all ten measured non-binders.**

The reason that works here is a mechanism and not a curve: the two false positives are not
*nearly* excluded by the mouse leg, they are at **exactly zero on five of five seeds**, which is
what this instrument returns when it finds no interface at all. A filter that depends on a
margin could erode with a larger panel; this one does not depend on a margin. The dual-species
criterion came from the brief, not from us, and it is the only specificity filter in this
pipeline that measured data supports at all — which is a statement about the *absence* of
evidence for the others, not a validation of this one on n=1.

**gitter-yolo10 is the entire problem in one molecule.** Pooled over 5 refold poses it reads a
**5.27× pH ratio** — **8th of the 132 molecules eligible to rank** on our primary objective —
alongside a human ipSAE above every positive control we have. (An earlier version of this
document said "rank 8 of 2,009", pairing a rank computed among molecules with ≥5 poses with the
denominator of every sequence ever scored. That overstated it by about 15× and the derived
"top 0.4%" was wrong; 2,009 is the scored universe, 132 is the rankable one.) A human-leg-only pipeline would have submitted a
molecule already measured not to bind. One number excludes it, and it is mouse 0.0000.

So the pH ratio does not discriminate binders, and we report it as the primary objective
anyway because it is the objective. **A 5× switch is not evidence of binding, on measured
data, in this pool.** Any ranking that reads the ratio without the affinity columns beside it
is reading a number that a known non-binder scores in the top 6% of everything eligible to rank.

Three limits. The positive class is **one molecule** of 53 aa, shorter than every design here,
so nothing in this section estimates sensitivity — the strongest honest reading is "the one
measured binder outranks all ten measured non-binders when both species are required."
"No binding detected" bounds affinity from below; it does not prove no interaction. And ten
negatives do not characterise a tail, so the 0.5893 we found is a floor on how high a measured
non-binder can score here, not a ceiling.

### 4.5 The matched null is degenerate, so the screening flag does almost nothing.

12 shuffle nulls — 4 designs × 3 independent shuffles, composition and length preserved
exactly, 5 seeds nested within each shuffle — score **0.0000 on both legs with 5 of 5 dead
seeds, all twelve of them.** A separate, larger shuffle run (22 molecules, 110 poses) agrees:
**20 of 22 are dead on all 5 seeds and the highest median in it is 0.0110.** Across 34 shuffled
sequences the null does not produce a single score worth a threshold. This was predicted on
review: full-sequence shuffling destroys the fold, so the null bounds the instrument's noise
floor at exactly zero, and a "95th percentile of a matched null" flag admits anything above
0.0000. We still report `affinity_above_null` in the
submission, and it carries almost no information. **The 10 measured non-binders are the only
negative class with a usable tail, and that tail reaches 0.5893.**

4 framework-preserving CDR decoys — the real 5 nM VHH scaffold with its three CDRs randomised,
located by conserved anchors rather than fixed indices — score 0.0000, 0.0117, 0.0000, 0.0000.
The real 5 nM VHH scores 0.0000 as well. **The instrument cannot separate a validated nanobody
from a randomised-CDR decoy on that nanobody's own framework.** §4.2 called the nanobody blind
spot partial on a 0/2 panel; on a matched decoy class it is total. The two VHH-format rows in
this submission are reported on that basis and are not claimed to be scored.

---

## 5. Three pH predictors, and how much they agree

| predictor | what it computes |
|---|---|
| PROPKA on the **generator** pose | pKa shift |
| PROPKA on an independent **ESMFold2 refold** | pKa shift |
| **Proton-PottsMPNN** ΔΔG | ΔG_bind(His protonated) − ΔG_bind(His neutral) |

**Our two PROPKA readings agree at odds ratio 1.55** — weak concordance. We previously wrote
that at most one of them measures something real; that does not follow, and we withdraw it.
Two noisy predictors can disagree while both carry signal. Of 21 designs both call switches,
only 14 agree on *which* histidine.

We therefore added a third, from a different family. Proton-PottsMPNN (Jacobsen, Ovchinnikov
et al., 30 Sep 2026) represents protonated and deprotonated His/Asp/Glu as distinct sequence
tokens. Its binding-energy term gives a direct observable rather than an inferred pKa:

    ΔΔG = G_bind(H433 = HIS-P) − G_bind(H433 = HIS-S)

computed over all 1,693 scored complexes by overriding one target-chain token and holding
every other residue's protonation state fixed.

**Potts vs the PROPKA refold gate.** Our first analysis split at the pool median, a threshold
chosen after seeing the data, and reported OR 3.70. Re-done at the **prespecified, physically
meaningful** threshold — the sign of ΔΔG, i.e. does protonation favour the complex at all:

|  | Potts ΔΔG < 0 | ΔΔG ≥ 0 |
|---|---|---|
| PROPKA switch (n=182) | 86 | 96 |
| PROPKA no-switch (n=1,475) | 183 | 1,292 |

**odds ratio 6.32, 95% CI [4.55, 8.79].**

An important correction to our own earlier reading. Median ΔΔG is +0.11 for PROPKA-called
switches against +3.20 for non-switches — **both positive**. That is a lower relative penalty,
not demonstrated acid-favoured binding. Only **47% of PROPKA-called switches, and 17% of the
whole pool, actually have ΔΔG < 0.**

Three further limits. A learned Potts score is not automatically a calibrated physical free
energy, so the units are model-internal. Fixing H433's protonation state does not give the
population-weighted affinity change between pH 6.5 and 7.4. And the 1,693 complexes are not
1,693 independent designs — they include repeated backbones and sequence families, so the
interval above is optimistic.

Agreement between two methods is evidence, not validation. Neither has experimental ground
truth on this target.

The controls validate it independently of us: all six positive controls fall in the lower
half of the design distribution, and both EGF-derived "nonbinders" fall in the **top 14%** —
a model that never saw this target says protonating H433 is maximally bad for an EGF-like
complex, which is exactly right for a neutral-pH agonist.

---

## 6. Reproducibility

**76% of the switches in our pool rest on a single pose** — 164 of 217 designs reading ≥ 1.20×
have n = 1. An earlier version of this section said **93%**, which is the fraction of the *whole
design pool* at a single pose (1,827 of 1,982), not of the switches. The switch-specific figure
is the relevant one and it is lower. Of the 12 single-pose switches that later received five
poses, **5 of 12 fell below the previous threshold on re-evaluation.**
That establishes instability under re-measurement. It is **not** a biological false-positive
rate: there is no experimental ground truth, and the threshold itself was derived from a
control we have since retired. We therefore require **n ≥ 5 poses** before a ratio may rank a design, and we
report the per-pose spread, because a median of 1.5 over poses of
[0.78, 1.01, 1.03, 1.99, 4.39, 4.81] is not a measurement.

**The non-switch floor is not zero.** Designs that contact a histidine but have no
carboxylate within 6 Å have a median ratio of **0.702** (n=570, only 0.4% above 1.20). This
is the desolvation signature of a binder sitting on a histidine — a constant of the method,
not a per-design result.

---

## 7. Errors we found and corrected

Listed because they bound the confidence of everything above.

| | |
|---|---|
| Apo-reference free leg | 10 apparent switches → 1 (§2) |
| ipSAE_min taking the min across chain *pairs*, not alignment *directions* | wrong score on 3-chain complexes |
| `best_ratio` as a max over all 5 histidines | untouched sites return exactly 1.000, a hard floor on 1,263 of 1,584 designs |
| active-site selection by largest ΔpKa | a histidine can shift 3 units entirely below pH 6.5 and change nothing |
| target construct identified by chain length | silently dropped 32 designs |
| chain order assumed binder-first | **29 of 69 run directories have the target in chain A** |
| `product` composed over helping sites only | selection on the outcome; manufactured our only claim above the 7.94× ceiling |
| joins on design basename | `boltzgen_egfr_d3_00` exists in four arms; **337 designs were never scored at all** |
| control bars matched by `startswith("POS_")` | dropped all 20 construct-matched controls |
| a 20-design confirmation run | folded a 125-aa molecule for a 129-aa design and "killed" a live candidate |
| novelty tested as a single TM < 0.50 bar | missed the "high structural similarity ALONE" clause; reported 0 of 238 passing where 114 clear (§9) |
| the refold pH gate run on a third of the run directories | the four-pair S15D replication of §3.7 sat unscored for a day; 170 switches missed, 179 generator false positives left standing |
| controls folded but never scored | `expctrl`/`ctrl2` wrote raw PAE matrices and no ipSAE output, so §4.4 — the only measured-outcome test in the project — was unavailable until the scoring step was backfilled |
| **§3.5's "central result" asserted in prose, never computed** | computed, the correlation is **+0.34**, the opposite sign. Retracted (§3.5) |
| **affinity columns taken from one run directory, chosen by name sort** | 5 of 10 shipped cells matched neither the pooled median nor max; now read from `master_rank.json` (§11) |
| an ROC-AUC reported from a control panel with **one** positive | removed; §4.4 reports ranks. The weakness was already recorded weeks earlier as "the gate is n=1 positive" and shipped anyway |
| a rank among 132 rankable molecules quoted against a denominator of 2,009 | overstated by ~15× (§4.4) |
| "the four highest-affinity **and** four highest-switching designs" | false on both halves, and it was the justification for 6 of 10 slots on one backbone (§3.7) |
| a measured pH-calibration ladder folded and never scored | the G532 series, the only molecules here with a published pH ratio, is absent from every analysis file (§12) |

We also retracted two published-in-log claims on re-measurement: a length effect that
disappeared at n=260, and a "pinning is 15× worse" conclusion that reverses on a like-for-like
comparison (pinned+cropped is 15/60 = 25%, the highest-reaching configuration we ran).

---

## 8. Context from the literature

Our window sits **above every published pH-switch transition.** Proton-PottsMPNN recovered
hundreds of pH-dependent binders from 8,407 designs with transition pH **4.0–5.8**. The Baker
lab's designs compare pH 7.4 to 5.4 and reach up to 1000× with three or more interface
histidines, from 4 experimental hits in 12,000 designs. The only published result near our
window is a **>2× at pH 6.0**.

### 8.1 RETRACTED — "nobody reports a design that binds more strongly at low pH." G532 does.

An earlier version of this section said: *"Both bodies of work design binders that release in
acid; neither reports a design that binds more strongly at low pH. Our direction is the one
neither achieved."* **A counterexample was in our collaborator's review a day before we wrote
that, and we had not read far enough into the email to find it.** It is the most important
citation in this document and we missed it.

**G532** is a published antibody that binds EGFR **more tightly at pH 6.5 than at 7.4**, on
both species, by SPR:

| target | KD at pH 6.5 | KD at pH 7.4 | KD(7.4)/KD(6.5) |
|---|---|---|---|
| human EGFR | 294 nM | 3,900 nM | **13.26×** |
| mouse EGFR | 547 nM | 1,810 nM | **3.31×** |

And the mechanism is **ours**: the study proposes antibody carboxylate interactions with the
target's **H433 and H370** — the same two histidines, and the same direction of linkage, that
this entire project is built on. (Caveat supplied with the citation: those contacts rest on
modelling and mutational data, not an experimentally solved complex.)

**Three things follow, and they revise §1, §3.4 and §3.5.**

**1. Mechanism B is not self-defeating, and we should not have described it that way.** A real
molecule exploits a target-side histidine and gets 13.26×. Our §3.5 retraction already withdrew
the pool-wide trade-off claim on our own data; G532 independently removes the structural argument
for it.

**2. The two-site H433+H370 route is demonstrated, not closed.** §3.4 concluded the route was
closed because 3 of 1,944 of *our* designs engaged both sites, all accidentally, and our one
deliberate attempt (`S60D`) bought 1.11× while destroying the interface. That remains an accurate
account of **our** failure. It is not evidence the route is unavailable: **13.26× exceeds the
7.94× single-proton upper bound**, so G532 must be linking more than one proton, and the sites it
names are the pair we measured at 8.5 Å. The honest statement is that we failed to build it, not
that it cannot be built.

**3. Our window is not above every published transition.** G532 is measured at exactly our pH
pair, 6.5 against 7.4. The "above every published pH-switch transition" framing stands only
against Proton-PottsMPNN (4.0–5.8) and the Baker lab (7.4 vs 5.4).

**What this costs us.** Our best design reads 5.46× against G532's 13.26×, and G532 is a
*measured* number against our *predicted* one. It is also an antibody, where CDR loops make
multi-site placement tractable, which is the asymmetry §8.2 describes. We are reporting a
single-site design at 98% of its site's ceiling in a problem where a two-site solution is known
to exist and is 2.4× better.

### 8.2 Where the histidines sit, and why it capped us

Their histidines are on the **binder**, where they are placeable and multipliable. Ours is on
the **target**, where evolution fixed its position. That choice, made before we read either
paper, is what capped this project at a single site — and G532 shows the cap is a property of
*our* search, not of the mechanism.

---

## 9. Novelty

The organisers reject at upload on their own novelty levels, so this is a hard gate on the
submission and not a scoring preference. Our implementation of the general-protein rule had one
clause missing — **high structural similarity alone (TM ≥ 0.80) is Level 2, with no sequence
condition attached** — and in its place we tested a single TM < 0.50 bar and called everything
above it a failure. That is where the "0 of 120 pass novelty" result recorded earlier in this
project came from. It is void.

Re-levelled on FoldSeek-vs-PDB results for every design we hold structures for:

| level | n |
|---|---|
| 1 Essentially Known | 82 |
| 2 Familiar | 42 |
| **3 Partly Novel — clears the gate** | **114** |
| 4 De Novo | 0 |

**114 of 238 clear Level ≥ 3. The old single bar passed 0 of 238** — including the entire
BindCraft pool, which supplies ranks 1–4, 8 and 9 of this submission. Nothing in the pool
reaches Level 4 under the general-protein rule; every design that clears does so through the
Level 3 clause (moderate structural similarity *or* >30% sequence identity, exactly one of
them), which is worth saying plainly: **this is a pool of partly novel designs, not de novo
ones.**

Two limits we state rather than hide. Adaptyv run MMseqs2 against SwissProt, PDB, patent
sequences, the therapeutic-antibody database and PLAbDab; we search **PDB only**, via FoldSeek,
so our sequence identities are lower bounds and a design clean here may still hit a patent or a
SwissProt entry. And the antibody branch needs CDRH3 identity from ANARCI numbering, which we
do not run; we implement it from conserved framework anchors instead
(`bin/antibody_novelty.py`, self-tested).

**The two nanobody rows carry upload risk, and we are stating it rather than discovering it.**
`h370_020_vhh` clears at antibody-rule Level **4** — the only de novo design in the set — but
under the general-protein rule it is Level 2 (TM 0.914) and would be auto-rejected.
`rimA02_d3_rimA_14_vhh` is general-rule Level 1. Both depend on Adaptyv's ANARCI calling them
antibodies. If it does not, those two rows fail at the gate.

---

## 10. The exclusion ledger

**Read this section with two corrections in mind.** First, the "45 designs" below are **45 run
names and 38 distinct molecules** — the gate is name-keyed on the discard side, so seven alias
pairs are double-counted. The ledger built to catch this project's keystone trap is itself keyed
on names. Second, and more seriously, the gate reads 69 rows out of the `phgate_*.tsv` files
rather than the 2,009-row sequence-keyed `master_rank.json` that this document elsewhere calls
canonical — **and no `phgate_*.tsv` exists for the `sd`/`sd2` arms that supplied five of the ten
shipped designs.** Joined properly by sequence, **28 molecules with n ≥ 5 poses outrank the
weakest shipped design** and were never candidates, two of them outranking the best shipped
design (5.819× and 5.630×). The gate prints PASS. Those 28 are *not* accounted for below; they
are a known, unclosed gap, recorded here rather than left for a reader to find.

Our pre-submission gate (`bin/check_discards.py`) fails the build if any design that beats a
submitted one on the primary objective was never measured on the ranking instrument — that is
exclusion by proxy rather than by evidence, and it had already happened twice in this project
at scales of 14 and 337 designs. It **passes**: every design outranking a submitted row on the
pH ratio was also scored on the instrument. It then warns on **45 run names — 38 distinct
molecules — that were measured and still rejected**, and requires the reason be recorded. Here it
is.

All 45 ratios in that list are **generator-pose** readings. Joined back to the pooled refold
table by binder sequence:

| why it was rejected | n |
|---|---|
| fewer than 5 refold poses **and** the refold pool kills the switch (<1.20×) **and** binds neither species | 20 |
| fewer than 5 poses + refold kills the switch + scores on one species only | 6 |
| fewer than 5 poses + refold kills the switch | 5 |
| fewer than 5 poses + binds neither species | 4 |
| refold kills the switch + binds neither species | 3 |
| refold kills the switch | 3 |
| refold kills the switch + one species only | 1 |
| scores on one species only (`bg04_r03`, 2.47× generator → 1.94× pooled, mouse 0.0000) | 1 |
| no refold pool at all (`rank2_noC_0`) | 1 |
| **alias of a design that IS submitted** | 1 |

**35 of 45 were never confirmed on 5 or more refold poses**, which §6 makes a precondition for a
ratio to rank anything. **39 of 45 drop below the 1.20× gate once the refold pool is pooled** —
the largest, `domIII_1His_ctrl_r17`, goes from 6.60× on one generator pose to 0.67× on the
pooled refold. This is §3.6 operating as a filter rather than as a table: the generator-pose
ratio is largely a readout of the generator's own geometry, and these 45 are what that looks
like when you stop trusting it.

The last row is an artifact of our own naming, not a rejection.
`rank15_boltzgen_egfr_d3_rimA_20` and `rimA01_r15_boltzgen_egfr_d3_rimA_20` are
**byte-identical sequences** under two run names; the second is submitted at rank 6. The gate
is name-keyed on the discard side and so flags an alias of a submitted design as a rejection.
We are recording it rather than silencing it, because name-vs-sequence confusion has broken
five analyses in this project and this is the sixth instance, caught by a check instead of by
luck.

Four designs in the list keep a pooled ratio above 1.20× but bind neither species at n < 5
(`domIII_1His_ctrl_r20` 5.72×, `bg01_r18` 3.48×, `domIII_1His_r15` 2.97×). Those are the only
entries in the 45 we would want back, and what they need is poses, not an argument.

---

## 11. The submission

**10 designs, not the 20 Track 3 allows.** The primary column is the pH ratio, pooled as the
**median over every ESMFold2 refold pose of that exact binder sequence**
(`bin/ph_pool_by_sequence.py`), joined to affinity **by sequence** (`bin/master_rank.py`)
because the same molecule appears in this project under as many as three run names. Both
affinity columns are ipSAE_min on ESMFold2, median of 5 seeds. `affinity_above_null` is
`max(ipSAE_human, ipSAE_mouse) ≥ 0.2218`, a **retained legacy reporting threshold** that came
from the retired EGF-derived control — not a percentile of the matched null, which §4.5 shows is
degenerate at 0.0000 and admits of no percentile. It excludes nothing (`MIN_AFFINITY = None`)
and should be read as a provenance marker, not a claim about binding.

**One provenance correction to the affinity columns themselves.** Until 2026-10-04 the emitter
took them from a helper that selects **one run directory** per molecule ("most seeds wins, ties
broken lexicographically") and medians within it. Because the same molecule appears here under
as many as three run names, the reported affinity was therefore a function of which run name
sorted first — the keystone error of this project, committed on the submission path. Five of the
ten cells matched neither the pooled median nor the pooled max as a result; the five single-run
molecules were unaffected. All ten now come from `master_rank.json`, pooled over every pose of
the binder sequence, which is what this section already claimed they were.

| # | design | class | pH ratio | human | mouse | poses |
|---|---|---|---|---|---|---|
| 1 | bc_s831683_mpnn8_S15D | protein | **5.46×** | 0.776 | 0.764 | 5 |
| 2 | bc_s831683_mpnn19_S15D | protein | **5.43×** | 0.808 | 0.786 | 5 |
| 3 | bc_s831683_mpnn9_S15D | protein | **5.43×** | 0.803 | 0.803 | 5 |
| 4 | bc_s831683_mpnn6_S15D | protein | **5.40×** | 0.780 | 0.751 | 5 |
| 5 | rimA02_d3_rimA_14_vhh | nanobody | 5.19× | 0.228 | 0.451 | 6 |
| 6 | rimA01_r15_boltzgen_egfr_d3_rimA_20 | protein | 4.58× | 0.577 | 0.560 | 6 |
| 7 | d2c_mpnn13_S88D_serasp | protein | 4.57× | 0.603 | 0.528 | 5 |
| 8 | bc_d3acid_l65_s831683_mpnn11 | protein | 4.01× | 0.795 | 0.784 | 6 |
| 9 | **bc_s831683_mpnn9_WT** | protein | 3.52× | 0.783 | 0.783 | 11 |
| 10 | h370_020_vhh | nanobody | 2.29× | 0.441 | 0.709 | 11 |

**The order of ranks 1–4 is not a claim.** Those four ratios span 0.22 pKa units (§3.7), inside
PROPKA's own error, so the table's ordering among them is noise. The CSV has to be submitted in
*some* order; read the top four as one result with four replicates, not as a preference.

**We gave back half the allocation on purpose.** The 20-design build existed and is in the
history. Ranks 11–20 of it did not stand on a measurement: eight read **below 1.0×** — no
switch at all — sitting on the 0.702× steric floor of §6, which is what *any* design that
touches a histidine with no carboxylate near it returns. Those rows reported the method back to
itself.

`bg04_r03` was described in an earlier version as "the weakest, 1.94× with **0.0000 on both
species**." **Both halves were wrong.** It is 1.94× over 15 poses with **human 0.0000 and mouse
0.4406** — it binds one species, not neither — and it was the **strongest** of the ten cut rows,
not the weakest. The 0.0000/0.0000 reading came from the run-name-keyed affinity path corrected
below; the "weakest" label was simply not checked. It remains out of the submission, because one
species plus three expression-QC flags including an unpaired cysteine is not a candidate, but the
reason we gave for cutting it was wrong twice over.

Also note the 0.702× floor cited above is the single-site model's own analytic minimum of 0.699×
(§1), so "sitting on the floor" means "indistinguishable from no linkage at all", which is the
stronger statement. Filling the allocation would have bought ten more wet-lab wells and no more
evidence, so we did not.

**Rank 9 is a control, deliberately.** `bc_s831683_mpnn9_WT` is the matched wild-type of rank 3,
one residue apart. It costs a slot and buys the one thing predictions cannot: if the four S15D
designs switch and this one does not, §3.7 is established in the laboratory independently of
every threshold in this document. If all five switch equally, the mutation is not the cause, and
we will know that too. That is the result we most want from this round, whichever way it goes.

**Ranks 1–4, 8 and 9 are six of ten slots on one backbone** (`d3acid_l65_s831683`) — six MPNN
sequences and one point mutant. It is a deliberate trade of diversity for replication of the
project's only causal result, made with §3.7's caveat in full view: those four ratios span
**0.22 pKa units** and are not distinguishable from one another. We are not claiming rank 1
beats rank 4. We are claiming the mutation does something, four times.

**Ranks 5 and 10 are nanobodies this instrument cannot score** (§4.5). They are included
because excluding them would be scoring them at 0.0000, which is exactly the error §4.2
documents. They also carry the upload risk described in §9.

One consistency note we are recording rather than papering over: `bin/check_discards.py` reads
the 29-design candidate JSON, not the emitted 10-row CSV, so the bar it tests against is looser
than the one that ships. It therefore flags a superset of what the shipped bar would flag — the
gate errs toward warning, which is the safe direction — and we left it alone rather than edit a
gate at submission time.

---

## 12. Limitations

1. **Every number here is a prediction of ours except the 11 in §4.4.** Nothing in this
   submission has been validated experimentally.
2. **There is no working affinity bar.** The original one derived from a control that is the
   agonist (§4.1). Its replacement, a composition-matched shuffled null, has now landed and is
   degenerate — 12 of 12 at exactly 0.0000 (§4.5) — so the `affinity_above_null` flag means
   "above zero". The only negative class with a real tail is the 10 measured non-binders, and
   that tail reaches 0.5893, above our own measured positive. Affinity is reported, not gated.
3. **The primary objective does not discriminate binders.** A measured non-binder reads 5.27×,
   8th of the 132 molecules eligible to rank (§4.4). The pH ratio ranks this submission because it is the stated
   objective, not because we have shown it selects for binding.
4. ESMFold2 does not model pH. Our ipSAE scores are pH-agnostic structural compatibility, not
   affinity at 6.5 or 7.4. The two legs are not an independent ensemble.
5. The assay uses glycosylated, tethered EGFR from HEK293. We designed against an
   unglycosylated crystal structure; domain III carries sequons.
6. **Nanobody-format designs cannot be scored on this instrument at all** — it reads a
   validated 5 nM VHH and a randomised-CDR decoy on that VHH's own framework as the same
   number, 0.0000 (§4.5). Two rows in this submission are that format.
7. The pool has been re-screened repeatedly on thresholds derived from the same measurements,
   so surviving designs are selected on the data that would evaluate them. §4.4 is the one
   exception, because its outcomes were measured by someone else before we folded anything.
8. **n = 1 on the positive class.** §4.4, the one measured-outcome test here, rests on a
   single measured binder of 53 aa; ten of the eleven published molecules are negatives. Nothing
   in that section is a discrimination estimate, and we deliberately report ranks rather than an
   ROC-AUC, which from one positive would restate one comparison ten times.
9. The two highest-ranked formats we cannot score (§4.5) are in the submission anyway, and
   two of ten rows depend on Adaptyv's ANARCI calling them antibodies (§9).
10. **The pH gate has never been validated against a measured pH outcome, and the control that
    would do it was supplied to us and not used.** G532 is a published antibody with **SPR
    ratios of 13.26× human / 3.31× mouse** at exactly our pH pair, switching on H433 and H370
    (§8.1). Our collaborator named it as the calibration pair, twice, and instructed that it be
    used for calibration only and never as a starting sequence. We did the inverse: G532 appears
    in this project as a design-generation arm and is absent from the control panel. It also
    carries published **ELISA** ratios of 8.08 / 1.64 / 0.76, which must not be compared against
    SPR KD ratios — an earlier internal analysis anchored on the ELISA figure. It was folded into 20 poses and never scored or gated; recovered, the gate
    returns 0.677 / 0.689 / 0.677 / 0.697 — wrong direction on a measured 8× switch and no
    separation from its own negative comparator. We do **not** present that as falsification,
    because the structures cannot support the mechanism being measured: the nearest carboxylate
    sits 8.97–16.84 Å from H433 across those poses, and a gate reading on such a structure is
    uninformative. What it does mean is that **the primary objective of this submission rests on
    an instrument with no validation against measured pH data** — not merely on one that "does
    not model pH" (limitation 4). Separating "the structure is wrong" from "the protonation model
    is wrong" requires scoring that series' interface, which we have not done.
11. **Effective n is 5, not 10.** Six of the ten submitted designs sit on one backbone
    (`d3acid_l65_s831683`), four of those differing only in the ProteinMPNN sequence and two by a
    single residue. Any hit rate or interval computed over designs rather than sequence families
    overstates n by up to six-fold on the arm carrying our only causal claim.
12. **The reproducibility claim has a boundary** — see the Repository note above. The pose cache
    is not published, and the documented emit command returns zeroed affinity columns without it.
