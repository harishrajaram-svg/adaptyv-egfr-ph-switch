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

**What a ratio of 1.0 means, and what 0.699× does not.** *Corrected 2026-10-05 after
review; the earlier wording in this section was wrong.* The general one-proton upper bound over
this pH pair is **7.94×** per site. For a site with pKa_free = 6.22 the single-site model spans
**0.699× to 5.55×**, where 0.699× is the limit as pKa_bound → −∞ (proton binding abolished on
complex formation) and 5.55× the limit as pKa_bound → +∞. **No linkage gives a ratio of exactly
1.0, not 0.699×** — `link(pKa_free, pKa_free) = 1` for every pKa, which the selftest now
asserts. 0.699× is therefore not a "no-switch floor": it is the opposite extreme, maximal
*acid-weakening* linkage, a reverse switch that binds more weakly as pH falls.

This matters for how §6 reads. Our empirical floor of **0.702×** over n=570 designs sits within
0.003 of an analytic limit. We previously offered that agreement as reassurance that the method
was well-calibrated. It is better read as a warning: a large group of designs piling up on an
analytic extreme is the signature of **protonation-model saturation** — PROPKA driving
pKa_bound off-scale so the linkage term collapses to its bound — rather than evidence that
those designs genuinely occupy a physical extreme. We cannot distinguish the two with the data
we have, so the 0.702× group is reported as **uninterpretable on this gate**, not as
"no switch detected". Over a 0.9 pH-unit window two sites give 63× and three give 500×; a
single site cannot exceed 7.94× by any design.

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
not: the generator repacks the target per design, and H370's own neighbouring aspartate —
**ASP344 in PDB numbering, D368 canonical**, the convention H370 itself is quoted in — sits a
few Å from its ring. (The two numbering systems were mixed inside one sentence here, and the
2.4–3.8 Å range quoted previously does not reproduce from the apo structure; the per-design
range is what the gate actually uses, and it varies because the generator repacks the target.) A rotamer shift in ASP344 then appears as a "switch" in designs
whose binder is 27 Å away. That artifact turned 10 apparent switches into 1. Deleting the
binder from the *same file* holds every other coordinate fixed, so the difference is the
binder and nothing else.

**The pH gate runs on an independent refold, not the generator's pose.** See §5.

---

## 3. What we measured

### 3.1 Cropping the target fixes epitope selection. Three arms.
Against the full ectodomain, 4/30 designs touched the declared patch. Against a 170-residue
domain III crop, **27/30 and 29/30** did. (The header previously said "Confirmed four times"
while the body reports three arms. Three.) The generator was not ignoring the
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

### 3.4 We could not build the two-site route from scratch. One point mutation may have. See §8.1 and §11.6.
**This section has been wrong twice.** It first read "the two-site route is closed" — it is not.
It was then corrected to "we could not build it" — which was true of our *generative* search and
is no longer true of the submission: `rimA01_r15_L133E`, a single Leu→Glu at position 133 of our
rank-1 design, reads **5.659× all-site, above H433's 5.55× single-site ceiling**, with H370
contributing 0.921–2.584× while H433 holds steady (§11.6). That is one molecule with a wide pose
spread, not a demonstration — but it is no longer nothing, and this section should not be read as
saying the route failed.
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
`NEG_nonbinder`, the molecule that set our "does it bind" floor, carries **43 of the 53 residues
of mature human EGF verbatim**, all six cysteines in correct linear order. Stated precisely,
because an earlier version said "81% mature human EGF" and that is the percentage **of EGF**, not
of the molecule: `NEG_nonbinder` is **134 residues**, so 43 are EGF-derived and **91 are not** —
32% of the molecule, not 81%. "Two insertions in the loops" also undersold 81 extra residues.
What matters is unchanged and does not depend on the denominator: the EGF epitope is present and
intact. Human EGF's high-affinity contacts are on domain III — the surface every design
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

*(This project writes the full ectodomain as 609 residues in some places and 621 in others —
including §3.1 above, against this section's own argument that a bar must come from the matching
construct. **621 is correct** for what was actually folded, and it is also what the organisers
assay: their construct is Met1–Ser645, and with EGFR's 24-residue signal peptide removed that is
621 residues of mature protein. The 609 figure is a stale earlier crop and appears nowhere in the
scoring path.)*
| domain III crop (170 aa) | **0.2218** | **0.4005** |

The strict ECD bar of 0.6224 is unreachable on the crop, where a construct-matched 1 nM
binder only reaches 0.4183.

---

### 4.4 Control recovery: the instrument graded against measured outcomes.

Adaptyv published 11 molecules expressed and measured against EGFR on the assay platform this
challenge uses — human EGF, and **10 designs for which no KD was reported**
(`proteinbase.com/collections/egfr-round1-second-submission`). We folded all 11 against both
target constructs, 5 seeds each, and scored them with the same ipSAE_min that ranks our
designs. **This is the only test in the project where the instrument is graded against an
experimental outcome rather than against itself.**

**How the ten are labelled, corrected 2026-10-05.** Earlier versions of this section called
them "measured non-binders". That is not what the data says and the distinction was put to us
directly: *"Distinguish assay failure and missing KD from no binding detected. Weak binding
beyond the quantifiable limit is right-censored, not below."* A molecule with no reported KD is
one whose affinity is **right-censored** — known only to be weaker than the assay's
quantifiable limit — and that class also contains molecules that failed to express or failed QC
for reasons unrelated to affinity. It is not a measurement of zero affinity. Throughout this
section they are therefore **"no KD reported"**, and every claim below is a claim about
separating one quantified binder from ten censored observations, which is weaker than
separating binders from non-binders. We do not know the true affinity of any of the ten.

| molecule | measured | ipSAE_min human | ipSAE_min mouse |
|---|---|---|---|
| Human EGF | **KD 55 nM** (n=15) | 0.3549 | 0.3700 |
| gitter-yolo10 | no KD reported | **0.5893** | 0.0000 |
| gitter-yolo9 | no KD reported | **0.4326** | 0.0000 |
| deepsatflow-design7 | no KD reported | 0.2168 | 0.1823 |
| gitter-yolo5 | no KD reported | 0.0264 | 0.0108 |
| gitter-yolo4 | no KD reported | 0.0109 | 0.0000 |
| gitter-yolo3 | no KD reported | 0.0000 | 0.0121 |
| gitter-yolo2, 6, 7, 8 | no KD reported | 0.0000 | 0.0000 |

**What this panel can and cannot support, before the numbers.** It is **ten right-censored
observations and one quantified binder.** An AUROC *is* computable here — it is **0.80**, the
fraction of censored molecules ranking below EGF — but it is a monotone function of that single
molecule's rank and carries no information beyond it, so reporting it as a discrimination
estimate would dress one comparison up as ten. We report the rank. *An earlier version of this
paragraph said there was "no discrimination estimate to be had", which overstated it; the
statistic exists, it just adds nothing.*

**On the human leg, 8 of the 10 no-KD molecules rank below the quantified binder and 2 rank
above it.** The highest-scoring molecule in the whole measured panel is one with no reported KD:
gitter-yolo10 at 0.5893 against EGF's 0.3549. Note what this does and does not show: gitter-yolo10
may bind EGFR more weakly than the assay can quantify, so "the instrument ranks a non-binder
first" is not established — what is established is that the instrument ranks a molecule of
*unknown, weaker-than-quantifiable* affinity above a 55 nM binder.

**The cross-reactivity requirement is what rescues it.** Both molecules that outrank the binder
score exactly **0.0000 on the mouse leg**. Scored the way the submission is scored — requiring
both species — **the one quantified binder outranks all ten no-KD molecules.** Stated at the
strength the data supports: on this panel, requiring both species is sufficient to rank the only
molecule with a measured KD above every censored one. With n = 1 positive this is consistent
with the instrument working and does not demonstrate that it does, and a lack of significance
here would be inconclusive rather than a falsification.

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
is reading a number that a molecule with no reported KD scores in the top 6% of everything
eligible to rank.

Three limits. The positive class is **one molecule** of 53 aa, shorter than every design here,
so nothing in this section estimates sensitivity — the strongest honest reading is "the one
quantified binder outranks all ten no-KD molecules when both species are required." A missing
KD **bounds affinity from above, not below**: it says the interaction is weaker than the assay
can quantify, and says nothing else, so it does not prove no interaction. And ten censored
observations do not characterise a tail, so the 0.5893 we found is a floor on how high a
no-KD molecule can score here, not a ceiling.

### 4.5 The matched null is degenerate, so the screening flag does almost nothing.

12 shuffle nulls — 4 designs × 3 independent shuffles, composition and length preserved
exactly, 5 seeds nested within each shuffle — score **0.0000 on both legs with 5 of 5 dead
seeds, all twelve of them.** A separate, larger shuffle run (22 molecules, 110 poses) agrees:
**20 of 22 are dead on all 5 seeds and the highest median in it is 0.0110.** Across 34 shuffled
sequences the null does not produce a single score worth a threshold. This was predicted on
review: full-sequence shuffling destroys the fold, so the null bounds the instrument's noise
floor at exactly zero, and a "95th percentile of a matched null" flag admits anything above
0.0000. We still report `affinity_above_null` in the
submission, and it carries almost no information. **The 10 no-KD molecules are the only
comparison class with a usable tail, and that tail reaches 0.5893** — noting that they are
right-censored, so this is a tail of unknown-affinity molecules, not of non-binders.

4 framework-preserving CDR decoys — the real 5 nM VHH scaffold with its three CDRs randomised,
located by conserved anchors rather than fixed indices — score 0.0000, 0.0117, 0.0000, 0.0000.
The real 5 nM VHH scores 0.0000 as well. **The instrument cannot separate a validated nanobody
from a randomised-CDR decoy on that nanobody's own framework.** §4.2 called the nanobody blind
spot partial on a 0/2 panel; on a matched decoy class it is total.

**And on the one antibody in this project with a measured KD, the column is inverted.** G532
(§8.1) binds EGFR at 294 nM and switches 13.26×. Scored on its own ESMFold2 poses, by chain pair
because `ipsae_min` refuses on a 3-chain complex, medians over 5 seeds:

| molecule | measured | best target:Fv ipSAE | intra-Fv ipSAE |
|---|---|---|---|
| **G532** | **294 nM, 13.26×** | **0.0135** | 0.8510 |
| G532V | ELISA 1.64 | 0.3713 | 0.8634 |
| G532Ctrl | ELISA 0.76 | 0.2503 | 0.8433 |
| G5V2 | — | 0.4289 | 0.8605 |

**The real binder scores lowest of the four and its non-switching comparator scores 18× higher,
while the Fv itself folds at 0.84–0.86 in every molecule.** ESMFold2 builds the antibody and then
fails to dock it. So the affinity reading for an antibody format here is not merely uninformative
— on the single measured example available it points the **wrong way**.

This is n=1 molecule with 3 comparators from one published series, run without an MSA, so it does
not establish a general docking failure rate. It is enough to stop us reading the affinity column
for the two VHH-format rows in this submission at all, which is how they are reported.

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

**odds ratio 6.32, 95% CI [4.55, 8.79]** — Woolf's log-normal interval, which is what that
figure is and was not previously named. **Note the 2×2 totals 1,657, not the 1,693 quoted in the
sentence above it:** 36 complexes have a Potts score and no PROPKA refold verdict, so they enter
the paragraph's denominator and not the table's. Both numbers are right about different things
and the mismatch was unflagged.

**Woolf's interval assumes independent observations, and the next paragraph says they are not.**
That is the one assumption the method makes and the one this data violates — the 1,657 complexes
include repeated backbones, sequence families and poses. The interval is therefore **narrower
than the truth**, and we report it as the arithmetic of the 2×2 rather than as an inference about
designs. A family-clustered interval would be the correct object and we have not computed one.

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
| a measured pH-calibration ladder folded and never scored | the G532 series, the only molecules here with a published pH ratio, is absent from every analysis file (§13) |

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
submission and not a scoring preference. **Sourcing, stated because an earlier version asserted
the Level ≥ 3 bar as established policy without one:** the rule itself is published
(`adaptyvbio.com/blog/novelty`); the **specific gate level** comes from organiser messages in the
competition Slack channel `#anthropic_adaptyv_competition`, not from any fetchable page, and
Adaptyv confirmed on 2026-10-04 that a **self-service novelty pipeline ships 2026-10-05** so
submitters can check before uploading. A competitor independently reported on 2026-10-03 that the
live platform flags lightly modified VHH frameworks, and Adaptyv have said the antibody threshold
is being re-tested — **so the bar our two VHH rows sit near may move before this deadline.** We
treat Level ≥ 3 as the working gate and will verify against their pipeline before upload. Our implementation of the general-protein rule had one
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

**Provenance of this count, because an internal review disputed it.** The 292 rows across the
general-rule FoldSeek TSVs deduplicate to **238 distinct designs, 114 clearing**, keying on the
design basename. Stripping a trailing `_modelN` collapses two more, giving 236/113; no key gives
233/109. That figure is reproducible but stale — it is this same count computed **before
`novelty_s15d.tsv` existed**, and that file holds the five S15D designs supplying ranks 1–4 and 9
of this submission. 22 designs appear in more than one TSV, and 2 of those carry slightly
different TM or identity values between FoldSeek runs (`bg01_r02` at fid 0.127 vs 0.143;
`bg04_r03` at TM 0.7169 vs 0.7114) — **no design's level differs between rows**, so the
levelling is stable even where the underlying search is not.

Three limits we state rather than hide. Adaptyv run MMseqs2 against SwissProt, PDB, patent
sequences, the therapeutic-antibody database and PLAbDab; we search **PDB only**, via FoldSeek,
so our sequence identities are lower bounds and a design clean here may still hit a patent or a
SwissProt entry.

**And our levels are computed whole-chain, while the published rule is domain-wise and
coverage-weighted.** We implement the level thresholds correctly but not the segmentation:
Adaptyv split a protein into domains and consider structural similarity together with sequence
coverage, and **splitting can only raise TM**, never lower it. So our TM values are lower bounds
on theirs, and a design we place at Level 3 on a whole-chain TM just under 0.80 could be Level 2
domain-wise. Six of our designs are 65 aa single-domain miniproteins where a split cannot move
anything; the exposure is on the longer rows. We did not implement segmentation and the gap runs
against us, not for us. And the antibody branch needs CDRH3 identity from ANARCI numbering, which we
do not run; we implement it from conserved framework anchors instead
(`bin/antibody_novelty.py`, self-tested).

**The two nanobody rows carry upload risk, but not equally — an earlier version of this section
treated them as the same case and they are not.**

`rimA02_d3_rimA_14_vhh` is a **canonical 129 aa VHH**: it carries the standard framework and will
almost certainly be annotated as an antibody by any ANARCI-class tool. Its exposure is the
*threshold*, not the classification — CDRH3 13.6% clears comfortably, global identity 84.8%.

`h370_020_vhh` is the row actually at risk. It is **98 aa with a non-canonical FR2**, which is
precisely the shape a framework-anchored annotator can fail to call. It clears at antibody-rule
Level 4 — the only Level 4 in the project — but under the general-protein rule it is Level 2
(TM 0.914) and would be auto-rejected. **If it is not classified as an antibody, it fails.**

A second mechanism may rescue both: the published novelty page describes classification running
on **structural** similarity as well as sequence, and on that route both rows are safe. Our
earlier framing pinned their survival on ANARCI alone, which is the more pessimistic of the two
readings and was stated as though it were the only one.

---

## 10. The exclusion ledger

*Rebuilt 2026-10-05 as `bin/exclusion_ledger.py`, keyed on binder sequence. The prose ledger
this replaces is preserved in git history at commit `efe1fa4^`. No outcome has been examined.*

**What was wrong with the old ledger.** It held two lists and never reconciled them. One was
45 run names that `bin/check_discards.py` warns on — measured, outranking a submitted design,
rejected anyway — described here as "38 distinct molecules" because the gate is name-keyed on
the discard side and alias pairs double-count. The other was a set of high-ranking molecules
that were never candidates at all, reported first as 28 and then as 59 after a re-fold changed
the denominator. Adding 38 and 59 into "97 excluded" would have been wrong if the two lists
intersect, and nobody had checked. The ledger built to catch this project's keystone trap was
itself keyed on names.

**The reconciliation.**

| | |
|---|---|
| warn run names | 45 |
| of those, resolvable to a binder sequence and not themselves shipped | 23 |
| **distinct sequences (List A)** | **15** — 8 alias collapses |
| **sequences outranking the weakest finalist, n ≥ 5, never candidates (List B)** | **60** |
| **overlap** | **0 — the lists are disjoint** |
| **union: distinct excluded molecules** | **75** |

Two of the 45 warn names resolved to `rimA01_r15`, which is **shipped at rank 1** — including
its own exact name. That is the same error §10 previously caught once by hand (citing 5.630×,
which is `bc_s360518_mpnn9_A22D` at rank 2) now found systematically. The gate itself was also
comparing against the wrong set: it read `submission_final.json`, a **31-design candidate
pool**, not the 12 designs shipped, so its threshold came from a design we did not submit and
the 19 unshipped candidates were skipped by a name-stub test and never checked against the
shipped set. Fixed; it now reads the graded CSV and matches by sequence, and reports 23 warns
against a 2.289× bar rather than 45 against 1.263×.

**Classification.** Every molecule in the union falls in exactly one category, in precedence
order. These are not equivalent and were previously pooled.

| category | n | meaning |
|---|---|---|
| **ELIGIBILITY** | 4 | fails a rule that disqualifies it whatever it scored — all four fail the novelty gate. Not reopenable. |
| **ASSESSMENT** | 8 | never adequately measured: all eight have 1–2 pH poses against a 5-pose minimum. Excluded for want of evidence, not on evidence. |
| **GATE-ONLY** | 63 | eligible and adequately measured, excluded only on a pH threshold or a judgement call. |

The 63 gate-only exclusions are **reopened**. The reason is the reviewer's, and we accept it:
pH-gate-dependent rejections are *unsupported by a validated selection rule*, so a gate-only
exclusion is not a finding. Note that the highest-ratio molecule in the whole warn list,
`domIII_1His_ctrl_r17` at 6.60×, is an ASSESSMENT case with a single pose — it was never a
measurement.

**One assessment across both pools, which closes the gap this section used to admit.** The
previous version of §10 said: *"all 59 are ranked on the superseded target-only basis — the
all-site product was computed only for the twelve submitted designs, so we cannot say how they
would rank on the basis this submission actually uses. That is the more serious half of the
gap."* It is now closed. All 63 reopened molecules were re-scored over their own human-leg
poses with the **same two-partner histidine-only gate that ranks the finalists**, same code
path, 402 poses. Results in `analysis/01-egfr/exclusion_ledger.json`.

**The result does not favour the submission.** Applying the submission's own tier-1 rule —
ratio ≥ 1.20, n ≥ 5 poses, pose spread ≤ 1.0× the median — and the submission's own ranking
basis:

- **36 reopened molecules clear tier 1 and outrank the weakest shipped tier-1 design**
  (`bc_s831683_mpnn19_S15D`, 1.774×). One is an Adaptyv control molecule and not a candidate;
  **35 are our own designs.**
- Several beat shipped designs on both axes at once.
  `ss_bc_s831683_mpnn6_S15D_S62H_routeA` reads **3.545×** with human 0.765 / mouse 0.744,
  against shipped `bc_s831683_mpnn6_S15D` at **1.835×** with human 0.780 / mouse 0.751 — the
  same backbone with one further mutation, essentially the same predicted affinity, and
  **1.9× the pH ratio**. `sd_d2c_101_l147_s144898_m_T65D` reads 4.735× with human 0.585.
- The highest, `c5_cf_short__boltzgen_egfr_cropfree_short_48` at **5.546×**, would rank second
  of everything on this basis, but scores human 0.241 / mouse 0.181.

Read plainly: **the stated primary objective prefers designs we did not submit.** The shipped
set was assembled with affinity carrying more weight than the objective licenses — this
document already states that "no design is excluded on affinity", and the ranking nonetheless
leaned on it. That is a selection effect in the submission, recorded here rather than left for
a grader to find.

**What blocks acting on it, and it is a real block.** **None of the 35 has a novelty
assessment.** Novelty Level ≥ 3 is an Adaptyv eligibility requirement, not a preference, and
every one of the 35 carries `NOT ASSESSED (no novelty record)`. Expression QC is likewise
unavailable: it was only ever run on the 20 submission candidates, so for the reconsidered pool
that axis is unassessed too, and the ledger records it as unassessed rather than defaulting it
to eligible. Novelty is also a severe filter on this target — §9 of this document measures
Adaptyv's own round-2 set at 2% of *binders* clearing the strict reading — so a high ratio is
no guarantee any of these is submittable.

**The novelty gate has now run over those 35, and it does not rescue the submission.**
FoldSeek against the PDB database, same `bin/novelty_gate.py` and same Level ≥ 3 bar applied to
the finalists (`analysis/01-egfr/novelty_reopened.tsv`):

| novelty level | n |
|---|---|
| Level 3 — **clears the eligibility gate** | **25** |
| Level 2 | 6 |
| Level 1 | 4 |
| Level 4 (fully de novo) | 0 — Adaptyv's round-2 baseline was 2.8% |

**25 of the 35 are eligible.** The blocker named above is removed for most of them, and the
submission is carrying **12 of the 20 permitted designs**, so eight slots are unused. Three of
the ten failures are flagged ANTIBODY and were scored by the general-protein rule, which §9
notes is *stricter* than the antibody rule — including `ss_rimA02_d3_rimA_14_vhh_T28H_routeA`
at 3.209×, a variant of shipped rank 6 — so those three need an ANARCI re-check before being
treated as excluded.

The cases that most directly contradict the current ranking, all Level 3:

| molecule | pH (his-only) | human | mouse | spread | compare |
|---|---|---|---|---|---|
| `ss_bc_s831683_mpnn6_S15D_S62H_routeA` | **3.545×** | 0.765 | 0.744 | 0.358 | shipped `mpnn6_S15D`: 1.835×, 0.780, 0.751 |
| `sd_d2c_101_l147_s144898_m_T65D` | **4.735×** | 0.585 | 0.366 | 0.944 | shipped `d2c_mpnn13_S88D`: 3.526×, 0.603, 0.528 |
| `ss_bc_s831683_mpnn19_S15D_S62H_routeA` | **3.189×** | 0.642 | 0.601 | 0.591 | shipped `mpnn19_S15D`: 1.774×, 0.808, 0.786 |
| `bcr_d3acid_l65_s831683_mpnn3_S15D` | 2.035× | 0.759 | 0.754 | 0.512 | beats shipped `mpnn6_S15D` at matched affinity |
| `bcr_d3acid_l65_s831683_mpnn17_S15D` | 1.814× | 0.780 | 0.755 | **0.019** | the most reproducible pH measurement in the project |
| `c5_cf_short__boltzgen_egfr_cropfree_short_48` | **5.546×** | 0.241 | 0.181 | 0.262 | would rank 2nd of everything; weak on both species |

`ss_bc_s831683_mpnn6_S15D_S62H_routeA` is the clearest: it is a shipped design plus one further
mutation, with predicted affinity indistinguishable from it (0.765/0.744 against 0.780/0.751)
and **1.9× the pH ratio**, at an acceptable spread and Level 3. On this submission's own
criteria it dominates a design we shipped.

Note also the family composition: **15 of the 25 are the `d3acid3_l60_s647537` backbone**, a
family with no representation in the submission at all, while six of the twelve shipped designs
sit on a single other backbone.

**What remains unresolved, stated as such.** This section does not change the submission, and
the decision whether to use the eight unused slots is not made here. Three things are known and
recorded: the exclusion of 63 molecules rested on a gate that cannot carry that weight; 25 of
them are eligible and outrank a shipped tier-1 design on the submission's own basis; and all of
these comparisons inherit the provisional status of that basis, since §11.7 measures Kendall
τ = +0.000 between it and the partnered alternative. Expression QC remains unrun for every one
of the 25. The original exclusion list and this amendment are both preserved and timestamped,
and every number here was produced before any outcome was examined.


## 11. The submission

**12 designs, ranked on the all-titratable-site pH product.** Track 3 allows 20.

*The tables and counts in this section are GENERATED from the emitted CSV by `bin/gen_methods_submission.py` (self-tested). They were hand-maintained through three submission changes in one evening and drifted badly — an audit found this section still describing eleven designs, its rank table omitting the twelfth, every rank above 7 off by one, and §12 attesting review of "all ten" sequences. The interpretive text is still written by hand; the numbers are not.*

### 11.1 The ranking basis changed, and it reordered everything

Until 2026-10-04 we estimated the pH ratio with a gate that measures only the **target's**
histidines. It never measured our own binders' titratable groups — and **eight of the twelve
submitted designs carry between one and three histidines of their own**. (An earlier version of
this sentence said "six of the eleven", which was wrong on both numbers; an audit then proposed
"seven of twelve", which was also wrong. Counted from `multisite_pooled.json`: ranks 2, 3, 4, 5,
9, 10, 11 and 12 carry at least one, and ranks 1, 6, 7 and 8 carry none.) Those get buried at the interface
and lose 1.5–2.5 pKa units, and by the same thermodynamic linkage of §1 that **opposes**
acid-tightening. We were counting the target's sites and ignoring ours.

`bin/ph_gate_multisite.py` composes over the **histidines** of both partners, with the free
leg taken by deleting the other chain in place — the mirror of the argument §2 makes for the
target leg. *Corrected 2026-10-05: an earlier version of this sentence said "every titratable
site", which the code did not do — it parsed HIS/ASP/GLU but added only histidines to the
product. The shipped column is a **two-partner histidine-only approximation**. §11.7 gives the
all-site and partnered alternatives, the resulting ranking instability, and why the
histidine-only order is nevertheless the one retained.* Measured over 75 poses, n = 5–11 per
design:

<!-- GENERATED:BASIS-TABLE -- do not edit by hand; `bin/gen_methods_submission.py --write` -->
| design | target-only | **all-site** | binder histidines | worst drag |
|---|---|---|---|---|
| rimA01_r15_L133E | 4.582 | **5.656** | 0 | — |
| rimA02_d3_rimA_14_vhh | 5.186 | **4.838** | 0 | — |
| rimA01_r15_boltzgen_egfr_d3_rimA_20 | 4.582 | **4.256** | 0 | — |
| bc_s360518_mpnn9_A22D | 5.630 | **3.738** | 1 | 0.776 |
| d2c_mpnn13_S88D_serasp | 4.572 | **3.526** | 2 | 0.983 |
| h370_020_vhh | 2.289 | **2.101** | 0 | — |
| bc_s831683_mpnn6_S15D | 5.397 | **1.835** | 3 | 0.661 |
| bc_s831683_mpnn19_S15D | 5.435 | **1.774** | 3 | 0.66 |
| bc_s831683_mpnn9_S15D | 5.428 | **1.062** | 3 | 0.344 |
| bc_s831683_mpnn8_S15D | 5.461 | **1.023** | 3 | 0.333 |
| bc_d3acid_l65_s831683_mpnn11 | 4.010 | **0.737** | 3 | 0.359 |
| bc_s831683_mpnn9_WT | 3.522 | **0.627** | 3 | 0.338 |
<!-- /GENERATED:BASIS-TABLE -->

**Every binder histidine moves down — 0.33 to 0.98, none up.** PROPKA noise would scatter both
ways. The gate's counter-charge guard fires on nearly all of them (nearest opposite charge
6.9–8.8 Å), so these are desolvation shifts with no electrostatic partner: the same mechanism as
the 0.702× steric floor of §6, and the same physics that defeated mechanism A (§3.5).

**This is not a different objective. It is a less wrong estimate of the same one** —
KD(7.4)/KD(6.5). The superseded number ships as its own CSV column so the change is auditable
rather than silent.

**Four of the twelve are no longer switches.** The bottom four rows fall below the 1.20× bar that §6 sets
as PROPKA's noise floor, so they carry no pH claim and are ordered by mouse affinity instead.
`bc_s831683_mpnn8_S15D` led this submission at 5.461× before the correction and is now rank 11
at 1.023×.

**What we are NOT claiming.** The binder's free leg comes from deleting the target in place, so
the isolated binder is not relaxed — a histidine buried in the complex may be solvent-exposed in
the real free binder, which would make its free pKa wrong. PROPKA on a buried histidine is its
hardest case. The *direction* is consistent across 76 poses and mechanistically coherent; the
*magnitude* is not established.

### 11.2 Final ranks

<!-- GENERATED:RANK-TABLE -- do not edit by hand; `bin/gen_methods_submission.py --write` -->
| rank | design | class | family | aa | **pH his-only (ranked)** | all-site | partnered | rank range | pose spread | target-only | poses | human | mouse | affinity assessable |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `rimA01_r15_boltzgen_egfr_d3_rimA_20` | protein | rimA01_r15_d3_rimA_20 | 150 | **4.256** | 34.534 | 4.843 | 3-11 | 0.12 | 4.582 | 6 | 0.594 | 0.567 | yes |
| 2 | `bc_s360518_mpnn9_A22D` | protein | d3acid3_l65_s360518 | 65 | **3.738** | 88.593 | 33.416 | 1-4 | 0.45 | 5.630 | 5 | 0.451 | 0.474 | yes |
| 3 | `d2c_mpnn13_S88D_serasp` | protein | d2c_101_l147_s144898 | 147 | **3.526** | 6.681 | 6.663 | 4-8 | 0.48 | 4.572 | 5 | 0.603 | 0.528 | yes |
| 4 | `bc_s831683_mpnn6_S15D` | protein | d3acid_l65_s831683 | 65 | **1.835** | 5.949 | 12.995 | 2-7 | 0.04 | 5.397 | 5 | 0.780 | 0.751 | yes |
| 5 | `bc_s831683_mpnn19_S15D` | protein | d3acid_l65_s831683 | 65 | **1.774** | 5.486 | 11.930 | 3-8 | 0.50 | 5.435 | 5 | 0.808 | 0.786 | yes |
| 6 | `rimA02_d3_rimA_14_vhh` | nanobody | rimA02_d3_rimA_14 (VHH) | 129 | **4.838** | 4.976 | 5.183 | 2-10 | 0.30 | 5.186 | 6 | 0.219 | 0.447 | **no** |
| 7 | `h370_020_vhh` | nanobody | h370_020 (VHH) | 98 | **2.101** | 2.140 | 2.267 | 6-12 | 0.59 | 2.289 | 11 | 0.417 | 0.709 | **no** |
| 8 | `rimA01_r15_L133E` | protein | rimA01_r15_d3_rimA_20 | 150 | **5.656** | 52.181 | 7.288 | 1-5 | 4.38 | 4.582 | 20 | 0.598 | 0.434 | yes |
| 9 | `bc_s831683_mpnn9_S15D` | protein | d3acid_l65_s831683 | 65 | **1.062** | 1.804 | 6.801 | 6-10 | 1.21 | 5.428 | 20 | 0.804 | 0.804 | yes |
| 10 | `bc_s831683_mpnn9_WT` | protein | d3acid_l65_s831683 | 65 | **0.627** | 1.195 | 5.226 | 9-12 | 1.76 | 3.522 | 26 | 0.786 | 0.784 | yes |
| 11 | `bc_d3acid_l65_s831683_mpnn11` | protein | d3acid_l65_s831683 | 65 | **0.737** | 1.841 | 7.685 | 4-11 | 0.64 | 4.010 | 6 | 0.796 | 0.784 | yes |
| 12 | `bc_s831683_mpnn8_S15D` | protein | d3acid_l65_s831683 | 65 | **1.023** | 1.644 | 6.789 | 7-11 | 0.16 | 5.461 | 5 | 0.776 | 0.764 | yes |
<!-- /GENERATED:RANK-TABLE -->

**Assessable designs rank ahead of unassessable ones within tier 1.** On a pure pH ordering
`rimA02_d3_rimA_14_vhh` leads the submission at 4.838× — on a human ipSAE of 0.219 that we
cannot interpret, because §4.5 shows this instrument scores a measured 294 nM antibody **below
its own non-switching comparator**. We neither demote it on the pH axis nor score it at 0.0000
(the §4.2 error); we place it after the designs where both axes mean something.

**The cost, stated:** rimA02 carries the second-highest honest pH ratio in the submission and
sits at rank 6, below a design at 1.774×. If the organisers rank strictly on the primary
objective, this ordering costs us. It is a judgement that credible-interface-first is the more
defensible frame, following the reviewer instruction to apply eligibility and interface checks
before the challenge priorities — not a claim that rimA02 is worse.

### 11.3 Six families, twelve designs

`d3acid_l65_s831683` **×6** (ranks 4, 5, 9, 10, 11, 12) · `rimA01_r15_d3_rimA_20` **×2**
(ranks 1 and 8 — rank 8 is `L133E`, a single-residue mutant of rank 1, so it is the **same
family**, not a new one) · `d3acid3_l65_s360518` (A22D) · `d2c_101_l147_s144898` ·
`rimA02_d3_rimA_14` (VHH) · `h370_020` (VHH).

**Effective n is 6 clusters, not 12 designs.** The 6-member cluster occupies ranks 4, 5, 9, 10,
11 and 12 — the basis correction moved it out of the top three entirely. Any interval must be
computed on families, not designs. **Two pairs differ by a single residue**:
`bc_s831683_mpnn9_S15D` and `bc_s831683_mpnn9_WT` (Ser15 vs Asp15), and
`rimA01_r15_boltzgen_egfr_d3_rimA_20` and `rimA01_r15_L133E` (Leu133 vs Glu133). In both cases
the pair is deliberate — the unmutated parent is submitted alongside its mutant so the mutation's
effect is measured in the laboratory rather than inferred from our gate — but neither pair is two
independent tests.

### 11.4 Reading conventions

* **all-site pH ratio** — median over every ESMFold2 refold pose of that exact binder sequence,
  composed over all titratable sites on both partners. A prediction, not a measurement.
* **ipSAE human / mouse** — median of 5 seeds, pooled by sequence from `master_rank.json`.
  **pH-agnostic.** For ranks 6 and 7 not interpretable at all.
* **`bc_s831683_mpnn9_WT` is a control, not a candidate**: the matched wild-type of
  `bc_s831683_mpnn9_S15D`, one residue apart. Its all-site product of 0.593× is itself
  informative — the unmutated backbone is predicted to bind *worse* in acid, which is the
  baseline the Ser→Asp install has to beat. The same is true of
  `rimA01_r15_boltzgen_egfr_d3_rimA_20` against `rimA01_r15_L133E`.
  *(Designs are referred to here by NAME rather than by rank. Rank references inside a
  rank-ordered file drift every time the file changes, and an audit found six of the twelve
  assessment strings citing the wrong design by rank for exactly that reason.)*
* `affinity_above_null` means `max(hu, mo) ≥ 0.2218`, a retained legacy threshold. See §4.5.

### 11.5 Why not 20

The 20-design build is in the history. Ranks 11–20 of it did not stand on a measurement: eight
read below 1.0× even on the old, more generous basis. Filling the allocation would have bought
ten more wet-lab wells and no more evidence. A22D was added on 2026-10-04 because it is a
backbone not otherwise represented, carries the largest causal swing in the set (0.826 → 5.630
target-only; 0.574 → 3.738 all-site, so the swing survives the correction), and has the tightest
seed reproducibility in the project.

One consistency note recorded rather than papered over: `bin/check_discards.py` reads the
30-design candidate JSON, not the emitted CSV, and still compares on the target-only ratio. It
therefore flags a superset of what the shipped bar would flag — the safe direction — and we left
it alone rather than edit a gate at submission time.

### 11.6 `rimA01_r15_L133E` — a specific mechanism with an unreproducible magnitude

*Rewritten 2026-10-05 after a 15-seed-per-variant L133 triad landed (`runs/esmfold2/w3_triad`).
Three numbers in the previous version of this section were measured on five poses and did not
survive deeper sampling. They are corrected below and the old values named.*

A single Leu→Glu at position 133 of rank 1, which is submitted unmodified alongside it as the
matched parent. We folded **Gln** and **Asp** at the same position as controls, 15 seeds each,
human and mouse legs. The triad's L133E binder sequence is byte-identical to the submitted one,
so its poses are pooled with the original five: **n = 20**, not 5.

| at position 133 | side chain | residue 133's own ratio | nearest counter-charge | H433 | H370 | his-only pH | human | mouse |
|---|---|---|---|---|---|---|---|---|
| — (parent) | Leu | not titratable | — | — | 0.92–0.94× | 4.256× | 0.594 | 0.567 |
| **Gln** | neutral, Glu-sized | not titratable | — | 4.668× | **0.926×** (0.912–0.968) | **4.229×** | 0.714 | 0.720 |
| **Glu** (submitted) | acid, reaches | **1.141×** | **4.26 Å** (2.98–5.35) | 4.633× | **1.327×** (0.921–**6.229**) | **5.656×** | 0.598 | 0.434 |
| **Asp** | acid, one CH₂ shorter | 1.005× | **7.33 Å** (5.15–25.91) | 2.267× | 0.955× (0.890–1.121) | 2.221× | **0.000** | 0.420 |

*Per-site and H433/H370 figures are from the 15 triad poses of each variant; the his-only
column for L133E is the pooled n = 20.*

**The mechanism is specific, on two independent controls.** Gln is the same size as Glu and
carries no titratable group: it leaves H370 at **0.926×**, exactly where the parent and every
other design in this submission sit, and leaves the overall ratio at 4.229× against the parent's
4.256×. So the gate is not responding to "a substitution happened at 133" — the position
tolerates substitution with no pH effect at all. Asp is an acid but one methylene shorter: its
nearest counter-charge sits at **7.33 Å** where Glu's sits at **4.26 Å**, it contributes
**1.005×** — nothing — and the overall ratio *falls* to 2.221×. Only the variant whose side
chain can physically reach H370 moves H370, and it is the only one of the three that does.
Length specificity across a one-methylene difference is the strongest mechanistic evidence in
this submission.

**It is still a two-site reading.** Over the 15 triad poses, **9 of 15 exceed the 5.55×
single-site ceiling** for H433 alone (§1) and **3 of 15 exceed the 7.94× one-proton bound**.
With two sites moving the ceiling is 7.94² = 63×, so these are not ceiling violations — they are
values one site cannot produce. H433 holds steady at 4.633× while H370 does the swinging.

**What did not survive the deeper sampling.** Three claims in the previous version are
withdrawn or corrected:

1. *"affinity held: human went **up**, 0.594 → 0.616."* At n = 20 the human leg is **0.598**
   against the parent's 0.594 — **flat, not up**. The apparent gain was five-pose noise. Note
   also that Gln reads 0.714, *better* than both, so position 133 is not where this binder's
   predicted affinity is limited.
2. *`L133D` his-only "0.723×".* At n = 15 it is **2.221×**. The earlier figure came from too
   few poses to estimate; L133D's own spread is 1.78.
3. *The Glu/Asp reach argument was geometric inference* — "Glu reaches ~3.9 Å from CB,
   aspartate ~2.5 Å", from LEU133's 6.28 Å to H370's ring nitrogen. It is now **measured** per
   pose: 4.26 Å median for Glu, 7.33 Å for Asp. The conclusion is unchanged and the basis is
   better.

**Why it ranks 8th and not 1st, now more firmly.** Its pose spread was reported as **1.31×** its
median on five poses. On twenty it is **4.38×**, range **3.838 to 28.629**. Quadrupling the
sampling left the median almost exactly where it was (5.659 → 5.656) and more than tripled the
measured scatter — the central estimate is reproducible, the magnitude is not. What varies pose
to pose is whether the glutamate reaches H370 at all, and the effect over its parent (1.33×)
sits well inside that scatter. `SPREAD_BAR = 1.0` therefore keeps it below every reproducible
tier-1 design, which is where it belongs: a real and specifically-controlled mechanism whose
size we cannot quote.

**What this is not.** The ipSAE differences across the triad are changes in a *confidence*
score, not measured retained affinity, and L133D's human 0.000 means the predictor found no
confident interface — not that the design was measured not to bind. The pH numbers remain
PROPKA estimates on fixed ESMFold2 conformations with the free leg taken by partner deletion.
Per §11.7 this design's all-site value (52.2×) is dominated by `binder:ASP33` at 7.37× with its
nearest counter-charge **9.84 Å** away, a desolvation artefact present in the parent too; the
designed GLU133 contributes 1.141×. The claim here rests on the Gln/Asp controls, not on the
size of any single number.

### 11.7 pH sensitivity analysis: the ranking basis is not established

*Added 2026-10-05 after review. This section exists because the basis we ranked on is not
the one §11.1 said it was.*

**What the shipped column actually measures.** `bin/ph_gate_multisite.py` was described
here and in the preregistration as composing "every titratable site on both partners". It
does not. It computes PROPKA pKa values for HIS, ASP and GLU in both legs and then adds
**only histidines** to the product; the acids are parsed and discarded. The shipped
`ph_ratio_*` column is therefore a **two-partner histidine-only approximation**, and the
column is now named that way. This gap mattered more than a missing term usually would,
because the designed intervention in most submitted families *is* an acid — A22D, S88D,
S15D, L133E, T65D, S60D. The gate was blind to the residue each design was built around,
and Mechanism B (a binder carboxylate reading a target histidine) is invisible to a
histidine-only product by construction.

**Three bases, same poses, same code path.** `bin/ph_sensitivity_multisite.py` recomputes
every submitted design on the same 75 human-leg poses through one code path, so the
differences below are attributable to the composition rule alone:

- **his-only** — histidines on both partners. The shipped basis.
- **all-site** — every HIS/ASP/GLU on both partners, which is what §11.1 claimed.
- **partnered** — every site whose nearest opposite charge on the other chain is within
  6 Å. `PARTNER_CUT` was fixed before this analysis, not tuned to it.

The histidine-only values reproduce the shipped CSV exactly on all twelve designs, which
confirms the join and the characterisation above.

**The ordering is not stable.** Kendall τ between the shipped basis and the partnered basis
is **+0.000** — the two orderings are uncorrelated. Designs move by up to **8 ranks**
(`rimA02_d3_rimA_14_vhh`: 2nd on his-only, 10th on partnered). Per-design rank ranges are in
the §11 table. **Every tier in this submission is therefore marked `provisional`, and no
order here should be read as established.**

**Why the shipped order is nevertheless retained.** The histidine-only value is the
**minimum of the three bases for all twelve designs**. Ranking on it is ranking on the
conservative envelope `min(his-only, all-site, partnered)` under a single uniform rule,
rather than on a basis chosen after seeing which order it produced. We did not revert to
the target-only ratio, and we did not promote either wider basis.

**Why neither wider basis can rank.** The only matched negative control we have —
`bc_s831683_mpnn9_WT`, the parent of `mpnn9_S15D`, carrying no designed acid — reads
**5.344× on the partnered basis**, above two shipped designs. A basis on which the
do-nothing parent looks like a 5× switch does not discriminate designs from their parents.
The all-site basis has a related defect: the largest single contributor to
`rimA01_r15_L133E`'s 53× is `binder:ASP33` at ratio 7.37 with its nearest counter-charge
**9.84 Å** away — a desolvation shift with no electrostatic partner, the same artefact class
this project used to rule out Mechanism A — while the actual designed `GLU133` contributes
only **1.30×**. This is independent support for the reviewer's point that L133E's
5.659 vs 5.554 is not evidence of a second site.

**Noise does not accumulate.** With 30–123 titratable sites per complex, a multiplicative
product invites the objection that it compounds PROPKA noise. Measured: the product over
sites the gate itself calls unmoved (|ratio−1| < 0.05) sits at **0.946–1.020** across all
twelve. Sub-threshold noise cancels; the products are driven by genuinely shifted sites.

**The one design with a mechanism you can point at.** On the all-site and partnered bases
`bc_s360518_mpnn9_A22D` ranks 1st, and it is the only submitted design whose two dominant
sites are mutually partnered: `binder:ASP22` (pKa 6.05 → 8.76, ratio 5.90) and
`target:H433` (6.19 → 9.61, ratio 5.63), **2.78 Å apart**, both shifting in the
switch-favouring direction with implied pKa_bound values inside the plausible window. That
is the designed Mechanism B, and the histidine-only gate scored it 3.738× because ASP22 —
the designed residue — never entered the product. It ranks 4th here only because we rank on
the conservative envelope. We flag it as the design most worth wet-lab attention, and we are
not reordering the submission on the strength of that judgement.

**Unassessed sites.** The gate now records, per pose, any titratable site for which a pKa
is unavailable in either leg, with the reason, instead of skipping it with a bare
`continue`. Across all 75 poses of the twelve submitted designs the count is **0** — every
site entered or was accounted for. Site-level detail for every pose is retained in
`analysis/01-egfr/ph_sensitivity.json` so no later question requires a re-run.

**What this is not.** These three numbers are a sensitivity analysis, not a confidence
interval. They bound how much the composition rule moves the answer; they say nothing about
whether PROPKA's pKa values are right, and the partner-deletion free leg remains a
fixed-conformation diagnostic rather than a measurement of the apo state.

## 12. Declarations

Stated because the organisers ask for them and an earlier version of this document made none.

**AI assistance.** This submission was produced by one person working with Claude (Anthropic)
throughout: design generation, scoring, analysis code, and the drafting of this document. Every
number here was computed by code in the published repository, and the code was written in that
collaboration. The errors in §7 were found the same way.

**Human review.** The submitting researcher has reviewed all **twelve** submitted sequences, their
molecule_class labels, their lengths, and the claims made about them in this document and in the
CSV. The sequences are de novo designs from this project's own generation runs; none is a
modification of a previously submitted design or of an existing characterised binder. The two
published binders that appear anywhere in this project — G532 and cetuximab-derived controls —
were used for calibration only and are not ancestors of any submitted sequence.

**Tools and licences.** BoltzGen and BindCraft for generation; ESMFold2 (via Modal) for structure
prediction; the Dunbrack `ipsae.py` v4 reference, MIT, commit `6174cf9e` for interface scoring;
PROPKA 3.5.1 for pKa; Proton-PottsMPNN (`potts_v6_afdb_edge_his0.3_acid0.06`) for the protonation
model; FoldSeek against PDB for novelty; HBPLUS v3.06. All open-source and used within their
licences. The vendored reference implementation in `outbox/ipsae-fixtures/vendor/` carries its
MIT licence file.

**Structures.** Predicted complexes are published at `submissions/structures/`, one
median-ipSAE pose each — not the best pose, which would be selection on the outcome.
**Coverage is 10 of 12**: the two designs added latest, `bc_s360518_mpnn9_A22D` and
`rimA01_r15_L133E`, have no published structure yet. Their poses exist and are scored; they are
simply not exported. Stated rather than implied, because an earlier version of this line claimed
full coverage of "all ten designs" when the submission held twelve.

**Funding and compute.** Self-funded. $503 of personal Modal spend on this target. No
institutional affiliation, no grant, no commercial interest in the outcome.

**Prior work by others that this rests on.** The G532 result (§8.1) and the two histidine sites
it implicates were supplied by a collaborator reviewing this work, not discovered here.

---

## 13. Limitations

1. **Every number here is a prediction of ours except the 11 in §4.4.** Nothing in this
   submission has been validated experimentally.
2. **There is no working affinity bar.** The original one derived from a control that is the
   agonist (§4.1). Its replacement, a composition-matched shuffled null, has now landed and is
   degenerate — 12 of 12 at exactly 0.0000 (§4.5) — so the `affinity_above_null` flag means
   "above zero". The only comparison class with a real tail is the 10 right-censored no-KD
   molecules, and that tail reaches 0.5893, above our own quantified positive. Affinity is
   reported, not gated.
3. **The primary objective does not discriminate binders.** A molecule with no reported KD
   reads 5.27×,
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
11. **Effective n is 6, not 12.** Six of the twelve submitted designs sit on one backbone
    (`d3acid_l65_s831683`) and a further two are a parent/point-mutant pair
    (`rimA01_r15_d3_rimA_20` and its `L133E`). Any hit rate or interval computed over designs
    rather than sequence families overstates n by up to six-fold on the arm carrying our only
    causal claim. See §11.3 for the partition.
12. **The reproducibility claim has a boundary** — see the Repository note above. The pose cache
    is not published, and the documented emit command returns zeroed affinity columns without it.
13. **The design family with the best measured prior is the one we scored least.** 60 of the 71
    recovered BindCraft sequences sit at **n = 1 pose**, below the n ≥ 5 floor §6 requires, so
    they are structurally ineligible for tier 1 regardless of merit. 22 of those read ≥ 2.0× at
    n = 1. BindCraft also beat BoltzGen on every axis we measured, on 11 invocations against
    1,944 designs. We did not resolve them, and the submission is poorer for it.
14. **The organisers rank outcomes partly on affinity at pH 6.5, and two of our ten rows have no
    usable affinity reading at all** (§4.5). We submitted them anyway, because excluding them
    would mean scoring them at 0.0000, which is the error §4.2 documents — but it means a fifth
    of the submission cannot compete on one of the stated criteria.
