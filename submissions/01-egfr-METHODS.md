# Methods — EGFR pH-switch binder design (Challenge 1, Track 3)

**Objective.** A de novo binder to human EGFR whose affinity is higher at pH 6.5 than at
pH 7.4, cross-reactive with mouse EGFR.

**Scale.** 1,948 designs generated across 46 arms; 4,192 complexes folded and scored;
1,944 designs pH-gated on generator poses; 1,657 re-gated on independent refolds; 1,693
scored by a third, protonation-aware model. All compute on Modal (ESMFold2, BoltzGen) and
locally (PROPKA, Proton-PottsMPNN).

This document leads with what did not work and with the errors we found in our own
instrument, because those are the results we are most confident in.

---

## 1. Mechanism

We used **thermodynamic linkage**. A titratable site with pKa_free in the unbound state and
pKa_bound in the complex changes the observed affinity with pH:

    ratio = K(6.5) / K(7.4),    K(pH) = (1 + 10^(pKa_bound − pH)) / (1 + 10^(pKa_free − pH))

Two mechanisms are available:

- **Mechanism B** — a carboxylate on the *binder* raises the pKa of a histidine on the
  *target*, so the salt bridge only forms once that histidine protonates. We used this.
- **Mechanism A** — the inverse: histidine on the binder, carboxylate on the target.

**The ceiling, and why it set the strategy.** As pKa_bound → ∞ the ratio tends to
10^(7.4−6.5) = **7.943×** per site. For EGFR's H433, whose free pKa we measure at 6.22,
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

## 3. What we measured — the negative results

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

### 3.4 The two-site route is closed, and no design in the pool beats the single-site ceiling.
3 of 1,944 designs engage H433 and H370 simultaneously, all by accident. The best physically
composed multi-site product in the entire pool is **6.58×** — below the 7.943× single-site
maximum and 15% of the 43.1× two-site ceiling. Of 23 designs that reached H370 incidentally,
**9 of 10 that switched did not bind.** H370 has 153 heavy atoms within 10 Å of its ring
against H433's 57; reaching into that cleft appears to cost the interface area binding needs.

### 3.5 Binding and switching trade off, on three independent predictors.
This is the central result. Designs with the tightest carboxylate–histidine contacts switch
most and bind least. It reproduces on the generator pose, on the refold, and on a trained
protonation-aware model that shares no machinery with either.

### 3.6 Contact distance predicts switching — on the generator pose only.
Measured on all 1,944 gated poses (H433 ∪ H370), generator pose:

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

## 4. The instrument, and two errors in it

Our ranking metric is **ipSAE_min** on ESMFold2, the minimum over both alignment directions
of one interface. Two findings about its calibration matter more than any design result.

### 4.1 The negative control is EGF-derived, and its activity is unknown.
`NEG_nonbinder`, the molecule that set our "does it bind" floor, is **81% mature human EGF**:
43 of 53 residues verbatim, all six cysteines in correct linear order, two insertions in the
loops. Human EGF binds EGFR at ~2 nM and its high-affinity contacts are on domain III — the
surface every design targets. There is no provenance for this sequence anywhere in the project.

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
candidates are marked as such rather than rejected for low scores. (Two are separately
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

**93% of the switches in our pool rest on a single pose.** Of the 12 single-pose switches that
later received five poses, **5 of 12 fell below the previous threshold on re-evaluation.**
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

**Both bodies of work design binders that release in acid; neither reports a design that
binds more strongly at low pH.** Our direction is the one neither achieved — which is a
reason for humility about our predictions, not confidence about our designs.

Their histidines are on the **binder**, where they are placeable and multipliable. Ours is on
the **target**, where evolution fixed its position. That choice, made before we read either
paper, is what capped this project at a single site.

---

## 9. Limitations

1. Every number is a prediction. No experimental validation.
2. The binding bar derives from a control that is the agonist. A composition-matched
   shuffled null is the replacement; until it lands, all affinity claims are provisional.
3. ESMFold2 does not model pH. Our ipSAE scores are pH-agnostic structural compatibility, not
   affinity at 6.5 or 7.4.
4. The assay uses glycosylated, tethered EGFR from HEK293. We designed against an
   unglycosylated crystal structure; domain III carries sequons.
5. Nanobody-format designs cannot be scored on this instrument (§4.2).
6. The pool has been re-screened repeatedly on thresholds derived from the same measurements,
   so surviving designs are selected on the data that would evaluate them.
