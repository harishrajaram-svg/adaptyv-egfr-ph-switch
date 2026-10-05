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

**What is and is not reproducible from it, measured rather than asserted.** The analysis code,
the submission, every artifact under `analysis/`, and the full ranking table are published. The
**6,011 cached pose outputs under `runs/` are not** — they are gitignored and run to tens of
gigabytes.

**MEASURED 2026-10-05, in a fresh `git clone` of the public repository with no local state:**

| | |
|---|---|
| `python3 bin/emit_submission_csv.py` | reproduces `submissions/01-egfr.csv` **byte-identically** |
| `python3 bin/gate_sweep.py` | **12 of 13 gates pass**, including `run_fixtures --check` at 12/12 and all ten fail-closed regressions |
| the one failure | `check_discards`, which scans for pose files to confirm a design was measured and cannot do that without the cache |

**This paragraph used to say the opposite, and the correction is the point.** Until 2026-10-05
it read that running the documented emit command on a fresh clone produced "a complete,
plausible CSV in which every affinity column reads 0.0000", with exit code 0 and no warning —
this project's own signature failure, a `0.0000` meaning *absent* read as *measured*,
reproduced inside its reproducibility claim. That was true when written and had never been
retested. It is now false, because the claim was tested instead of repeated.

Getting there needed two fixes, both found by actually cloning rather than by reasoning about
it. `bin/ipsae_min.py` hardcoded `ROOT/.venv/bin/python`, so the production scorer died with
`FileNotFoundError` in any checkout; it now falls back to the running interpreter. And it looked
for the reference only at `ROOT/ipsae/ipsae.py`, which is gitignored and therefore **not
published** — so the scorer failed closed on all 12 fixture cases for every reader, while
passing locally. It now falls back to the pinned, sha256-guarded copy vendored with the
fixtures, which is the copy that *is* published. That is the same fault the external reviewer
reported as a 404, one level deeper: the file had been published and the code still could not
find it.

**So the boundary is now one gate, not the submission.** The graded CSV, the methods document's
generated blocks, the scorer, the fixtures and the fail-closed regressions all reproduce from
the clone alone. The pose cache is still needed to re-derive the affinity columns from
structures and to run `check_discards`; ask and we will supply it.

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
switch at n ≥ 5, the engaged site is **H433 in 48, H370 in 1 and H383 in 1** — zero tag.
<!-- GENERATED:SWITCH-SITE do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
15 of the 16 submitted designs switch on **H433**. The remaining 1: **H370** -- `c5_cf_short__boltzgen_egfr_cropfree_short_48` (rank 1).
<!-- /GENERATED:SWITCH-SITE --> We designed against a tag-free crystal structure, so the tag
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

> **NUMBERING CONVENTION — READ THIS BEFORE ANY RESIDUE NUMBER IN THIS DOCUMENT.**
> This document uses **two** conventions and did not say so until 2026-10-05.
>
> | | histidines | glycan sequons | domain boundaries | published structures |
> |---|---|---|---|---|
> | convention used | **UniProt P00533 canonical** (precursor) | **mature** | **mature** | **mature** |
>
> EGFR's signal peptide is 24 residues, so **mature = canonical − 24**. The five
> domain-III histidines are therefore:
>
> | this document calls it | the structures call it (mature) |
> |---|---|
> | **H433** (the switch site) | **H409** |
> | **H370** (the second site) | **H346** |
> | H383 | H359 |
> | H418 | H394 |
> | H358 | H334 |
>
> Every `.cif` in `submissions/structures/`, `targets/egfr/egfr_d3_6aru.pdb` (numbered
> 311–480), and the contact footprints in `analysis/01-egfr/finalist_footprints.json` use
> the **mature** numbers. So the footprint of every design contains **409**, not 433, and
> that is the same residue. A reviewer checking "H433" against our own coordinates would
> find a different residue there — mature 433 is not a histidine at all.
>
> The glycan sequon is quoted the other way round: **Asn420 is mature** (canonical Asn444),
> and 420 appears in the structures as 420. The mismatch is recorded as limitation 34 rather
> than corrected in place, because renaming ~110 occurrences hours before a deadline is a
> larger risk than the mislabel.

The EGFR ectodomain carries **17 histidines** (6ARU, apo). **Five of them lie in the
domain-III crop** our binders are designed against (mature 311–480); those five are the census
below. The other twelve are not tabulated because none of them can participate: the nearest
histidine outside the crop is **H507 (canonical; mature H483) at 37.5 Å** from H433's ring
nitrogens, against the 8.5 Å
of the H433–H370 pair, and the farthest is H183 at 83.9 Å. Measured over all sixteen
non-H433 histidines in the ECD structure, so the "only pair close enough" conclusion below is a
statement about all 17, not only about the five shown. (This line read "a full PROPKA census of
all 17 histidines" above a five-row table, which described the scope wrongly in the direction
of overstating it.)

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
rank-1 design, reads **5.656× on the histidine-only basis, above H433's 5.55× single-site
ceiling**, with H370 contributing 0.921–6.229× while H433 holds steady (§11.6). *(Two
corrections in this sentence: the value was 5.659× on five poses and is 5.656× on the twenty
now pooled, and the basis is the two-partner histidine-only gate, not the all-site product —
§11.7 explains why those are different quantities.)* That is one molecule with a wide pose
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
1.11× — inside PROPKA's noise — while dropping human ipSAE from 0.654 to **0.012** (§3.4; this
cited §3.7, which does not mention `S60D` — the design appears at this line, in §8.1 and in
§11.7, and nowhere else). H370
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

`d2c_mpnn13` is a **predicted** binding candidate with no supported pH switch (0.70×) — the reviewer's
own phrasing, 2026-10-03, which this sentence previously shortened to "binds and does not
switch". A single Ser→Asp at position 88 gives
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
| domain III crop (170 aa) | **0.2218** | **0.4005** |

*(This project writes the full ectodomain as 609 residues in some places and 621 in others,
against this section's own argument that a bar must come from the matching construct.
**621 is correct** for what was actually folded, and it is also what the organisers assay:
their construct is Met1–Ser645, and with EGFR's 24-residue signal peptide removed that is 621
residues of mature protein. The 609 figure is a stale earlier crop and appears nowhere in the
scoring path. An earlier note here blamed §3.1 for it; §3.1 gives patch-recovery counts and no
residue figure at all, so the pointer was wrong as well as the number.)*

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

**The cross-reactivity requirement does NOT rescue this control.** *Corrected 2026-10-05; the
previous heading here read "what rescues it" and the paragraph below it argued the opposite of
what the data supports.* The reviewer's correction was explicit: *"The mouse scores improve
separation in this panel, but a zero predicted interface is not an experimentally demonstrated
specificity mechanism. Without matched mouse outcomes, this does not validate mouse binding or
rescue the failed human control criterion."*

What is true: both molecules that outrank the binder on the human leg score exactly **0.0000 on
the mouse leg**, so scored as the submission is scored — requiring both species — the one
quantified binder ranks above all ten no-KD molecules on this panel.

What that does not establish, and the distinction is the whole point:

- **There are no matched mouse outcomes.** Adaptyv measured these molecules against human EGFR.
  We have no experimental mouse result for any of them, so a mouse prediction of 0.0000 is not
  corroborated by anything. The dual-species filter is being validated against data that does
  not exist for the species doing the work.
- **A zero predicted interface is not a specificity mechanism.** §4.4b of this document shows a
  pose that reproduces 72% of a crystallographic interface scoring ipSAE_min 0.0000. On this
  instrument 0.0000 means "no confident interface found", which is a statement about the
  predictor, not about the molecule. Reading it as "this molecule does not bind mouse EGFR" is
  the same error as reading a missing KD as a measured non-binder.
- **The human criterion failed and remains failed.** Two of ten censored molecules outrank the
  only quantified binder on the human leg. Adding a second species on which we have no outcomes
  cannot convert that into a pass; it can only show that a different, unvalidated filter happens
  to order this panel differently.

The earlier paragraph argued that this was "a mechanism and not a curve" because the exclusions
sit at exactly zero on five of five seeds rather than at a narrow margin. That argument is
withdrawn. Margin-independence would matter if the zeros were known to mean no interaction; §4.4b
shows they are not known to mean that. The dual-species criterion came from the brief rather than
from us, and the honest statement is that **no specificity filter in this pipeline is supported by
measured data** — including this one.

**gitter-yolo10 is the entire problem in one molecule.** Pooled over 5 refold poses it reads a
**5.27× pH ratio** — **18th of the 246 molecules eligible to rank** at the n ≥ 5 bar §6 sets, i.e. the top 7.3% — on our primary objective —
alongside a human ipSAE above every positive control we have. (An earlier version of this
document said "rank 8 of 2,009", pairing a rank computed among molecules with ≥5 poses with the
denominator of every sequence ever scored. That overstated it by about 15× and the derived
"top 0.4%" was wrong; 2,009 is the scored universe and 246 is the rankable one at n ≥ 5. A later version said "8th of 132", which does not reproduce from `master_rank.json` either — the rankable pool grew as poses were added, and the figure was never recomputed. 18 of 246 is what the current artifact gives.) A human-leg-only pipeline would have submitted a
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

### 4.4b rAC1 against its own co-crystal (4UIP): the instrument scores a correct interface zero

*Added 2026-10-05. The reviewer asked for this comparison to be done by structural contact
recovery and specifically not by attaching a predicted PAE to crystallographic coordinates —
a PAE describes a predictor's uncertainty about its own output, and a crystal has none.
`bin/rac1_contact_recovery.py` is purely geometric; no ipSAE or PAE value is attached to the
crystal anywhere in it.*

rAC1 is a published EGFR binder solved in complex with the receptor (PDB **4UIP**). We folded
it against both target constructs, 5 seeds each, and asked a geometric question: does the
prediction put the same residues in contact as the crystal? An interface contact is a residue
pair whose closest heavy atoms are within 5.0 Å. Residue numbering is mapped through a global
sequence alignment of each chain to its crystal counterpart, so the cropped construct still
lines up; 23 predicted contacts across all ten poses had no crystal counterpart and were
dropped and counted rather than treated as non-contacts.

The crystal interface is **65 residue–residue contacts over 28 epitope and 27 paratope
residues**.

| pose set | contact recall | precision | epitope recall | paratope recall |
|---|---|---|---|---|
| `ecd` seed 1 | **0.723** | 0.618 | **0.929** | 0.926 |
| the other 9 poses | **0.000** | 0.000 | 0.000–0.036 | 0.000–0.963 |

**One pose in ten reproduces the interface, and it reproduces it well** — 72% of the crystal
contacts and 93% of the epitope. The other nine recover **no** crystal contact at all, while
still predicting 31–121 contacts each: they dock the binder somewhere else. The paratope-recall
column is the interesting part of the failure. Several of the misdocked ECD poses present
**63–96% of the correct binder surface** while recovering **0%** of the correct epitope. The
predictor largely knows which face of rAC1 does the binding and puts it against the wrong face
of EGFR.

**And the instrument cannot tell the difference. All ten poses score ipSAE_min = 0.0000** —
including the one that is substantially correct. Ranked by ipSAE_min the correct pose sits at
position 6 of 10, which is an artefact of ties, not of discrimination: every value is exactly
zero.

Two consequences, both of which constrain how this submission may be read.

1. **A score of 0.0000 on this instrument does not mean "no interface".** It can mean "a
   correctly reproduced crystallographic interface that the predictor is not confident about".
   This matters directly, though not in the way this paragraph first claimed: it read
   *"several submitted designs carry 0.0000 on one species, and §12 of this document already
   warns that an absent measurement must not read as a measured zero."* **No** shipped design
   carries a 0.0000 on either species — the weakest assessable legs are **0.132** (human, `c5_cr_crop_patch…_05`) and **0.168** (mouse, `bcr_d3acid3_l60_s647537_mpnn3`); this read "0.132 and 0.181", which skipped the 0.168 leg — and
   §12 is Declarations and contains no such warning; it is **§4.2** that documents the error,
   cited by limitation 22. What stands is the point itself: here is a *measured,
   crystallographically-solved* interface reading 0.0000, so a zero on this instrument is not
   evidence of no interface. That is why the two antibody-format rows are reported as
   inadequately assessed rather than as low.
2. **It is a second instance of the G532 pattern** (CONTROL-TABLE §5b, a different
   document), now with
   geometry attached rather than inferred. A real binder, folded against its real target,
   scored at the instrument's floor. The earlier claim that the G532 result was "fully
   attributed to the pose" cannot be made here, because in this case the pose is **right** in
   one of ten tries and the score is zero anyway.

**The independent check: Chai-1 reproduces the failure exactly, and the CONSTRUCT turns out to
matter more than the predictor.** *Added 2026-10-05. Chai-1 (`runs/chai1/w4_indep/rac1_d3`,
5 models) folded the same rAC1 complex against the same domain-III crop, and the identical
geometric test was run on its structures. Chai emits no residue-level PAE, so no score
comparison is possible; this is contacts against contacts.*

| | contact recall | epitope recall | **paratope recall** |
|---|---|---|---|
| ESMFold2, d3 crop (5 poses) | 0.000 ×5 | 0.000–0.036 | **0.000–0.074** |
| **Chai-1, d3 crop (5 models)** | **0.000 ×5** | **0.000 ×5** | **0.000–0.074** |
| ESMFold2, full ECD (5 poses) | **0.723**, then 0.000 ×4 | **0.929**, then 0.000 ×4 | **0.370–0.963** |

**ESMFold2 is not the weak link here.** An architecturally independent model — different
weights, different training, diffusion co-folding rather than a folding trunk — reproduces its
crop failure to within rounding on every column. Whatever is wrong is not specific to
ESMFold2, which removes the most convenient explanation for the zero-scoring positives. (This
read "ESMFold2-Fast". `biomodals/modal_esmfold2.py:41` defaults `ESMFOLD2_HF_REPO` to
`biohub/ESMFold2` — the **Full** model — and `bin/score-esmfold2.sh` sets no override, so every
pose in this project is Full. Naming the Fast model understated the arm.)

**Chai-1 also folded six of the shipped finalists, and all eleven complexes it ran are
reported here.** The arm produced results for six submitted designs and five calibration
complexes; an earlier version of this document reported one of the six, which is the selective
reporting §13 warns about.

<!-- GENERATED:CHAI-TABLE do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
| complex | Chai-1 ipTM (median of 5) | interface residues | clashing models |
|---|---|---|---|
| **`ss_bc_s831683_mpnn6_S15D_S62H_routeA`** (rank 4) | **0.838** | 39 | 0 |
| **`d2c_mpnn13_S88D_serasp`** (rank 5) | **0.818** | 30 | 0 |
| **`rimA02_d3_rimA_14_vhh`** (rank 11) | **0.440** | 22 | 0 |
| **`rimA01_r15_boltzgen_egfr_d3_rimA_20`** (rank 3) | **0.332** | 24 | 0 |
| **`c5_cf_short__boltzgen_egfr_cropfree_short_48`** (rank 1) | **0.201** | 22 | 0 |
| *— calibration and reference complexes —* | | | |
| `barnase_barstar` | 0.877 | 22 | 0 |
| `cetuximab_scfv_ecd` | 0.793 | 29 | 0 |
| `fin_bc_s360518_mpnn9_A22D` | 0.788 | 39 | 0 |
| `egf_hu` | 0.500 | 25 | 0 |
| `g532_ecd` | 0.340 | 35 | 0 |
| `nano2_ecd` | 0.167 | 44 | 0 |

Six of the 16 shipped designs were folded by Chai-1, and it does **not** rate them alike: `ss_bc_s831683_mpnn6_S15D_S62H_routeA` reads 0.838 against `c5_cf_short__boltzgen_egfr_cropfree_short_48` at 0.201, a spread of 0.637 ipTM across designs our own pH objective orders quite differently. 2 of the six sit at or above the cetuximab scFv positive control (0.793): `ss_bc_s831683_mpnn6_S15D_S62H_routeA`, `d2c_mpnn13_S88D_serasp`. 3 sit **below human EGF** (0.500): `rimA02_d3_rimA_14_vhh`, `rimA01_r15_boltzgen_egfr_d3_rimA_20`, `c5_cf_short__boltzgen_egfr_cropfree_short_48`. Chai emits no residue-level PAE, so ipSAE cannot be computed on these and ipTM is not comparable to our ranking metric. It is a second opinion on whether an interface forms at all, not a second measurement of the objective.

**And the calibration set says not to over-read it.** `g532_ecd` is a PUBLISHED, experimentally-confirmed pH-switchable EGFR binder, and Chai-1 scores it **0.340** -- below 3 of our six designs and well below human EGF. `nano2_ecd` reads 0.167 on 44 interface residues, the largest interface in the set and the lowest score. So a low Chai ipTM is **not** evidence that a design does not bind: on the one molecule here with a real measured answer, this metric is wrong. The table supports the positive direction only -- three designs form an interface an independent predictor rates at or near the level of the cetuximab control -- and it cannot be used to argue against the designs at the bottom, including rank 1. Reporting it the other way round would be the single most tempting over-read available in this submission.
<!-- /GENERATED:CHAI-TABLE -->

**What separates the arms is the target construct, and the separation is 5-for-5 rather than a
lucky draw.** On the full ectodomain every pose presents **37–96% of the correct binder face**,
and one pose gets the whole interface right. On the crop, paratope recall never exceeds **7%** in
either predictor — the binder is docking by a different face entirely. The single 72%-recall
success is the weaker evidence here; the consistent, non-overlapping paratope-recall separation
between crop and ECD is the stronger.

**And it is not because the epitope is missing from the crop.** Measured: the 4UIP epitope is
28 residues spanning mature 411–489, of which **26 (92%) lie inside the crop's 311–480 window**
and 2 fall in domain IV. The crop contains almost the whole epitope and both predictors still
miss it.

**Why this matters to the submission: 15 of the 16 submitted designs were scored against that
crop** — every design except `d2c_mpnn13_S88D_serasp`, which is scored against the full
ectodomain and sees all 17 target histidines rather than the crop's 5. (This read "16 of the
18"; the count is derived from each design's own target-histidine census.) The one molecule in this project with a solved complex is never docked correctly on the
crop by either predictor, while the full ectodomain recovers it once in five and gets the binder
face approximately right every time. We are not able to rescore the submission on full ECD
before the deadline — one design already uses it, the other 17 would need refolding and
re-gating — so this is recorded as a limitation on the construct rather than fixed.

What none of this establishes: **n = 1 molecule.** One co-crystal cannot measure how often either
predictor docks correctly, cannot establish a general crop effect, and rAC1 is a 238 aa
non-antibody scaffold unlike most designs here. The paratope-recall separation is consistent
across ten poses but they are ten poses of one molecule.

### 4.5 The matched null is degenerate, so the screening flag does almost nothing.

12 shuffle nulls — 4 designs × 3 independent shuffles, composition and length preserved
exactly, 5 seeds nested within each shuffle — score **0.0000 on both legs with 5 of 5 dead
seeds, all twelve of them.** A separate, larger shuffle run (22 molecules, 110 poses) agrees:
**20 of 22 are dead on all 5 seeds and the highest median in it is 0.0110.** Across 34 shuffled
sequences the null does not produce a single score worth a threshold. This was predicted on
review: full-sequence shuffling destroys the fold, so the null bounds the instrument's noise
floor at exactly zero, and a "95th percentile of a matched null" flag admits anything above
0.0000. **The `affinity_above_null` column has therefore been REMOVED from the emitted CSV
(2026-10-05).** A threshold taken as the 95th percentile of a null that is a point mass at zero
is a percentile of nothing, and a column of that kind sitting in the one file Adaptyv grades
invites exactly the reading it cannot support. Affinity is reported in full as the two
continuous `ipsae_min` columns instead, and `AFFINITY_FLAG_DEPRECATED` survives in the code
only to preserve the record of what was once reported. **The 10 no-KD molecules are the only
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

**WITHDRAWN 2026-10-05: the control-validation claim.** This paragraph read: *"The controls
validate it independently of us: all six positive controls fall in the lower half of the design
distribution, and both EGF-derived 'nonbinders' fall in the top 14% — a model that never saw this
target says protonating H433 is maximally bad for an EGF-like complex, which is exactly right
for a neutral-pH agonist."* It is withdrawn for three reasons, each sufficient on its own:

1. **The reviewer ruled it out by name**, 2026-10-04 §5: *"The EGF-derived sequences' Potts
   ranks cannot independently establish their binding or agonism."* And §1: *"Relabel both
   EGF-derived controls as activity-unknown, document their provenance gap, and withdraw claims
   based on their supposed negative status."* §4.1 carried out the relabelling; this paragraph
   then kept arguing from the status that had just been withdrawn.
2. **It contradicts the two sentences immediately above it**, which say *"Agreement between two
   methods is evidence, not validation. Neither has experimental ground truth on this target."*
   A paragraph cannot disclaim validation and then claim it four lines later.
3. **No published artifact supports the numbers.** `potts_ddg_pool.json` and
   `potts_ddg_extra.json` hold 2,092 records between them and contain **no control molecule** —
   no positive control, no cetuximab, no EXPNEG molecule. The only `g532` labels are our own
   mimic designs. So neither "lower half" nor "top 14%" is reproducible from the repository,
   which is the standard every other number in this document is held to.

What remains, and it is the weaker claim: the Potts model and the PROPKA refold gate agree at
an odds ratio of 6.32, on our own designs, with no experimental ground truth on either side.

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
| a rank among rankable molecules quoted against a denominator of 2,009 | overstated by ~8×; the rankable pool is 246 at n ≥ 5, and the "132" that replaced 2,009 does not reproduce either (§4.4) |
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
BindCraft pool, which supplies 10 of the 16 submitted designs. Nothing in the pool
reaches Level 4 under the general-protein rule; every design that clears does so through the
Level 3 clause (moderate structural similarity *or* >30% sequence identity, exactly one of
them), which is worth saying plainly: **this is a pool of partly novel designs, not de novo
ones.**

**Provenance of this count, because an internal review disputed it.** The 292 rows across the
general-rule FoldSeek TSVs deduplicate to **238 distinct designs, 114 clearing**, keying on the
design basename. Stripping a trailing `_modelN` collapses two more, giving 236/113; no key gives
233/109. That figure is reproducible but stale — it is this same count computed **before
`novelty_s15d.tsv` existed**, and that file holds the five S15D-family designs (`bc_s831683_mpnn19_S15D`, `bc_s831683_mpnn6_S15D`, `bc_s831683_mpnn8_S15D`, `bc_s831683_mpnn9_S15D`, `bc_s831683_mpnn9_WT`)
that are in this submission. 22 designs appear in more than one TSV, and 2 of those carry slightly
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
to eligible. Novelty is also a severe filter on this target — of the 53 round-2 designs that actually BOUND,
**exactly one** would have passed the strict reading, i.e. **2%**, against 4% for non-binders,
so novelty is mildly anti-correlated with binding here. (That calibration is stated and sourced
in `bin/novelty_gate.py`, measured against Adaptyv's own pre-computed FoldSeek results for the
393 tested round-2 EGFR designs. This sentence attributed it to §9, which contains no round-2
measurement and no 2% figure of any kind.) So a high ratio is no guarantee any of these is
submittable.

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
submission was carrying **12 of the 20 permitted designs** when this was written, so eight slots
were unused. It now carries **18 of 20** — two of the molecules tabulated below were
subsequently shipped, which is what this section was for. Three of the ten failures are flagged
ANTIBODY and were scored by the general-protein rule, which §9 notes is *stricter* than the
antibody rule — including `ss_rimA02_d3_rimA_14_vhh_T28H_routeA` at 3.209×, a variant of
`rimA02_d3_rimA_14_vhh` (shipped, now rank 11) — so those three
need an ANARCI re-check before being treated as excluded.

The cases that most directly contradicted the ranking as it stood on 2026-10-04, all Level 3.
**Two were acted on and are now in the submission**, marked ✔ with their current rank; the
other four remain excluded and the contradiction they pose stands:

| molecule | pH (his-only) | human | mouse | spread | compare |
|---|---|---|---|---|---|
| ✔ `ss_bc_s831683_mpnn6_S15D_S62H_routeA` **(now shipped, rank 4)** | **3.545×** | 0.765 | 0.744 | 0.358 | shipped alongside its parent `mpnn6_S15D` (1.835×) as a declared pair |
| `sd_d2c_101_l147_s144898_m_T65D` † | **4.735×** | 0.585 | 0.366 | 0.944 | shipped `d2c_mpnn13_S88D`: 3.526×, 0.603, 0.528 |
| `ss_bc_s831683_mpnn19_S15D_S62H_routeA` | **3.189×** | 0.642 | 0.601 | 0.591 | shipped `mpnn19_S15D`: 1.774×, 0.808, 0.786 |
| `bcr_d3acid_l65_s831683_mpnn3_S15D` | 2.035× | 0.759 | 0.754 | 0.512 | beats shipped `mpnn6_S15D` at matched affinity |
| `bcr_d3acid_l65_s831683_mpnn17_S15D` | 1.814× | 0.780 | 0.755 | **0.019** | the most reproducible pH measurement in the project |
| ✔ `c5_cf_short__boltzgen_egfr_cropfree_short_48` **(now shipped, rank 1)** | **5.546×** | 0.241 | 0.181 | 0.262 | weak on both species, and §11 records four independent strikes against it |

`ss_bc_s831683_mpnn6_S15D_S62H_routeA` was the clearest: a shipped design plus one further
mutation, with predicted affinity indistinguishable from it (0.765/0.744 against 0.780/0.751)
and **1.9× the pH ratio**, at an acceptable spread and Level 3. On this submission's own
criteria it dominated a design we shipped — **so it was added**, at rank
5, with its parent retained beside it as a declared
parent/mutant pair (§11.3). This paragraph argued for a change that has since been made; it is
left in because the reasoning is the record of why.

**† Both of these were added and then removed the same day, and the reason matters.** They
were the top two of the five additions on the pH objective, and both are **near-identical to a
design already shipped** — 0.985 and 0.986 respectively. The review's instruction was *"avoid
filling available slots with nearly identical variants"*, so ranking the reopened pool on the
objective produced precisely what we had been told not to do. They were replaced by
`bcr_d3acid3_l60_s647537_mpnn3` (2.914×) and `bcr_d3acid3_l60_s647537_mpnn11` (2.747×), which
open a backbone family that had no representation. **Correction 2026-10-05:** these two are **0.867 identical to each other** — both 60 aa, differing at 8 positions — so they are MPNN redesigns of one backbone, not two independent designs. 0.467 is their identity to everything *else* submitted. The same false figure shipped in both CSV rows and is corrected there. The rows above are left in place because they are what the reopened ledger
found; they are no longer what the submission contains. §11.3 records the swap.

Note also the family composition: **15 of the 25 are the `d3acid3_l60_s647537` backbone**, a
family with no representation in the submission at all, while six of the twelve then-shipped
designs sat on a single other backbone.

**What remains unresolved, stated as such.** This section does not change the submission, and
the decision whether to use the then-unused slots is not made here — 2 remain unused
now, and six of the molecules below were subsequently added. Three things are known and
recorded: the exclusion of 63 molecules rested on a gate that cannot carry that weight; 25 of
them are eligible and outrank a shipped tier-1 design on the submission's own basis; and all of
these comparisons inherit the provisional status of that basis, since §11.7 measures Kendall
τ = +0.046 between it and the partnered alternative. Expression QC remains unrun for every one
of the 25. The original exclusion list and this amendment are both preserved and timestamped,
and every number here was produced before any outcome was examined.


## 10b. Finalist epitope footprints: the four checks, and one finding they surface

*Added 2026-10-05 at the reviewer's request: "apply full-ECD, glycan, receptor-state and
human/mouse contact checks to the actual finalist footprints — my earlier domain-II assessment
does not clear these designs." He was right that it does not. `bin/finalist_footprints.py`
computes each design's footprint from its own human-leg poses — a residue is in the footprint if
it contacts the binder within 5.0 Å in a **majority** of that design's poses — and maps it to
mature ECD numbering (canonical P00533 25–645 → mature 1–621) by sequence alignment, so
d3-crop and full-ECD poses land on one coordinate system.*

| check | result across all 18 |
|---|---|
<!-- GENERATED:FOOTPRINT-TABLE do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
| **domain** | **every design, 100% of contacts, in domain III (L2)** -- no domain-II contact anywhere |
| **full-ECD** | **0 of 16** have any contact outside the 170 aa domain-III crop (mature 311-480), so the crop is adequate and no footprint required the full ECD to assess |
| **glycan** | **3 of 16** touch an N-glycosylation sequon, all of them Asn420: `c5_cf_short__boltzgen_egfr_cropfree_` (rank 1), `bcr_d3acid3_l60_s647537_mpnn3` (rank 7), `bcr_d3acid3_l60_s647537_mpnn11` (rank 8). Of 11 sequons in the construct, only 1 is contacted |
| **human/mouse** | median identity **at the contacted positions** is **0.86**; range 0.77-0.92 over 16 designs |
<!-- /GENERATED:FOOTPRINT-TABLE -->

**Domain II is not in play.** The earlier assessment concerned domain II; these binders do not
touch it. That resolves the question in the designs' favour but by irrelevance, not by passing.

**The glycan flag is on three designs, including rank 1.**
`c5_cf_short__boltzgen_egfr_cropfree_short_48` (rank 1), `bcr_d3acid3_l60_s647537_mpnn3`
(rank 7) and `bcr_d3acid3_l60_s647537_mpnn11` (rank 8) all contact **Asn420**, one of eleven
N-X-S/T sequons in the human ectodomain and one of four in domain III (N328, N337, N389, N420).
It is the only sequon any design touches. None of our folded structures carries a glycan, so
that contact is made against a surface that is glycosylated in a real cell and bare in every
structure we scored. This is a liability on the design the pH objective ranks first and on two
more — though the two `bcr` rows are an 86.7%-identical pair, so the three are two independent
exposures, not three. It was not visible before the footprints were computed, and this
paragraph said "the glycan flag is on rank 1" for a day after the count became three.

**Cross-reactivity is weaker at the epitope than whole-protein identity suggests.** The two
ECDs are highly similar overall, but **14% of contacted positions differ between human and
mouse** at the median, and the worst two are `cons_gap_h370_only__boltzgen_egfr_h370_018`
(0.77) and `bc_s360518_mpnn9_A22D` (0.78). Since the submission requires both species,
conservation at the epitope is the relevant quantity, and these are the designs most exposed
to a species difference the ipSAE columns will not show.

**Receptor state: the construct matches the assay, and I had this backwards for an hour.**
An earlier version of this section said 6ARU was "a ligand-bound receptor in the extended
conformation" and listed receptor state as unassessed. That was wrong, and wrong in the
direction of understating the submission's position. **6ARU is the cetuximab-Fab–EGFR
ectodomain complex** (`HEADER ... 6ARU`, X-ray, 3.20 Å), and cetuximab's mechanism is to bind
domain III and hold the receptor in the **tethered, autoinhibited** arrangement rather than to
permit the extended one. Verified from the construct itself: all four ectodomain modules are
present (domain I 162, II 145, III 170, IV 132 residues, span mature 4–612). The organisers'
own measurement spec gives the assay targets as the **full ectodomains, tethered** (Sino
Biological 10001-H08H and 51091-M08H). So the conformational state we folded against is the
state the assay measures, not a mismatch.

**What remains a real caveat is subtler, and it is new.** The domain-III surface in this
construct is a surface whose conformation was determined **with an antibody bound to domain
III**. Every one of the 18 designs binds domain III (above). So all 18 were docked onto a
geometry templated by a bound antibody in the same region they target, and the Fab was stripped
before folding without any relaxation of the surface it had been in contact with. We have not
quantified how much the domain-III backbone in 6ARU differs from an unliganded tethered
ectodomain, and we are not asserting the overlap between our designs' shared core epitope and
cetuximab's own epitope at residue level, because the Fab coordinates were removed from the
construct and we did not re-derive them. Both are checkable and neither was checked.

**The finding the checks surfaced: this submission has ten backbone families and one epitope.**
Across all 18 designs the union of contacted residues is only **62 distinct positions
(mature 316–474)**, and **19 residues are contacted by at least 80% of the
designs**: 325, 348, 349, 350, 353, 355, 357, 382, 384, 408, 409, 411, 412, 417, 418, 438, 440, 465, 467. (This read "all 17 designs ... 58 distinct
positions ... 20 residues", computed before the eighteenth design was added.) The
backbone diversity reported in §11.3 is real and the epitope diversity is close to nil — every
design is a different scaffold presented to the same patch of domain III.

That is a correlated-failure risk the family counts conceal. If this patch is the wrong patch —
glycan-shielded in vivo, occluded in the tethered state, or simply not a site where a pH switch
can be built — the submission does not fail in nine partly-independent ways, it fails once. We
are stating it rather than diversifying, because the deadline does not permit generating and
assessing a second epitope, and because the one thing worse than a concentrated submission is a
concentrated submission presented as a diverse one.

## 10c. The second-site arm: 12 paired attempts, and it does not work

*Added 2026-10-05. `runs/esmfold2/w2_ss_d3` and `w2_ss_ecd` (120 poses) were run to answer one
question — can a second titratable site be engineered into a design we already have? — and the
paired analysis had not been done. Twelve designs, each a single point mutation on a parent
that is in this submission. **routeA adds a histidine to the binder; routeB adds an acid.** Both
parent and child are scored on the same histidine-only gate over their own human-leg poses, so
the fold change is attributable to the mutation.*

| route | mutation | parent | child | fold | design |
|---|---|---|---|---|---|
| A | S62H | 1.835 | **3.545** | **1.93** | `ss_bc_s831683_mpnn6_S15D_S62H_routeA` |
| A | S62H | 1.774 | **3.189** | **1.80** | `ss_bc_s831683_mpnn19_S15D_S62H_routeA` |
| A | N65H | 2.101 | 2.014 | 0.96 | `ss_h370_020_vhh_N65H_routeA` |
| A | T45H | 2.101 | 1.628 | 0.77 | `ss_h370_020_vhh_T45H_routeA` |
| A | T28H | 4.838 | 3.209 | 0.66 | `ss_rimA02_d3_rimA_14_vhh_T28H_routeA` |
| B | L133E | 4.256 | **5.656** | 1.33 | `ss_rimA01_r15_..._L133E_routeB` |
| B | M42E | 0.627 | 0.814 | 1.30 | `ss_bc_s831683_mpnn9_WT_M42E_routeB` |
| B | M42E | 0.737 | 0.949 | 1.29 | `ss_bc_d3acid_l65_s831683_mpnn11_M42E_routeB` |
| B | S62E | 1.774 | 1.739 | 0.98 | `ss_bc_s831683_mpnn19_S15D_S62E_routeB` |
| B | S62E | 1.835 | 1.781 | 0.97 | `ss_bc_s831683_mpnn6_S15D_S62E_routeB` |
| B | L133D | 4.256 | 2.116 | 0.50 | `ss_rimA01_r15_..._L133D_routeB` |
| B | S60E | 3.526 | **0.747** | **0.21** | `ss_d2c_mpnn13_S88D_serasp_S60E_routeB` |

**The headline is negative. Median fold change 0.98 over 12 attempts; 5 of 12 improved at all.**
Adding a titratable site to a finished design is, on this evidence, a coin flip centred on no
effect. routeA: median 0.96, 2 of 5. routeB: median 0.98, 3 of 7. Neither route is better than
the other and neither is better than doing nothing.

**Two pieces of real structure inside that null.**

*Histidine addition works on protein scaffolds and failed on every antibody scaffold.* The two
successes are **the same mutation at the same position on two sister sequences** — S62H on
`s831683` `mpnn6` and `mpnn19`, giving **1.93×** and **1.80×**. That is a replicated effect, the
strongest form available without a wet lab. All three routeA failures are VHH scaffolds (0.66,
0.77, 0.96), consistent with §4.2 and §4.5 on antibody blindness, though here the failure is in
the pH gate rather than the affinity instrument.

*Acid addition can destroy a working switch.* `S60E` on `d2c_mpnn13_S88D_serasp` took a 3.526×
switch to **0.747×** — a 79% loss, and the single largest effect in the arm. `L133D` halved its
parent. Burying a charge with no counter-charge in reach is the failure mode of §11.6, and it
appears here twice.

**The selection effect this creates, stated because it is ours.** This submission contains
**two of the twelve** second-site designs — `rimA01_r15_L133E` and
`ss_bc_s831683_mpnn6_S15D_S62H_routeA` — and **both are among the five that improved**. We do not
submit the seven that did not. That is selection on the outcome, and the correct reading of those
two rows is *two successes out of twelve attempts at the same strategy*, not *a strategy that
works*. The arm's own median says it does not.

**So read each of the two against its own siblings.** Whoever grades `L133E` should know that
`L133D` — the same position, one methylene shorter — **halved** the parent (§11.6). Whoever grades
the S62H row should know that the same substitution **reduced** the ratio on three other
scaffolds, all of them VHH (0.66, 0.77, 0.96), and that the acid route at the identical position
(`S62E`) did nothing on either sister sequence (0.97×, 0.98×). S62H is the arm's one replicated
success — 1.93× and 1.80× at the same position on two sister sequences — and it is submitted
*with its parent* precisely so that claim is tested rather than taken on our gate's word (§11.3).

*Both of these rows moved during 2026-10-05.* `ss_bc_..._S62H_routeA` was added on pH rank,
removed as a near-identical variant, and restored as a declared parent/mutant pair; §11.3 records
the distinction that makes the third of those defensible where the first was not.

## 11. The submission

**16 designs, ranked on the two-partner histidine-only pH product.** Track 3 allows 20.
Twelve were submitted on 2026-10-04; **five were added on 2026-10-05 from the reopened
exclusion pool of §10**, by a rule fixed before the result was examined (this file's own
`rank_key` over the 25 eligible reopened molecules, capped at 2 additions per backbone and a
ceiling of 7 designs per backbone across the submission). They occupy ranks 1, 2, 5, 7, 8 and 9.
(The cap read "7 of 17 per backbone", which stated the then-current submission size rather than
the rule; the rule is the 7-design ceiling, and `d3acid_l65_s831683` sits exactly at it.) Nothing was displaced — the
submission was at 12 of 20 and the six additions use free slots. **Two** slots remain unused.
(This read "five additions at ranks 1, 2, 3, 6 and 8" with "three slots" left: there are six
additions, rank 3 is one of the original twelve, and §11.3 of this same document already said
two.)

*The tables and counts in this section are GENERATED from the emitted CSV by `bin/gen_methods_submission.py` (self-tested). They were hand-maintained through three submission changes in one evening and drifted badly — an audit found this section still describing eleven designs, its rank table omitting the twelfth, every rank above 7 off by one, and §12 attesting review of "all ten" sequences. The interpretive text is still written by hand; the numbers are not.*

### 11.1 The ranking basis changed, and it reordered everything

Until 2026-10-04 we estimated the pH ratio with a gate that measures only the **target's**
histidines. It never measured our own binders' titratable groups — and <!-- GENERATED:BINDER-HIS -- do not edit by hand; `bin/gen_methods_submission.py --write` -->
**eight of the sixteen submitted designs carry at least one histidine of their own**: `ss_bc_s831683_mpnn6_S15D_S62H_routeA` (4); `bc_d3acid_l65_s831683_mpnn11`, `bc_s831683_mpnn19_S15D`, `bc_s831683_mpnn6_S15D`, `bc_s831683_mpnn8_S15D`, `bc_s831683_mpnn9_S15D`, `bc_s831683_mpnn9_WT` (3 each); `d2c_mpnn13_S88D_serasp` (2). The other eight carry none.
<!-- /GENERATED:BINDER-HIS -->
(This count was wrong five times by hand — "six of the eleven", "seven of twelve",
"eight of the seventeen" from a pre-addition count carried forward, "ten of seventeen",
then "eight of the seventeen" again after the swap — so it is generated from
`ph_sensitivity.json` and stated by **name** rather than by rank.) Those get buried at the interface
and lose 1.5–2.5 pKa units, and by the same thermodynamic linkage of §1 that **opposes**
acid-tightening. We were counting the target's sites and ignoring ours.

`bin/ph_gate_multisite.py` composes over the **histidines** of both partners, with the free
leg taken by deleting the other chain in place — the mirror of the argument §2 makes for the
target leg. *Corrected 2026-10-05: an earlier version of this sentence said "every titratable
site", which the code did not do — it parsed HIS/ASP/GLU but added only histidines to the
product. The shipped column is a **two-partner histidine-only approximation**. §11.7 gives the
all-site and partnered alternatives, the resulting ranking instability, and why the
histidine-only order is nevertheless the one retained.* Measured over 165 poses, n = 5–26 per design (this read "75 poses, n = 5–11", which was the 12-design pose set):

<!-- GENERATED:BASIS-TABLE -- do not edit by hand; `bin/gen_methods_submission.py --write` -->
| design | target-only | **his-only (graded)** | binder histidines | worst drag |
|---|---|---|---|---|
| rimA01_r15_L133E | 4.619 | **5.656** | 0 | — |
| c5_cf_short__boltzgen_egfr_cropfree_short_48 | 5.819 | **5.546** | 0 | — |
| rimA02_d3_rimA_14_vhh | 5.186 | **4.838** | 0 | — |
| c5_cr_crop_patch__boltzgen_egfr_crop_patch_05 | 5.265 | **4.812** | 0 | — |
| rimA01_r15_boltzgen_egfr_d3_rimA_20 | 4.582 | **4.256** | 0 | — |
| ss_bc_s831683_mpnn6_S15D_S62H_routeA | 5.386 | **3.545** | 4 | 0.66 |
| d2c_mpnn13_S88D_serasp | 4.572 | **3.526** | 2 | 0.979 |
| cons_gap_h370_only__boltzgen_egfr_h370_018 | 3.478 | **3.180** | 0 | — |
| bcr_d3acid3_l60_s647537_mpnn3 | 3.154 | **2.914** | 0 | — |
| bcr_d3acid3_l60_s647537_mpnn11 | 2.997 | **2.747** | 0 | — |
| bc_s831683_mpnn6_S15D | 5.397 | **1.835** | 3 | 0.661 |
| bc_s831683_mpnn19_S15D | 5.435 | **1.774** | 3 | 0.66 |
| bc_s831683_mpnn9_S15D | 5.428 | **1.062** | 3 | 0.339 |
| bc_s831683_mpnn8_S15D | 5.461 | **1.023** | 3 | 0.333 |
| bc_d3acid_l65_s831683_mpnn11 | 4.010 | **0.737** | 3 | 0.359 |
| bc_s831683_mpnn9_WT | 3.522 | **0.627** | 3 | 0.338 |
<!-- /GENERATED:BASIS-TABLE -->

**Every binder histidine moves down except the one we designed to move up.** Of **236**
binder-histidine site-pose readings across all eighteen designs, **231 sit below 1.0** (0.33 to
0.98) and **5 are above it — all five of them `ss_bc_s831683_mpnn6_S15D_S62H_routeA`'s
`binder:HIS62`, in every one of its five poses, at 1.647 to 2.303×.** That is the site the S62H
install was built to create, and it is the only binder histidine in the submission that favours
the acid state. An earlier version of this sentence read "0.33 to 0.98, none up", which was
wrong in the submission's own disfavour: it suppressed the one engineered site that worked. PROPKA noise would scatter both
ways. The gate's counter-charge guard fires on nearly all of them (nearest opposite charge
6.9–8.8 Å), so these are desolvation shifts with no electrostatic partner: the same mechanism as
the 0.702× steric floor of §6, and the same physics that defeated mechanism A (§3.5).

**This is not a different objective. It is a less wrong estimate of the same one** —
KD(7.4)/KD(6.5). The superseded number ships as its own CSV column so the change is auditable
rather than silent.

**Four of the eighteen are no longer switches — on the deletion free leg, and that
qualification matters.** §11.7 shows this count falls to **one** on a separately-folded apo free
leg, so it is a property of how the free leg is estimated rather than of the designs. Read the
four as the weakest on every basis computed, not as designs shown not to switch. The bottom four rows fall below the 1.20× bar that §6 sets
as PROPKA's noise floor, so they carry no pH claim and are ordered by mouse affinity instead.
`bc_s831683_mpnn8_S15D` led this submission at 5.461× on the superseded target-only basis and
is now last, at 1.023× on the graded his-only basis.

**What we are NOT claiming.** The binder's free leg comes from deleting the target in place, so
the isolated binder is not relaxed — a histidine buried in the complex may be solvent-exposed in
the real free binder, which would make its free pKa wrong. PROPKA on a buried histidine is its
hardest case. The *direction* is consistent across 76 poses and mechanistically coherent; the
*magnitude* is not established.

### 11.2 Final ranks

<!-- GENERATED:RANK-TABLE -- do not edit by hand; `bin/gen_methods_submission.py --write` -->
| rank | design | class | family | aa | **pH his-only (ranked)** | all-site | partnered | rank range | pose spread | target-only | poses | human | mouse | affinity assessable |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `c5_cf_short__boltzgen_egfr_cropfree_short_48` | single_chain | cf_cropfree_short (c5) | 70 | **5.546** | 5.685 | 5.779 | 2-10 | 0.26 | 5.819 | 6 | 0.241 | 0.181 | yes |
| 2 | `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05` | single_chain | cr_crop_patch (c5) | 66 | **4.812** | 4.849 | 5.285 | 4-11 | 0.44 | 5.265 | 11 | 0.132 | 0.204 | yes |
| 3 | `rimA01_r15_boltzgen_egfr_d3_rimA_20` | single_chain | rimA01_r15_d3_rimA_20 | 150 | **4.256** | 34.534 | 4.843 | 3-14 | 0.12 | 4.582 | 6 | 0.594 | 0.567 | yes |
| 4 | `ss_bc_s831683_mpnn6_S15D_S62H_routeA` | single_chain | d3acid_l65_s831683 | 65 | **3.545** | 3.643 | 8.303 | 4-11 | 0.36 | 5.386 | 5 | 0.765 | 0.744 | yes |
| 5 | `d2c_mpnn13_S88D_serasp` | single_chain | d2c_101_l147_s144898 | 147 | **3.526** | 6.681 | 6.663 | 5-9 | 0.48 | 4.572 | 5 | 0.603 | 0.528 | yes |
| 6 | `cons_gap_h370_only__boltzgen_egfr_h370_018` | single_chain | h370_018 (gap) | 90 | **3.180** | 27.774 | 3.859 | 4-15 | 0.17 | 3.478 | 11 | 0.457 | 0.215 | yes |
| 7 | `bcr_d3acid3_l60_s647537_mpnn3` | single_chain | d3acid3_l60_s647537 | 60 | **2.914** | 3.037 | 3.106 | 10-16 | 0.09 | 3.154 | 6 | 0.574 | 0.168 | yes |
| 8 | `bcr_d3acid3_l60_s647537_mpnn11` | single_chain | d3acid3_l60_s647537 | 60 | **2.747** | 2.918 | 2.960 | 11-17 | 0.06 | 2.997 | 6 | 0.443 | 0.234 | yes |
| 9 | `bc_s831683_mpnn6_S15D` | single_chain | d3acid_l65_s831683 | 65 | **1.835** | 5.949 | 12.995 | 2-13 | 0.04 | 5.397 | 5 | 0.780 | 0.751 | yes |
| 10 | `bc_s831683_mpnn19_S15D` | single_chain | d3acid_l65_s831683 | 65 | **1.774** | 5.486 | 11.930 | 3-14 | 0.50 | 5.435 | 5 | 0.808 | 0.786 | yes |
| 11 | `rimA02_d3_rimA_14_vhh` | nanobody | rimA02_d3_rimA_14 (VHH) | 129 | **4.838** | 4.976 | 5.183 | 3-13 | 0.30 | 5.186 | 6 | 0.219 | 0.447 | **no** |
| 12 | `rimA01_r15_L133E` | single_chain | rimA01_r15_d3_rimA_20 | 150 | **5.656** | 52.181 | 7.288 | 1-6 | 4.38 | 4.619 | 20 | 0.598 | 0.434 | yes |
| 13 | `bc_s831683_mpnn9_S15D` | single_chain | d3acid_l65_s831683 | 65 | **1.062** | 1.804 | 6.801 | 7-16 | 1.21 | 5.428 | 20 | 0.804 | 0.804 | yes |
| 14 | `bc_s831683_mpnn9_WT` | single_chain | d3acid_l65_s831683 | 65 | **0.627** | 1.195 | 5.226 | 12-18 | 1.76 | 3.522 | 26 | 0.786 | 0.784 | yes |
| 15 | `bc_d3acid_l65_s831683_mpnn11` | single_chain | d3acid_l65_s831683 | 65 | **0.737** | 1.841 | 7.685 | 5-17 | 0.64 | 4.010 | 6 | 0.796 | 0.784 | yes |
| 16 | `bc_s831683_mpnn8_S15D` | single_chain | d3acid_l65_s831683 | 65 | **1.023** | 1.644 | 6.789 | 8-17 | 0.16 | 5.461 | 5 | 0.776 | 0.764 | yes |
<!-- /GENERATED:RANK-TABLE -->

**Before reading this order, read this about rank 1.** `c5_cf_short__boltzgen_egfr_cropfree_short_48`
heads the table because `rank_key()` sorts tier 1 on the pH ratio with no affinity term. It
carries more independent problems than any other row, and they are documented separately in
five places, so they are collected once here:

- the **weakest assessable interfaces** in the submission — human 0.242, mouse 0.181 (limitation 25)
- it contacts the **Asn420 glycan sequon** (§10b; limitation 14), as do two other shipped designs
- the **lowest Chai-1 ipTM** of six finalists folded, 0.201 — though the calibration set forbids
  reading that against it, since a published pH-switchable binder scores 0.340 (§4.4b)
- **30.0% alanine**, a 7-residue poly-Ala run and GRAVY +0.46, the only design the expression
  gate flags (limitation 35)
- it is the **only design switching on H370** rather than H433 (§1), the site where 9 of 10
  incidental switchers did not bind and the one deliberate attempt destroyed the interface (§3.4)
- it has the **worst-folded binder monomer** in the slate, mean pLDDT 0.726, with 9.6% of its
  residues below 0.50 against a maximum of 2.2% anywhere else

It is not dropped, because the ratio is the stated objective and we rank on it rather than on a
post-hoc preference. It is flagged here because a reviewer reading only the order would
otherwise take it for our best design, and it is not.


**Assessable designs rank ahead of unassessable ones within tier 1.** On a pure pH ordering
`rimA02_d3_rimA_14_vhh` leads the submission at 4.838× — on a human ipSAE of 0.219 that we
cannot interpret, because §4.5 shows this instrument scores a measured 294 nM antibody **below
its own non-switching comparator**. We neither demote it on the pH axis nor score it at 0.0000
(the §4.2 error); we place it after the designs where both axes mean something.

**The cost, stated:** `rimA02_d3_rimA_14_vhh` carries the **third-highest** pH ratio in the
submission at 4.838× and sits at **rank 11**, below designs reading as low as 1.774×. (This
read "second-highest ... rank 6", both of which were true of a smaller submission.) If the organisers rank strictly on the primary
objective, this ordering costs us. It is a judgement that credible-interface-first is the more
defensible frame, following the reviewer instruction to apply eligibility and interface checks
before the challenge priorities — not a claim that rimA02 is worse.

### 11.3 Eight families, sixteen designs

<!-- GENERATED:FAMILY-LIST -- do not edit by hand; `bin/gen_methods_submission.py --write` -->
`d3acid_l65_s831683` **x7** (ranks 4, 9, 10, 13, 14, 15, 16) - `rimA01_r15_d3_rimA_20` **x2** (ranks 3, 12) - `d3acid3_l60_s647537` **x2** (ranks 7, 8) - `cf_cropfree_short (c5)` (rank 1) - `cr_crop_patch (c5)` (rank 2) - `d2c_101_l147_s144898` (rank 5) - `h370_018 (gap)` (rank 6) - `rimA02_d3_rimA_14 (VHH)` (rank 11)

**Effective n is 8 clusters, not 16 designs.** The largest cluster, `d3acid_l65_s831683`, holds 7 designs at ranks 4, 9, 10, 13, 14, 15, 16; 5 families contribute a single design each. Any interval must be computed on families, not designs.
<!-- /GENERATED:FAMILY-LIST -->

**Three pairs of submitted designs exceed 90% sequence identity, and all three are declared
parent/mutant comparisons.** Measured pairwise over all 18:

| identity | pair | mutation | why both ship |
|---|---|---|---|
| 0.993 | `rimA01_r15_boltzgen_egfr_d3_rimA_20` / `rimA01_r15_L133E` | L133E | measure the L133E effect in the laboratory (§11.6) |
| 0.985 | `bc_s831683_mpnn6_S15D` / `ss_bc_s831683_mpnn6_S15D_S62H_routeA` | S62H | measure the S62H effect — the second-site arm's only replicated result (§10c) |
| 0.985 | `bc_s831683_mpnn9_S15D` / `bc_s831683_mpnn9_WT` | S15D | measure the S15D effect; the WT is the matched parent |

In every case the parent ships alongside its mutant **so the mutation's effect is measured rather
than inferred from our gate**. Six of the eighteen designs sit in such a pair, and no pair is two
independent tests — that is what they are for.

**The distinction this rests on, because it was got wrong once today.** The review asked for both
halves of one sentence: *"Keep the experimental WT/mutant comparisons **and** avoid filling
available slots with nearly identical variants."* A near-identical variant is the first thing when
it ships as a declared comparison and the second thing when it is padding. On 2026-10-05 five
designs were added from the reopened pool of §10 ranked on the pH objective, and two of them —
`sd_d2c_101_l147_s144898_m_T65D` (0.986 to shipped `d2c_mpnn13_S88D_serasp`) and
`ss_bc_s831683_mpnn6_S15D_S62H_routeA` (0.985 to shipped `bc_s831683_mpnn6_S15D`) — were padding
by that test: they entered on rank alone, with no stated reason for the near-duplication. That
took the count to four pairs, two of them unintended. Both were removed and replaced by
`bcr_d3acid3_l60_s647537_mpnn3` and `_mpnn11`, which open a backbone family that had no
representation. **Correction 2026-10-05:** these two are **0.867 identical to each other** — both 60 aa, differing at 8 positions — so they are MPNN redesigns of one backbone, not two independent designs. 0.467 is their identity to everything *else* submitted. The same false figure shipped in both CSV rows and is corrected there.

`ss_bc_s831683_mpnn6_S15D_S62H_routeA` was then **restored on a different and stated basis**: as a
declared parent/mutant pair, for the same reason the other two pairs ship. The mutation is worth
a wet-lab test specifically — S62H is the only replicated result of the twelve-design second-site
arm, raising its parent **1.93×** here and **1.80×** at the same position on the sister sequence
`mpnn19_S15D`, against an arm median of 0.98× (§10c). It also carries the strongest independent
corroboration in the submission: Chai-1 ipTM **0.838** with a **39-residue** interface, the
joint-largest in the validation set and above both working positives (cetuximab scFv 0.793, human
EGF 0.500), while its predicted affinity is indistinguishable from its parent's (0.765/0.744
against 0.780/0.751). The submission therefore stands at **16 of the 20 permitted**, with four
slots deliberately unused rather than filled.

**On the eligibility rule.** The organisers state that *iterating on any previously submitted
design is explicitly disallowed*, stricter than the challenge page's "existing binder" wording.
Our reading is that "previously submitted" means an earlier round or upload — **nothing from
this project has been uploaded** — so all 17 are first-time submissions and the rule is not
engaged by the two remaining pairs. Under a stricter within-submission reading those two pairs
would be the exposure, and they predate every change made on 2026-10-05.

### 11.4 Reading conventions

* **all-site pH ratio** — median over every ESMFold2 refold pose of that exact binder sequence,
  composed over all titratable sites on both partners. A prediction, not a measurement.
* **ipSAE human / mouse** — median of 5 seeds, pooled by sequence from `master_rank.json`.
  **pH-agnostic.** Not interpretable at all for the two antibody-format designs,
  `rimA02_d3_rimA_14_vhh` and `h370_020_vhh` — see §4.5. (This read "for ranks 6 and 7", which
  were their positions in a 12-design submission; naming them instead of their ranks is the
  rule this document states and did not follow.)
* **`bc_s831683_mpnn9_WT` is a control, not a candidate**: the matched wild-type of
  `bc_s831683_mpnn9_S15D`, one residue apart. Its **his-only** product of 0.627× is itself
  informative — the unmutated backbone is predicted to bind *worse* in acid, which is the
  baseline the Ser→Asp install has to beat. The same is true of
  `rimA01_r15_boltzgen_egfr_d3_rimA_20` against `rimA01_r15_L133E`.
  *(Designs are referred to here by NAME rather than by rank. Rank references inside a
  rank-ordered file drift every time the file changes, and an audit found six of the twelve
  assessment strings citing the wrong design by rank for exactly that reason.)*
* `affinity_above_null` **was removed from the CSV on 2026-10-05.** It meant
  `max(hu, mo) ≥ 0.2218`, a threshold derived from a null that proved degenerate. It is
  deprecated historical metadata, not an interpretation. See §4.5.

### 11.5 Why not 20

The 20-design build is in the history. Ranks 11–20 of it did not stand on a measurement: eight
read below 1.0× even on the old, more generous basis. Filling the allocation would have bought
ten more wet-lab wells and no more evidence. A22D was added on 2026-10-04 because it is a
backbone not otherwise represented, carries the largest causal swing in the set (0.826 → 5.630
target-only; 0.574 → 3.738 all-site, so the swing survives the correction), and has the tightest
seed reproducibility in the project.

One consistency note, now resolved: this read *"`bin/check_discards.py` reads the 30-design
candidate JSON, not the emitted CSV, and still compares on the target-only ratio."* It was
changed on 2026-10-05 — `SUB_CSV = "submissions/01-egfr.csv"` — and now validates against the
shipped set on the graded basis, matching by sequence. The note is kept because it describes why
the earlier threshold was 1.263 instead of 2.289 and why 19 unshipped candidates were skipped.

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
every submitted design on the same 165 human-leg poses through one code path, so the
differences below are attributable to the composition rule alone:

- **his-only** — histidines on both partners. The shipped basis.
- **all-site** — every HIS/ASP/GLU on both partners, which is what §11.1 claimed.
- **partnered** — every site whose nearest opposite charge on the other chain is within
  6 Å. `PARTNER_CUT` was fixed before this analysis, not tuned to it.

The histidine-only values reproduced the shipped CSV exactly on all twelve designs then
submitted, which confirms the join and the characterisation above. The six designs added on
2026-10-05 were scored through the same three bases before being added, so all **18**
are on one footing; across the 18 the analysis covers **165 human-leg poses** with
**0 unassessed titratable sites**. (This read "the five designs ... all 17 ... 158 poses".)

**The ordering is not stable.** Kendall τ between the shipped basis and the partnered basis
is **+0.046** over the eighteen shipped designs (80 concordant pairs against 73 discordant of
153) — the two orderings are effectively uncorrelated. This read **+0.000**, which was the
figure over the twelve designs submitted at the time. Designs move by up to **12 ranks** (`bc_d3acid_l65_s831683_mpnn11`, 5-17 across the three bases)
(`rimA02_d3_rimA_14_vhh`: 2nd on his-only, 10th on partnered). Per-design rank ranges are in
the §11 table. **Every tier in this submission is therefore marked `provisional`, and no
order here should be read as established.**

**Why the shipped order is nevertheless retained.** The histidine-only value is the
**minimum of the three bases for all 18 designs** — re-checked by the emitter on every
run rather than remembered, and it reports 18 of 18. (This read "all seventeen
designs ... 17 of 17"; the emitter computes the count over every scored row, which is
18.) Ranking on it is ranking on the
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
only **1.141×** (its median ratio over all 20 poses; this read 1.30×, which is nearer H370's
own 1.327× swing and attributed the wrong quantity to the designed residue). This is independent support for the reviewer's point that L133E's
5.659 vs 5.554 is not evidence of a second site.

**Noise does not accumulate.** With 30–123 titratable sites per complex, a multiplicative
product invites the objection that it compounds PROPKA noise. Measured: the product over
sites the gate itself calls unmoved (|ratio−1| < 0.05) sits at **0.946–1.020** across the
twelve measured when this check was run. Sub-threshold noise cancels; the products are driven
by genuinely shifted sites.

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
`continue`. Across all **149 poses of the sixteen submitted designs** the unassessed count is **0** — though **five site-poses are flagged `implausible`**, a different guard: `target:ASP13` on `c5_cf_short…_48`, `binder:ASP45` on `bc_s360518_mpnn9_A22D`, `binder:ASP58` on `ss_bc_s831683_mpnn6_S15D_S62H_routeA`, and `binder:ASP101` on `rimA02_d3_rimA_14_vhh` in two poses. Each has an implied pKa_bound that could not be inverted inside the 0.5–9.0 window the gate allows for aspartate, so the site is reported and excluded rather than silently composed. So every site either entered the product
or was accounted for with a reason; none was dropped silently, which is the property this guard
exists to give. Site-level detail for every pose is retained in
`analysis/01-egfr/ph_sensitivity.json` so no later question requires a re-run.

**pKa-perturbation sensitivity: the pH ratio cannot order this submission.**
*Added 2026-10-05. The three-basis comparison above answers how much the **composition rule**
moves the answer. The reviewer asked a different question — "compare the deletion estimate with
consistently prepared apo/relaxed alternatives **and plausible pKa perturbations**" — and this
is the pKa half. `bin/ph_pka_perturbation.py` perturbs every site's stored pKa_free and
pKa_bound by independent Gaussian noise, recomputes each site's linkage, re-takes the
histidine-only product, re-takes the median over poses, and re-ranks all 18 (400 draws per
σ, re-run on 2026-10-05 at σ = 0.4, 0.8 and 1.2 and persisted as three artifacts). σ is PROPKA 3's
own reported RMSD (~0.8 pKa units, worse for buried residues), not a tuned value.*

<!-- GENERATED:SIGMA-TABLE do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
| sigma (pKa units) | keep baseline rank | span >= 5 ranks |
|---|---|---|
| 0.4 | 6 of 16 | 12 of 16 |
| **0.8 (PROPKA's own RMSD)** | **2 of 16** | **15 of 16** |
| 1.2 | 1 of 16 | 16 of 16 |
<!-- /GENERATED:SIGMA-TABLE -->

<!-- GENERATED:PERT-FINDINGS do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
**At PROPKA's own stated accuracy the ordering is not identifiable.** `d2c_mpnn13_S88D_serasp` spans ranks 1-18; `ss_bc_s831683_mpnn6_S15D_S62H_rout` spans ranks 1-17; `bcr_d3acid3_l60_s647537_mpnn3` spans ranks 2-15. The widest span is 17 of 16 ranks.

**The single most-stable design holds a top-three slot in 60% of draws** (`rimA01_r15_L133E`: rank 1 on the unperturbed pH ratio, which is the quantity being perturbed, and rank 12 in the shipped CSV, which also applies the pose-spread and antibody penalties; perturbed median 3). An earlier version of this section claimed no design exceeded 50%; it did, and it is the design the pH ratio puts first, so the error ran in the submission's favour. The conclusion does not depend on sigma: it already holds at the optimistic 0.4.

**What does survive.** Two things. First, the **bottom group is robustly at the bottom**: `bc_s831683_mpnn9_S15D` stays at rank 10 or worse in 95% of draws, `bc_d3acid_l65_s831683_mpnn11` stays at rank 10 or worse in 95% of draws, `bc_s831683_mpnn9_WT` stays at rank 14 or worse in 95% of draws -- 3 designs take a top-three slot in 0% of draws. "These are not switches" is stable under the noise. Second, a **top set exists even though its order does not**: 5 designs `rimA01_r15_L133E` (60%), `c5_cf_short__boltzgen_egfr_cro` (49%), `rimA02_d3_rimA_14_vhh` (43%), `d2c_mpnn13_S88D_serasp` (28%), `rimA01_r15_boltzgen_egfr_d3_ri` (26%) hold a top-three slot in at least 25% of draws, against 0-16% for the other 11.
<!-- /GENERATED:PERT-FINDINGS -->

**One thing this understates, in the submission's favour.** The perturbation moves the pH ratio
only, and the shipped order is produced by three sort keys, of which two — the pose-spread
penalty (`SPREAD_BAR`) and the antibody-unassessable penalty — do not depend on pKa at all.
The shipped order is therefore more stable than the table above, because those two keys pin
five rows regardless of the noise. The table is the right statement about *the pH ratio as a
ranking instrument*; it is not the full statement about the shipped order.

**The conclusion we draw, which is the reviewer's own instruction.** *"If rankings change
materially, use provisional tiers."* They change materially. Every tier in the CSV is already
marked `provisional`, and the defensible claim from this submission is **a top set and a bottom
set, not a rank order**. We have not collapsed the CSV to tiers because the platform takes an
ordered file, but no number in it should be read as placing one design above its neighbour.

**What remains unresolved.** The apo/relaxed half of the reviewer's request is not done. It
needs consistently-prepared unbound structures for both partners and new folding, and the
partner-deletion free leg remains a fixed-conformation diagnostic. Also unaddressed: multiplying
per-site ratios assumes **independent titration**, which is not a general treatment of coupled
sites, and these draws model pKa error as **noise** — a systematic PROPKA bias on buried
histidines would move every design together and this analysis would not detect it.

**Apo free leg vs partner-deletion free leg.** *Added 2026-10-05, the remaining half of the
reviewer's sensitivity request: "compare the deletion estimate with consistently prepared
apo/relaxed alternatives." `bin/ph_apo_freeleg.py` against `runs/esmfold2/w5_apo` — every
submitted binder folded ALONE by the same predictor, same 5 seeds, same pipeline (19 monomers,
95 structures). The **bound leg is identical in both**, so every difference below is
attributable to the free leg alone.*

This is a comparison of two approximations and the reviewer said so plainly: *"A separately
predicted apo structure is another approximation, not automatically the correct answer."*
Deletion holds the side chains in a conformation the free protein does not adopt; the apo fold
gives a plausible unbound conformation that has no particular relationship to the bound pose,
so pairing its pKa with the complex's bound pKa mixes two structures. Neither is the free
protein.

**Internal check first.** All **8** designs carrying binder histidines move; all **9** carrying
none move by exactly nothing. The free-leg choice affects only the designs whose own histidines
enter the product, which is what it should do and is evidence the comparison isolates what it
claims to.

**The ordering is largely preserved — Kendall τ = +0.868**, against τ = +0.000 for the
partnered composition basis. On this axis the submission is far more stable than on the
composition axis or under pKa noise.

**What is not preserved is the bottom of the table.**

| | deletion leg | apo leg |
|---|---|---|
| designs reading below 1.0× | 2 | **0** |
| designs reading below the §6 bar of 1.20× | 4 | **1** |

`bc_s831683_mpnn9_WT` goes **0.627 → 1.161** and `bc_d3acid_l65_s831683_mpnn11` goes
**0.737 → 1.260**: both cross from "not a switch" to "a mild switch" on nothing but the choice
of free leg. **So the statement elsewhere in this document that four designs are no longer
switches is not robust.** It holds on the deletion leg and largely dissolves on the apo leg,
and we cannot say which leg is right. What survives is weaker and should be read instead: those
four are *the weakest four on every basis we have computed*, not *designs shown not to switch*.

**And the direction of the effect is inconsistent within one backbone.** The six `s831683`
designs each carry three binder histidines. The apo leg moves four of them **up** by 1.71–2.02×
(`mpnn8_S15D` 1.023 → 2.065, `mpnn9_S15D` 1.061 → 1.863) and two of them **down** by 0.90–0.91×
(`mpnn6_S15D` 1.835 → 1.679, `mpnn19_S15D` 1.774 → 1.595). Same backbone, same number of binder
histidines, same intervention — opposite-signed response to the free-leg change. Whatever the
apo fold is doing to these histidines' burial, it is not doing it systematically, which is a
reason to distrust both legs at this level of resolution rather than to prefer one.

The largest single mover is `d2c_mpnn13_S88D_serasp`, **3.526 → 5.041 (+43%)**, which rises
from fifth to third. It carries two binder histidines and the full 621 aa ECD as its target.

**Note on scope.** The target's free leg is left on the deletion estimate throughout the table
above, because the target apo structure is one shared fold and mixing it per-pose would confound
the binder comparison. A bounded target-leg and side-chain-relaxed comparison is reported
separately in §11.8.

**What this is not.** These three numbers are a sensitivity analysis, not a confidence
interval. They bound how much the composition rule moves the answer; they say nothing about
whether PROPKA's pKa values are right, and the partner-deletion free leg remains a
fixed-conformation diagnostic rather than a measurement of the apo state.

### 11.8 Relaxed free leg, and the target's own free leg

*The third rung of the free-leg ladder the reviewer asked for. §11.7 compares the
partner-deletion estimate against a separately-folded apo structure; those two differ in
backbone as well as side chains, so neither isolates the effect he actually named —
**"side-chain relaxation ... not represented"**. This section isolates it.*

**Method.** `bin/relax_chain.py` takes each bound complex, extracts one chain, repairs it with
pdbfixer, and energy-minimises it under amber14 + GBn2 implicit solvent **with the backbone
harmonically restrained** (10 kcal/mol/Å² on N, CA, C, O). Two choices are deliberate:

- **The backbone is restrained.** Released, this becomes a slow refold and stops isolating
  side-chain relaxation — it would just be a worse version of the apo arm.
- **Implicit solvent, not vacuum.** In vacuum, surface polar side chains collapse onto the
  protein to satisfy their own electrostatics. That is precisely the burial change a pKa
  calculation responds to, so vacuum would manufacture the effect being measured.

So the three legs differ in exactly one controlled way each: deletion freezes everything, relaxed
frees the side chains only, apo changes the backbone too.

**Scope.** The binder chain of every human-leg pose of all 18 designs, plus the target chain for
a bounded three-design subset (`rimA01_r15_L133E`, `bc_s360518_mpnn9_A22D`,
`bcr_d3acid3_l60_s647537_mpnn3`) to test whether the target leg behaves like the binder leg.
191 relaxations in total. The target subset is bounded because a 170 aa target minimisation costs
~42 s against ~9 s for a 65 aa binder.

**Completed 2026-10-05 01:17.** 191 relaxations, 0 failures, 60.8 min wall-clock.

**Coverage, stated before the result.** Of the 18 shipped designs, the relaxed leg can only move
a design that carries a histidine on a chain that was relaxed:

| | designs | why |
|---|---|---|
| compared | 10 | a binder histidine, or a relaxed target, or both |
| relaxed, no movable site | 7 | **zero binder histidines** — the pH signal is carried by the target's own H370/H433, and a binder-only relaxation cannot touch it |
| no relaxed structure | 1 | `ss_bc_s831683_mpnn6_S15D_S62H_routeA` entered the submission after the relax queue was built — a real gap, not a filtered one |

The seven "no movable site" designs are not missing data and their deletion numbers are not in
doubt; they are simply outside this arm's reach. Reporting them as uncovered would overstate the
arm, and dropping them silently would overstate its coverage.

**Result 1 — the direction is one-way.** All ten designs move **down**:

| | fold (relaxed ÷ deletion) |
|---|---|
| median | **0.842** |
| range | 0.623 – 0.971 |
| moving ≥10% | 8 of 10 |
| moving up | **0 of 10** |

This is a bias, not scatter. Relaxing the free state lets a partially buried histidine's side
chain reorganise and recover part of its solvated pKa; the free-state pKa rises toward normal,
the bound-minus-free gap narrows, and the linkage ratio falls. **The shipped deletion basis is
therefore optimistic on the pH ratio — by about 16% at the median and up to 38% at the worst.**
It is stated here as a signed bias rather than a symmetric uncertainty, because that is what the
data show.

**Result 2 — the ordering survives.** Kendall τ between the deletion and relaxed orderings is
**+0.956**, and the largest single-design rank shift is **1**.

This is the opposite of §11.7's perturbation result, and the contrast is the point. Perturbing
the pKa values at PROPKA's own accuracy moves 17 of 18 designs by ≥5 ranks (§11.7's
generated σ table; this read "16 of 17", the pre-addition figure). Changing the
free-leg definition — a much larger conceptual change — barely moves the ordering at all. So the
ranking is fragile with respect to **pKa accuracy** and robust with respect to **free-leg
choice**. Those are two distinct axes, and only one of them scrambles the table. §11.7's
conclusion stands unchanged: the objective cannot order this submission. This section does not
rescue it.

**Result 3 — the target leg is not a special case.** The three designs whose target chain was
also relaxed move 0.85×, 0.86× and 0.89×, inside the binder-only range. There is no evidence that
relaxing the target behaves differently from relaxing the binder, which is the only claim the
bounded three-design subset can support. It cannot rule out a target-specific effect at designs
not in the subset.

**Result 4 — the three-rung ladder.**

| design | deletion | relaxed | apo | span |
|---|---|---|---|---|
| `rimA01_r15_L133E` | 5.656 | 4.816 | 5.656 † | 1.17× |
| `bc_s360518_mpnn9_A22D` | 3.738 | 3.209 | 3.681 | 1.16× |
| `d2c_mpnn13_S88D_serasp` | 3.526 | 2.938 | 5.041 | 1.72× |
| `bcr_d3acid3_l60_s647537_mpnn3` | 2.914 | 2.586 | 2.914 † | 1.13× |
| `bc_s831683_mpnn6_S15D` | 1.835 | 1.781 | 1.679 | 1.09× |
| `bc_s831683_mpnn19_S15D` | 1.774 | 1.683 | 1.595 | 1.11× |
| `bc_s831683_mpnn9_S15D` | 1.061 | 0.699 | 1.863 | 2.67× |
| `bc_s831683_mpnn8_S15D` | 1.023 | 0.704 | 2.065 | 2.93× |
| `bc_d3acid_l65_s831683_mpnn11` | 0.737 | 0.460 | 1.260 | 2.74× |
| `bc_s831683_mpnn9_WT` | 0.627 | 0.423 | 1.161 | 2.74× |

Span across the three free legs: **median 1.445×, maximum 2.934×**.

† **These two rows are a two-rung ladder, not a three-rung one.** The apo arm (§11.7) remaps only
*binder* histidines, and these designs have none, so their apo entry is the deletion number
restated rather than an independent estimate. Their 1.13–1.17× span comes entirely from the
relaxed column. Counting those as three-way agreement would be the same vacuous-agreement
artifact that `check_claims.py` was built to catch, so they are marked instead of averaged in.

For the s831683 family the two arms move in **opposite** directions — relaxed down to 0.42–0.70,
apo up to 1.16–2.07. The deletion basis sits between them rather than at an extreme, which is
mildly reassuring about the shipped choice and says nothing about which leg is right. None of the
three is the free protein.

**What this changes in the submission: nothing.** The graded column remains the deletion,
histidine-only basis, declared as such in its own column name. What this arm adds is a bound on
that basis's modelling error — **≈1.4× at the median, ≈2.9× at the worst, with the error signed
optimistic** — and the finding that the error does not reorder the table. Recorded as
limitations 23 (the signed bias) and 24 (the arm's reach).

## 12. Declarations

Stated because the organisers ask for them and an earlier version of this document made none.

**AI assistance.** This submission was produced by one person working with Claude (Anthropic)
throughout: design generation, scoring, analysis code, and the drafting of this document. Every
number here was computed by code in the published repository, and the code was written in that
collaboration. The errors in §7 were found the same way.

<!-- GENERATED:DECL-REVIEW do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
**Human review.** The submitting researcher has reviewed all **16** submitted sequences -- their `molecule_class` labels, their lengths, and the claims made about them in this document and in the CSV.

Of these, **10** were in the submission as it stood on 2026-10-04 and **6** were added on 2026-10-05. The additions are `c5_cf_short__boltzgen_egfr_cropfree_short_48`, `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA`, `cons_gap_h370_only__boltzgen_egfr_h370_018`, `bcr_d3acid3_l60_s647537_mpnn3`, `bcr_d3acid3_l60_s647537_mpnn11`.

Designs that were in the 2026-10-04 set and are **no longer submitted**: `bc_s360518_mpnn9_A22D`, `h370_020_vhh`.

What has been verified for the 6 additions by code, and is reproducible from the repository: each comes from this project's own generation runs (§10); each was re-scored on the same three pH bases over its own human-leg poses; and the provenance audit below covers them. What has **not** been done for them: expression QC. Measured rather than asserted -- `analysis/01-egfr/express_qc.tsv` joins to **16 of the 16** submitted designs, so 0 have no expression-QC row. Novelty IS established for all 16: `bin/check_novelty_coverage.py` is green, and the four designs that had no levelled record were re-run on 2026-10-05 (`analysis/01-egfr/novelty_gap4.tsv`). This sentence said the checker was RED, which it was for about an hour before the gap was closed.
<!-- /GENERATED:DECL-REVIEW -->

**Provenance.** All eighteen sequences are de novo designs from this project's own generation
runs; none is a modification of a previously submitted design or of an existing characterised
binder. This was checked by code over all eighteen, two ways: no submitted design has any
`g532` ancestry in its generation lineage, and the highest sequence identity of any submitted
design to any known binder (G532 heavy and light chains, G532V, G532Ctrl, G5V2, and rAC1 from
4UIP) is **39.7%**, a VHH framework match. The two published binders that appear anywhere in
this project — G532 and the cetuximab-derived controls — were used for calibration only and are
not ancestors of any submitted sequence.

**Tools and licences.** BoltzGen and BindCraft for generation; ESMFold2 (via Modal) for structure
prediction; the Dunbrack `ipsae.py` v4 reference, MIT, commit `6174cf9e` for interface scoring;
PROPKA 3.5.1 for pKa; Proton-PottsMPNN (`potts_v6_afdb_edge_his0.3_acid0.06`) for the protonation
model; FoldSeek against PDB for novelty; HBPLUS v3.06. All open-source and used within their
licences. The vendored reference implementation in `outbox/ipsae-fixtures/vendor/` carries its
MIT licence file.

<!-- GENERATED:DECL-STRUCT do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
**Structures.** Predicted complexes are published at `submissions/structures/`, one median-ipSAE pose each -- not the best pose, which would be selection on the outcome. **Coverage is 9 of 16.** Without a published structure: `c5_cf_short__boltzgen_egfr_cropfree_short_48`, `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA`, `cons_gap_h370_only__boltzgen_egfr_h370_018`, `bcr_d3acid3_l60_s647537_mpnn3`, `bcr_d3acid3_l60_s647537_mpnn11`, `rimA01_r15_L133E`. Their poses exist and are scored; they are simply not exported. Stated rather than implied.
<!-- /GENERATED:DECL-STRUCT -->

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
   degenerate — 12 of 12 at exactly 0.0000 (§4.5) — so the (now removed) `affinity_above_null`
   flag meant
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
11. **Ten backbone families, one epitope.** All 18 designs contact the same patch of domain
    III — **19** residues are shared by ≥80% of them, over a union of 62 positions (§10b;
    this read 20, the count before the eighteenth design was added). The backbone diversity in §11.3 does
    not buy epitope diversity, so a wrong epitope fails the whole submission at once rather
    than ten partly-independent times.
12. **The domain-III surface was templated by a bound antibody.** The construct is the full
    tethered ectodomain from 6ARU — which matches the organisers' assay spec — but 6ARU is the
    cetuximab-Fab complex, the Fab was stripped without relaxation, and all 18 designs bind the
    same domain III the Fab occupied. The deviation from an unliganded tethered ectodomain is
    unquantified (§10b).
13. **A 0.0000 on this instrument is not "no interface".** On its own co-crystal (4UIP), a
    pose reproducing 72% of the crystal contacts and 93% of the epitope scores ipSAE_min
    0.0000 — identical to the nine poses that recover no crystal contact at all (§4.4b).
14. **Three shipped designs contact a glycosylation sequon, not one.**
    `c5_cf_short__boltzgen_egfr_cropfree_short_48` (rank 1),
    `bcr_d3acid3_l60_s647537_mpnn3` (rank 7) and `bcr_d3acid3_l60_s647537_mpnn11` (rank 8)
    all contact
    **Asn420**; no structure we folded carries a glycan (§10b). This limitation read "rank 1
    contacts a glycosylation sequon" while §10b read "1 of 17", so the exposure looked like
    one unlucky row rather than a property shared by three designs — and the two bcr rows
    are an 86.7%-identical pair, so two of the three are not independent.
15. **The domain-III crop may be the wrong construct, and most of the submission uses it.**
    On the only molecule here with a solved complex, neither ESMFold2 nor an architecturally
    independent model (Chai-1) recovers a single crystallographic contact on the crop, while
    the full ectodomain recovers 72% of them once in five poses and gets 37–96% of the correct
    binder face in all five. 92% of that epitope is inside the crop, so absence is not the
    explanation. 17 of 18 submitted designs are scored on the crop (§4.4b); only
    `d2c_mpnn13_S88D_serasp` is scored on the full ectodomain.
16. **The "not a switch" verdict is free-leg dependent.** Four designs read below the 1.20×
    bar on the partner-deletion free leg and only one does on a separately-folded apo free
    leg; two cross 1.0× on that change alone (§11.7). The two legs order the submission
    consistently (τ = +0.868) but disagree about its floor.
17. **The pH ratio cannot order this submission.** Under PROPKA's own reported accuracy
    (±0.8 pKa units) 17 of 18 designs span five or more ranks, and the widest span is 17 of
    18 (§11.7). The most stable design holds a top-three slot in 60% of draws and the next
    two in 49% and 43%; everything below the top five holds one in at most 17%. A top set
    and a bottom set are defensible; a rank order is not. (This limitation previously read
    "no design holds a top-three slot in more than half of draws", which was false — one
    does, at 60%, and it is the top-ranked design.)
18. **Two shipped designs are selected successes from a failed arm.** The second-site strategy
    improved 5 of 12 attempts, median fold 0.98 (§10c). Two of the five successes are in this
    submission and none of the seven failures is. Read the S62H and L133E rows as two
    successes out of twelve attempts, not as a working method.
19. <!-- GENERATED:LIMIT-FAMILY do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
    **Effective n is 8, not 16.** 7 of the 16 submitted designs sit on one backbone (`d3acid_l65_s831683`, ranks 4, 9, 10, 13, 14, 15, 16), and 2 further families are two-design clusters: `d3acid3_l60_s647537` (ranks 7, 8); `rimA01_r15_d3_rimA_20` (ranks 3, 12). Any hit rate or interval computed over designs rather than sequence families overstates n by up to 7-fold on the arm carrying our only causal claim. See §11.3 for the partition.
<!-- /GENERATED:LIMIT-FAMILY -->
20. **The reproducibility boundary is one gate wide, measured rather than asserted.** In a
    fresh clone with no local state: the emit command reproduces the graded CSV byte-identically,
    and 12 of 13 gates pass. The exception is `check_discards`, which needs the unpublished pose
    cache to confirm a design was measured. This limitation previously read "the documented emit
    command returns zeroed affinity columns without it", which was true of an earlier state and
    had not been retested; two packaging faults in `bin/ipsae_min.py` were found and fixed by
    cloning the repository and running it (see the Repository note). The pose cache is still not
    published.
21. **The design family with the best measured prior is the one we scored least.** 60 of the 71
    recovered BindCraft sequences sit at **n = 1 pose**, below the n ≥ 5 floor §6 requires, so
    they are structurally ineligible for tier 1 regardless of merit. 22 of those read ≥ 2.0× at
    n = 1. BindCraft also beat BoltzGen on every axis we measured, on 11 invocations against
    1,944 designs. We did not resolve them, and the submission is poorer for it.
22. <!-- GENERATED:LIMIT-AFFINITY do not edit between these markers; python3 bin/gen_methods_submission.py --write -->
    **The organisers rank outcomes partly on affinity at pH 6.5, and 1 of the 16 submitted rows have no usable affinity reading at all** (§4.5): `rimA02_d3_rimA_14_vhh`. We submitted them anyway, because excluding them would mean scoring them at 0.0000, which is the error §4.2 documents — but it means 6% of the submission cannot compete on one of the stated criteria.
<!-- /GENERATED:LIMIT-AFFINITY -->
23. **The shipped pH basis is optimistic, by a measured and signed amount.** The relaxed free
    leg (§11.8) moves all ten testable designs **down** — median 0.842×, worst 0.623× — so
    freezing the free state in the bound conformation inflates the headline ratio. Across the
    three free legs the span reaches 2.934×. The graded column is not corrected for this,
    because none of the three legs is the free protein and picking one post hoc would be the
    reshaping §11.8 was published early to prevent. The ordering is unaffected (Kendall
    τ = +0.956); the magnitudes are not reliable to better than about 1.4×.
24. **The relaxed arm reaches 10 of 18 designs, and the gap is not random.** **Nine** of the
    eighteen designs carry no binder histidine (§11.1) — their switch is target-borne — and a
    binder-only relaxation cannot move them, so seven of those nine are outside this check's
    reach entirely. (This limitation said "seven designs" and then "the eight designs" for the
    same set, in one sentence; the count of designs with no binder histidine is nine.) One
    further design
    (`ss_bc_s831683_mpnn6_S15D_S62H_routeA`) entered the submission after the relax queue was
    built and has no relaxed structure at all. The 0.842× median therefore describes the
    binder-borne designs, not the submission.

25. **The pH objective puts the two weakest predicted interfaces at the top.** Ranks 1 and 2
    read ipSAE human/mouse of 0.242/0.181 and 0.132/0.204 — the two weakest assessable
    interfaces in the submission. `rank_key()` orders tier 1 by the pH ratio alone, with no
    affinity term, so this is the ranking rule working as specified rather than a bug. It
    means the submission's own top of table is where its structural evidence is thinnest.
26. **The one two-site causal claim rests on a shift inside PROPKA's own error.** L133E's
    H370 shift is **+0.540 pKa units at the median** over 20 poses (range −0.40 to +2.56),
    and only **5 of 20** poses exceed the ±0.8 unit accuracy §11.7 adopts for the
    perturbation study. The two-site reading (§11.6) is therefore supported by a minority of
    poses on a median effect smaller than the instrument's stated error. We report it because
    the mechanism is specifically controlled — Glu reaches, Asp does not — not because the
    magnitude is resolved.
27. **The reviewer's tethered-versus-extended footprint comparison was never run.** He asked for the
    complete binder footprint compared across tethered and ligand-bound extended assemblies,
    including the second receptor, glycans and membrane-facing orientation, and said
    explicitly that distance from one tether contact cannot settle it. §10b runs four other
    checks and does not run this one. The receptor state we model (6ARU, tethered) matches
    the assay construct, which is why the gap is tolerable, but it is a gap.
28. **"The crop is adequate" is argued from crop-docked poses.** §10b concludes no footprint
    required the full ECD because 0 of 18 contact outside the domain-III crop — but those
    footprints are computed on poses docked against that same crop, which cannot place a
    contact outside it. §4.4b's rAC1 result points the other way: on the full ectodomain the
    paratope recall is 0.37–0.96 and on the crop 0.00–0.07. The honest statement is that the
    crop is adequate *for the poses we generated*, which is not the same claim.
29. **No solvent-accessible surface area was computed.** The reviewer: *"The burial atom count is a
    useful proxy, not a substitute for solvent-accessible surface area."* It is still the sole
    support for the H370 burial conclusion. `biomodals/modal_sasa.py` exists and was never
    run; "SASA" appears in no deliverable.
30. **The seed-instability pilot the reviewer specified was never run.** He asked for 20–30 diverse
    candidates enriched near decision boundaries at ~10 seeds each, assessing rank changes,
    pose consistency and threshold crossings. What exists instead is 5 seeds per design on the
    shipped set and a 15-seed triad on one design (§11.6). The pKa-perturbation study (§11.7)
    answers a different question — instrument noise, not seed noise.
31. **Mechanism A was never tested, and was twice asserted to be ruled out.** the reviewer's 2026-09-29
    answer 3 was to give mechanism A most of the initial design effort conditional on finding
    a suitable local acidic surface, retaining B as a smaller branch. §1 records that we used
    mechanism B. No mechanism-A result is reported anywhere in this document, so the two
    statements that A was ruled out are not supported by an experiment we ran.
32. **Every shipped design clears novelty Level ≥ 3; the tightest margin is 0.0076.**
    `bin/check_novelty_coverage.py` is green at 18 of 18. It was RED earlier on 2026-10-05:
    four designs had no levelled record, and `bc_s360518_mpnn9_A22D`'s figure existed only in a
    commit message and in its own CSV prose. Commit 332b09e (2026-10-04 19:53) read *"Novelty
    re-measured on the MUTANT pose ... TM 0.792, 15.2% identity, Level 3 — clears by 0.008,
    so it is also the row most exposed to a domain-wise novelty rejection"* — and no novelty
    TSV anywhere in the repository contained it, so the number was unreproducible and that
    warning never reached this document.
    **Closed by re-running it rather than by declaring it unreproducible.** FoldSeek was
    reinstalled from the upstream static binary, the PDB database re-downloaded (2.2 GB
    transfer, 6.4 GB indexed), the four binder chains re-extracted from their **original
    unrelaxed** ESMFold2 poses by sequence, and `bin/novelty_gate.py` re-run. Results in
    `analysis/01-egfr/novelty_gap4.tsv`, committed: `bc_s360518_mpnn9_A22D` **TM 0.792**,
    the commit message's figure exactly; `bc_d3acid_l65_s831683_mpnn11` 0.778;
    `rimA01_r15_L133E` 0.644; `d2c_mpnn13_S88D_serasp` 0.598. All Level 3. Measured identity
    for `bc_s360518_mpnn9_A22D` is **0.141** against the 15.2% quoted, which does not change
    its level.
33. **Seven designs clear novelty with under 0.04 TM of margin; the tightest is 0.0076.**
    The level-2 cliff is TM 0.80 and HIGH structural similarity *alone* lands a design
    there. Margins: rank 4 **+0.0076**, rank 17 +0.0217, rank 15 +0.0283, rank 18 +0.0331,
    rank 9 +0.0337, ranks 8 and 16 +0.0377. Our whole-chain `qtmscore` is an approximation
    of the organisers' domain-segmented computation (§9), so these verdicts are not robust
    to their pipeline or to a larger database. `cf_short120_r031` at TM 0.811 is this
    project's own precedent for landing on the wrong side by 0.011. **`bc_s360518_mpnn9_A22D` is the
    submission's sharpest eligibility exposure and it is 0.0076 from rejection.** The two
    antibody-format designs sit *above* 0.80 (0.855, 0.871) and clear anyway, on the CDRH3
    branch, where TM margin is not the cliff.
34. **The document numbers histidines in canonical coordinates and everything else in mature
    coordinates, and did not say so for its whole life.** `H433` is UniProt P00533 canonical;
    the same residue is **H409** in the mature numbering used by the domain boundaries this
    document states (III = 311–480), by every structure in `submissions/structures/`, by
    `targets/egfr/egfr_d3_6aru.pdb`, and by the contact footprints — which is why every
    design's footprint contains 409 and none contains 433. Mature 433 is not a histidine, so
    a reviewer checking our central claim against our own coordinates would find the wrong
    residue. The glycan sequon runs the other way: `Asn420` is mature (canonical 444). The
    mapping is now stated at §1; the ~110 in-place occurrences were not renamed, because
    doing that hours before a deadline is a larger risk than the mislabel. **Nothing
    numerical depends on it** — every pKa, distance and ratio was computed on the structures,
    in mature coordinates, and only the printed labels use the other convention.
35. **Two designs carry measured expression liabilities, and `bin/express_qc.py` had never
    been run on the shipped set.** It was run on 2026-10-05 against
    `submissions/01-egfr.csv` for the first time — `analysis/01-egfr/express_qc.tsv` had
    been a stale 20-design file, 10 rows of which are not in this submission and 8 of this
    submission's designs absent from it. 15 of 18 carry no flagged liability. The three
    that do:
    - **`c5_cf_short__boltzgen_egfr_cropfree_short_48`**, which the pH objective ranks
      **first**, is **30.0% alanine** (21 of 70) with a **7-residue poly-alanine run**, an
      **8-residue hydrophobic run** and **GRAVY +0.46** — the highest in the submission,
      where every other design sits at +0.24 or below. Binders are produced by *E. coli*
      cell-free synthesis, which is the regime in which long hydrophobic runs aggregate. This
      is a fifth independent strike against the top-ranked design, alongside its weak
      affinity, its Asn420 glycan contact, its lowest-in-set Chai ipTM of 0.201 and the
      unorderability of the objective itself.
    - **`rimA02_d3_rimA_14_vhh`** has **pI 6.38** (6.41 Bjellqvist, 6.54 EMBOSS) and a net
      charge of **−0.0 at pH 6.5**, so it sits at its isoelectric point *in the assay buffer
      where the headline low-pH measurement is taken*. A protein at its pI is at minimum
      solubility. It is the only design with this property: every other sits at pI 3.57–4.72
      and carries −5 to −15 net charge at 6.5. Nothing in the project's own QC flag logic
      tests pI, so this was invisible until measured directly under three independent pKa
      sets.
    - `h370_020_vhh` carries 2 cysteines, the canonical VHH framework pair, which is
      expected rather than a defect.
36. **The uploaded CSV contains no route to this document.** Eighteen rows defer their
    caveats to "METHODS 11.7", "METHODS 4.5" and similar, and the file carries no URL, no
    repository reference and no methods column — `grep -oE 'https?://[^ ,"]+'` on it returns
    nothing. Track 3 submits the CSV alone, so a grader reading only the uploaded artifact
    cannot reach any of the qualifications those rows rely on. The repository is public at
    `https://github.com/harishrajaram-svg/adaptyv-egfr-ph-switch` and the pointer is now
    included in each row's assessment text; before 2026-10-05 it was not.
37. **Two designs were removed by the organisers' own novelty check, and our gate had
    cleared both.** The submission was uploaded on 2026-10-05 with 18 designs. Proteinbase
    runs its novelty filter at upload and scored 16 at 3/4 and two at **2/4**, below the
    required bar, blocking submission until they were removed:
    - **`bc_s360518_mpnn9_A22D`** — our gate: Level 3 at qTM 0.7924, clearing the level-2
      cliff by **0.0076**. Limitation 33 named this as the submission's sharpest eligibility
      exposure and said our whole-chain `qtmscore` only approximates their domain-segmented
      computation. The margin was real and it fell the wrong way. It was the
      best-corroborated design in the submission — the largest causal swing (0.826× →
      5.630× against its own matched wild-type), the tightest seed reproducibility, a
      2.78 Å ASP22–H433 pair, and Chai-1 ipTM 0.788.
    - **`h370_020_vhh`** — our gate: Level **4**, on the antibody branch (CDRH3 identity
      0.273, global 0.526). Theirs: 2/4. Its whole-chain qTM is 0.8714, above the 0.80 HIGH
      line, so their pipeline does not appear to apply the antibody rule here — but
      `rimA02_d3_rimA_14_vhh` at qTM 0.8552 **passed** at 3/4, so it is not a simple
      general-rule substitution either. Two data points do not characterise the difference
      and we do not claim to understand it. **Our antibody-branch levelling should be read
      as unvalidated.**
    Removing both broke no declared parent/mutant pair: A22D's wild-type was never shipped
    and `h370_020_vhh` was not half of a pair. The submission is now **16 of 20**, with one
    antibody-format design remaining. The removals are recorded in
    `bin/emit_submission_csv.py` as an `INELIGIBLE` map rather than by hand-deleting rows,
    so the reason travels with the code and a re-emit cannot quietly reinstate them.
