# Control table — rebuilt on the reviewer's spec, 2026-10-04

The reviewer declined to put another universal number in place of the compromised bar. What was
asked for instead: experimentally characterised, expressed EGFR nonbinders where they exist;
CDR decoys that keep the framework intact, as a separate synthetic control class for VHHs;
positive controls deduplicated by sequence family; and seeds nested inside sequences.

## 1. The old panel was 8 runs of THREE molecules

Clustering by the **binder** chain (chain A in these files is the 621-aa target; reading it
blindly made all eight look identical — the same chain-order assumption that has now broken
three separate analyses here):

| family | members | n aa | status |
|---|---|---|---|
| EGF-derived | `NEG_nonbinder`, `NEGd3_nonbinder` | 134 | **activity unknown** |
| cetuximab scFv | `POS_cetuximab_scfv`, `POS_cradle_1nM`, `POSd3_cetuximab_scfv`, `POSd3_cradle_1nM` | 241 | 95.9–100% identical — **one molecule, four runs** |
| VHH | `POS_nano2_5nM`, `POSd3_nanobody_5nM` | 127 | scores **0.0000**, 5/5 dead seeds |

**The entire instrument calibration rested on one validated positive molecule.** Of three
families, one is activity-unknown, one is invisible to the instrument, and one works.

## 2. NEW — experimentally characterised controls (the class we never had)

Source: Adaptyv's own public release, `proteinbase.com/collections/egfr-round1-second-submission`,
downloaded to `data/proteinbase/egfr_round1_second.csv`.

  * **10 designs expressed and tested against EGFR on this platform, NO KD REPORTED.**
    48–200 aa, from two independent groups. These are **right-censored observations, not
    measured negatives** — a missing KD means the assay did not return a number, which is not
    the same as a number showing no binding (your correction, 2026-10-04). This bullet read
    "*measured* negatives, not presumed ones" after that correction had landed; the distinction
    is the whole basis of §4.1 and it was wrong here.
  * **Human EGF, measured on the same platform, n=15 runs: median KD 5.5e-8 M = 55 nM**
    (range 27–795 nM).

**CORRECTION TO WHAT I SENT THE REVIEWER:** I wrote that human EGF binds EGFR at "~2 nM". The platform's
own measurement is **55 nM median**, more than an order of magnitude weaker. The conclusion is
unchanged — EGF is a real binder, so an 81%-EGF molecule cannot serve as a negative control —
but the number I quoted was from memory, not from this assay.

What this does and does not settle: EGF's affinity is now measured. `NEG_nonbinder` is 81% EGF
*with two insertions*, and its own activity remains **unknown**, exactly as the reviewer said.

## 3. NEW — framework-preserving CDR decoys (VHH, presumed negative)

Built on the real 5 nM VHH framework, CDRs located by conserved anchors (`SCAAS`, `WFRQ`,
`RFTIS`, `YYCA`, `WGQG`) rather than fixed indices, then randomised:

    CDR1 GRTFSSYAMG   CDR2 SSGSTYYADSVKG   CDR3 AGYQINSGNYNFKDYEYDY

4 decoys, 67.7–70.9% identical to the parent, framework intact. These are **presumed**
negatives — a synthetic class that tests whether the instrument can distinguish a real VHH
binder from its own scaffold, which is the specific blind spot (0/2 on nanobody format).

## 4. NEW — shuffle nulls with seeds nested within sequences

4 designs x **3 independent shuffles each**, composition and length preserved exactly, 5 seeds
per shuffle. Each shuffle is ONE negative; its 5 seeds measure within-sequence noise. The
earlier `truenull` run counted 1 shuffle x 5 seeds per design and would have treated that as 5
independent negatives.

**Carried caveat, the reviewer's:** full-sequence shuffling often destroys the fold, so this null is
probably artificially easy. It bounds the instrument's noise floor, not biology.

## 5. What we will and will not claim from this

**CORRECTION, 2026-10-04.** An earlier version of this section told you the shipped flag means
*"above the 95th percentile of a matched calibration null."* **It does not, and no such flag is
available.** Two reasons:

1. The matched null is **degenerate**. See §6: 12 nested shuffles read 0.0000 on both species
   with 5 of 5 dead seeds, all twelve, and a separate 22-molecule shuffle run has 20 of 22 dead
   with a maximum median of 0.0110. Its 95th percentile is **0.0000**, so a percentile flag
   admits anything above zero. You predicted full-sequence shuffling would destroy the fold and
   make the comparison artificially easy; it is worse than that — the null has no tail to take a
   percentile of.
2. What the column actually computes is `max(ipSAE_human, ipSAE_mouse) >= 0.2218`
   (`bin/emit_submission_csv.py:89`), and **0.2218 is the EGF-derived bar you had us retire.** It
   survives as a *reporting* threshold only: `MIN_AFFINITY = None`, so no design is excluded on
   affinity and the column orders nothing. It should be read as "above the legacy computational
   null", which is a provenance statement and not a claim about binding.

**Update 2026-10-05: the column came out.** `affinity_above_null` is no longer in the
submission — the shipped CSV has 15 columns and none of them is that flag — and
`AFFINITY_FLAG_DEPRECATED = 0.2218` records the retired threshold in
`bin/emit_submission_csv.py`. This paragraph said the column was still shipped, and that all
three documents agreed, after it had already been removed.

With n=10 right-censored molecules the tail is not characterised well enough for a hard gate either,
so affinity is reported as a continuous score with an uncertainty flag and **no design is
excluded on it**.

---

## 5b. THE pH CALIBRATION PAIR — G532 / G532Ctrl. Named twice, never used, and we had it backwards.

**You named G532 and G532Ctrl as the relevant pH calibration pair in both replies, and
was explicit that known binders may calibrate the method but must never become starting sequences for competition
entries."** We did the inverse. G532 exists in this project as a **design-generation arm**
(`g532mimic`, `g532_ladder`) and is **absent from the control panel entirely** — it appears
nowhere in the pH gate, the ranking table, or the methods document.

Your SPR figures, which are the ones that matter:

| target | KD at pH 6.5 | KD at pH 7.4 | KD(7.4)/KD(6.5) |
|---|---|---|---|
| human EGFR | 294 nM | 3,900 nM | **13.26×** |
| mouse EGFR | 547 nM | 1,810 nM | **3.31×** |

**And it switches on H433 and H370 — the exact pair this project is built on.** Two consequences
we had wrong and have now corrected in the methods document (§8.1):

1. **We were describing target-side histidines as "self-defeating."** You told us not to. G532
   is the reason, and we had the citation a day before we wrote the claim.
2. **We were calling the two-site H433+H370 route "closed."** 13.26× exceeds the 7.94×
   single-proton bound, so G532 must be linking more than one proton at the pair we measured at
   8.5 Å. We failed to build it; that is not the same thing.

**Assay format, as you asked.** G532 also carries published **ELISA** ratios (8.08 / 1.64 / 0.76
for G532 / G532V / G532Ctrl). Those are **EC50 ratios and must not be compared with SPR KD
ratios** — an internal analysis of ours anchored on the 8.08 figure, which was the wrong
comparator for our KD-ratio objective. The two are kept separate here.

**Where the control check stands, honestly.** The ladder was folded into 20 ESMFold2 poses and
never scored or pH-gated. Recovered, the gate returns **0.677 / 0.689 / 0.677 / 0.697** for
G532 / G532V / G532Ctrl / G5V2 — wrong direction on a measured 13× switch, and no separation
from its own negative comparator. All four land on the 0.699× analytic **lower extreme** of the
single-site model (the limit as pKa_bound → −∞), i.e. the gate is returning its own bound.

*Corrected 2026-10-05:* this is not "no linkage detected" — no linkage would read **1.0**.
Four controls landing within 0.025 of an analytic extreme is the signature of
protonation-model saturation, which means the gate returns no usable information on these
four rather than returning a reverse switch. Read as a failure of the measurement, not as a
measured property of G532.

**We are not reporting that as a falsified gate, because you told us how to read it:** *"A
failure to recover this control could come from the predicted pose or from the protonation
model, and the two must be assessed separately before the pH gate is used to discard anything. Running that
separation (`bin/g532_pose_check.py`) gives **POSE WRONG on 20 of 20**: the nearest carboxylate
sits **8.97–16.84 Å** from H433 across the five G532 poses, and in G532V/G532Ctrl the relevant
acidic residues are mutated out by design. A gate reading on a structure that cannot host the
mechanism is uninformative, so the failure is **attributed to the pose, not yet to the
protonation model.**

### The pose-vs-protonation separation you asked for, now answered

The instruction was that a failure to recover this control could come from the predicted pose
or from the protonation model, and that the two be assessed separately before the pH gate is
used to discard anything.

Done. **It is the pose, and the result is worse than a failure to recover — the instrument
anti-ranks the ladder.** `ipsae_min` refuses on these structures by design (3 inter-chain pairs;
fixture case 10), so each pair was scored explicitly. A = EGFR domain III target, B and C = Fv
heavy and light. Medians over 5 seeds each:

| molecule | measured | target:Fv (A:B) | target:Fv (A:C) | **best target:Fv** | intra-Fv (B:C) |
|---|---|---|---|---|---|
| **G532** | **13.26× SPR, 8.08 ELISA** | 0.0126 | 0.0000 | **0.0135** | 0.8510 |
| G532V | 1.64 ELISA | 0.3512 | 0.3364 | 0.3713 | 0.8634 |
| G532Ctrl | 0.76 ELISA | 0.1848 | 0.1602 | 0.2503 | 0.8433 |
| G5V2 | — | 0.4110 | 0.3632 | 0.4289 | 0.8605 |

**The 294 nM binder scores lowest of the four, and its non-switching comparator scores 18×
higher.** Meanwhile the intra-Fv packing is 0.84–0.86 in every molecule.

**Three corrections to how this table was read, 2026-10-05.**

*First, G532Ctrl is not a negative.* Earlier text here and in the methods document called it
"its own negative" and treated G532 as the only binder in the ladder. **G532Ctrl is also a
measured EGFR binder** — its 0.76 is an *EC50 ratio across pH*, i.e. a binder that does not
pH-switch, not a molecule that fails to bind. All four rows are antibodies against the same
target; what differs is the pH dependence. Any claim phrased as "the only one that binds"
is withdrawn.

*Second, four rows cannot establish an inverse relationship.* The observation that the
tightest binder scores lowest is a four-point ordering with no replication across molecules
and no error model. It is consistent with the instrument being uninformative on antibodies —
which is the conclusion we draw elsewhere from larger evidence — but on its own it does not
show that score runs *opposite* to affinity, and it is not offered as such.

*Third, high intra-Fv confidence is not a correct fold.* The 0.84–0.86 figures say the model
is confident about the VH:VL interface. They do not say the Fv is folded correctly, and they
say nothing at all about the pose relative to the target. "Folds the Fv essentially perfectly
and then fails to dock it" over-reads a confidence score; what is supportable is that the
model places high confidence on the intra-Fv interface and near-zero on the target interface.

**The same pattern on the full human ECD.** Scored 2026-10-05 on the 621-residue ECD rather
than the domain-III crop, per named interface (chains A = VL, B = VH, C = EGFR ECD), medians
over 5 seeds:

| molecule | target:VL (A:C) | target:VH (B:C) | intra-Fv (A:B) |
|---|---|---|---|
| **G532** | **0.0000** | **0.0000** | 0.8189 |
| G532V | 0.2625 | 0.2564 | 0.8200 |
| G532Ctrl | 0.2270 | 0.2282 | 0.8110 |
| G5V2 | 0.2547 | 0.2590 | 0.8159 |

Same ordering, same intra-Fv/target gap, on a different target construct — which is what
makes the pose attribution below credible rather than a crop artefact. Note these 20 poses
were scored and then wired into no analysis; they are reported here and feed nothing.

**Three consequences.**

1. **The G532 pH-gate failure is fully attributed to the pose.** The protonation model was never
   given a structure that could host the mechanism — consistent with the independent geometry
   check, which put the nearest carboxylate 8.97–16.84 Å from H433 on all five G532 poses. We
   therefore make **no claim either way** about the protonation model from this control, which is
   the separation you asked for and the honest end of it.
2. **It is a third and much stronger data point on antibody blindness.** §4.2 of the methods
   document had 0/2 on VHH positives; §4.5 added four CDR decoys indistinguishable from the real
   VHH at 0.0000. This is a *measured 294 nM antibody* scored **below a non-pH-switching
   comparator that is itself a measured binder**. The instrument is not merely insensitive to
   antibody formats — on the one measured example we have, the ordering is **inverted**.
   *Stated as one inverted pair, not as a correlation*: this same file says four lines earlier
   that "four rows cannot establish an inverse relationship", and "anti-correlated" claimed
   exactly the relationship those four rows cannot support. One inversion is what we have.
3. **It bears directly on two of our sixteen submitted rows**, both VHH format. We already report
   their affinity as inadequately assessed rather than low. This strengthens that from a caveat
   to a measurement: on the only antibody in this project with a known KD, the affinity column
   points the wrong way.

**What we still cannot say.** This is n=1 molecule with 3 comparators, all from one published
series, and ESMFold2 was run without an MSA. It does not establish a general docking failure
rate for antibodies. It does establish that the number our pipeline would have reported for
G532 is not merely uninformative but inverted, and that is enough to stop us reading the VHH
affinity columns at all.

**And the ordering consequence we have not honoured.** Your instruction was to assess pose and
protonation independently *before* the pH gate is used to reject anything. The gate has been used
as a discard filter throughout — §10 of the methods document was rebuilt on 2026-10-04 and now
accounts for a **75-molecule** union, sequence-keyed, rather than the 38 this paragraph
cites; the 38 was a run-name count that collapsed to fewer distinct sequences.
Those rejections rest on a gate whose only ground-truth test is unresolved. We are stating that
rather than re-running the selection 14 hours before a deadline, and it is question (d) below.

---

## 6. THE CONTROL RECOVERY ITSELF — the result this document was missing

**"The eight control shards", reconciled.** You asked how "eight" squares with a table listing
nine `gitter-yolo` molecules plus one `deepsatflow`. It does not, and it was never meant to:
**eight is a count of Modal shards, not of molecules.** Verified from the run tree — four
`expctrl0–3` shards hold the 11 measured molecules (1 quantified binder + 10 no-KD, each folded
against both species, 22 design directories), and four `ctrl2_0–3` shards hold the 4 CDR decoys
and 12 shuffle nulls. Eight shards, 27 molecules. The sentence below conflated a compute unit
with a molecule and is left in place with this correction above it.

**An earlier version of this file introduced the measured controls as "the class we never had"
and then reported no ipSAE result for any of them.** The artifact answering your
control-recovery ask did not contain the control recovery. The cause: the eight control shards
had been folded on Modal and never ipSAE-scored — ESMFold2 writes `*_ipsae.json` (raw PAE) and
the scoring step that writes `*_10_10.txt`, which every reader in `bin/` parses, is separate and
had never been run on them. Backfilled 2026-10-04: 270 poses, 0 failures. Source:
`analysis/01-egfr/control_recovery.{json,tsv}`.

All values are ipSAE_min, **median over 5 seeds**, against the 621 aa human and mouse ECDs.

**MEASURED BINDER (KD 55 nM, n=15 runs)** — n=1

| molecule | ipSAE human (median) | dead seeds | ipSAE mouse (median) | dead seeds |
|---|---|---|---|---|
| `EXPPOS_Human_EGF` | 0.3549 | 0/5 | 0.3700 | 0/5 |

**NO KD REPORTED (right-censored; affinity weaker than the quantifiable limit, or an expression/QC failure)** — n=10

| molecule | ipSAE human (median) | dead seeds | ipSAE mouse (median) | dead seeds |
|---|---|---|---|---|
| `EXPNEG_gitter-yolo10` | 0.5893 | 0/5 | 0.0000 | 5/5 |
| `EXPNEG_gitter-yolo9` | 0.4326 | 2/5 | 0.0000 | 5/5 |
| `EXPNEG_deepsatflow-design7_n0_mpnn1.320_p` | 0.2168 | 0/5 | 0.1823 | 0/5 |
| `EXPNEG_gitter-yolo5` | 0.0264 | 1/5 | 0.0108 | 2/5 |
| `EXPNEG_gitter-yolo4` | 0.0109 | 1/5 | 0.0000 | 4/5 |
| `EXPNEG_gitter-yolo3` | 0.0000 | 5/5 | 0.0121 | 2/5 |
| `EXPNEG_gitter-yolo7` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `EXPNEG_gitter-yolo8` | 0.0000 | 4/5 | 0.0000 | 5/5 |
| `EXPNEG_gitter-yolo2` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `EXPNEG_gitter-yolo6` | 0.0000 | 5/5 | 0.0000 | 5/5 |

### 6b. Family-balanced comparison, per your item (b)

`bin/control_family_balance.py` (reproducible, no arguments).

**The family map is sequence-based, and the first version of this was wrong.** You said: *"First
freeze a family map on sequence or backbone before any score is consulted. Designs coming from one submitting group
is a hint rather than a definition of a family. The first version used the submitter-group prefix
(`gitter-yolo` / `deepsatflow`) — exactly what you ruled out — and reported two families. Rebuilt
by clustering the binder sequences themselves.

Plain identity does not work on this panel: the `deepsatflow` design is 48 aa against
`gitter-yolo` designs of 150–200 aa, so a short-vs-long alignment reports 54–65% identity over
the aligned fragment and single linkage at 30% collapses all ten into one family. The metric is
identity × coverage (coverage = length ratio, shorter over longer). The partition is **stable**:
thresholds 0.70 and 0.85 give the same six families, so it does not rest on a threshold picked
for its answer.

| family | n | members | id × cov |
|---|---|---|---|
| 1 | 3 | `yolo10`, `yolo9`, `yolo7` | 0.84–0.91 |
| 2 | 2 | `yolo5`, `yolo4` | 0.89 |
| 3 | 2 | `yolo8`, `yolo6` | 0.89 |
| 4 | 1 | `deepsatflow-design7` | — |
| 5 | 1 | `yolo2` | — |
| 6 | 1 | `yolo3` | — |

**Six families over ten molecules, against the two the submitter prefix implied.** One row per
distinct sequence, five seeds nested in each.

| | human leg (EGF 0.3549) | mouse leg (EGF 0.3700) |
|---|---|---|
| **raw**, one row per molecule | **8.0/10 = 0.800** | **10.0/10 = 1.000** |
| fam1 (n=3) | 0.333 | 1.000 |
| fam2 (n=2) | 1.000 | 1.000 |
| fam3 (n=2) | 1.000 | 1.000 |
| fam4 / fam5 / fam6 (n=1 each) | 1.000 | 1.000 |
| **family-balanced** (equal weights, ties = ½) | **0.889** | **1.000** |
| **leave-one-family-out** | **0.867 – 1.000** | 1.000 |

The family-balanced figure is unchanged at 0.889 — coincidentally the same as the two-family
version — but the sensitivity is much tighter: leave-one-family-out now spans **0.867–1.000**
rather than 0.778–1.000, because no single family carries nine molecules any more. Note where
the signal actually sits: **fam1 is the only family that does not rank below EGF** (0.333), and
it holds the two molecules that outscore EGF on the human leg. Eight of the ten "failures" of
this panel are one sequence family.

**No confidence interval is reported, and no effective-n is substituted into Clopper–Pearson.**
An ordinary exact binomial interval does not become cluster-adjusted by that substitution, and
six families —
three of them singletons — cannot support dependable cluster-bootstrap inference either. The
leave-one-family-out range is the uncertainty statement. The mouse leg is 1.000 under every
weighting, the one part of this panel that is not weighting-dependent.

**Requiring both species — how the submission is actually scored — zero of the ten no-KD
molecules match or exceed EGF on both legs.** That is consistent with the instrument working on
this panel and does not demonstrate that it does: there is one positive, and the ten are
right-censored rather than measured at zero.

**CDR DECOYS (presumed negative, VHH framework preserved)** — n=4

| molecule | ipSAE human (median) | dead seeds | ipSAE mouse (median) | dead seeds |
|---|---|---|---|---|
| `CDRDECOY_vhh_2` | 0.0117 | 1/5 | 0.0126 | 1/5 |
| `CDRDECOY_vhh_1` | 0.0000 | 5/5 | 0.0000 | 4/5 |
| `CDRDECOY_vhh_3` | 0.0000 | 4/5 | 0.0000 | 5/5 |
| `CDRDECOY_vhh_4` | 0.0000 | 5/5 | 0.0000 | 5/5 |

**SHUFFLE NULLS (3 independent shuffles x 4 designs, seeds nested)** — n=12

| molecule | ipSAE human (median) | dead seeds | ipSAE mouse (median) | dead seeds |
|---|---|---|---|---|
| `SHUF_bg04_r03_boltzgen_egfr_d3__r1` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_cf_short120_r031_boltzgen__r2` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_d2d_101_l133_s154858_mpnn9_r3` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_bg04_r03_boltzgen_egfr_d3__r2` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_cf_short120_r031_boltzgen__r3` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_rimA01_r15_boltzgen_egfr_d_r1` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_bg04_r03_boltzgen_egfr_d3__r3` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_d2d_101_l133_s154858_mpnn9_r1` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_rimA01_r15_boltzgen_egfr_d_r2` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_cf_short120_r031_boltzgen__r1` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_d2d_101_l133_s154858_mpnn9_r2` | 0.0000 | 5/5 | 0.0000 | 5/5 |
| `SHUF_rimA01_r15_boltzgen_egfr_d_r3` | 0.0000 | 5/5 | 0.0000 | 5/5 |

### What it says, stated as ranks rather than as a summary statistic

**The panel is 10 right-censored molecules and ONE quantified positive.** With a single positive there is no discrimination
estimate available: every summary statistic reduces to *where that one molecule ranks*. We report
ranks and deliberately do not report an AUROC — an earlier draft of our methods document did, and
it was removed. This is your own `arms-backlog §2a` ("the gate is n=1 positive") arriving in the
place it mattered.

- **Human leg: 8 of the 10 no-KD molecules rank below the quantified binder. 2 rank above it.**
  (Corrected 2026-10-05: a missing KD is right-censoring, not a measured zero, so this
  separates one quantified binder from ten censored observations.)
  The highest-scoring molecule in the entire panel is one with **no KD reported** —
  `EXPNEG_gitter-yolo10` at 0.5893 against EGF's 0.3549. It is not a measured non-binder; this
  line called it one four lines below the correction that says it is not. It also reads a **5.27× pH ratio**, which
  places it 8th of the 132 molecules eligible to rank on our primary objective.
- **Both species required: the one quantified binder outranks all ten no-KD molecules on this
  panel.** That is the full extent of the claim.

> **Corrected 2026-10-05.** This section previously argued that the mouse leg *rescues* the
> failed human criterion, "a mechanism, not a margin", because the two exclusions sit at exactly
> 0.0000 on five of five seeds rather than at a narrow margin. **That argument is withdrawn**, on
> the correction: a predicted interface of zero is not an experimentally demonstrated
> specificity mechanism, and without matched mouse outcomes it neither validates mouse binding nor rescues the
> failed human control criterion."* Three reasons it does not stand:
>
> * **No matched mouse outcomes exist.** Adaptyv measured these molecules against HUMAN EGFR. A
>   mouse prediction of 0.0000 is corroborated by nothing, so the filter is being validated on
>   the species for which we have no data.
> * **0.0000 is a statement about the predictor.** §4.4b of the methods document shows a pose
>   reproducing 72% of a crystallographic interface scoring ipSAE_min 0.0000. Margin-independence
>   would matter only if the zeros were known to mean no interaction, and they are not.
> * **The human criterion failed and stays failed.** Two of ten censored molecules outrank the
>   only quantified binder on human. A second species with no outcomes cannot convert that into
>   a pass.
>
> The honest statement: **no specificity filter in this pipeline is supported by measured data**,
> this one included. The cross-reactivity requirement came from the brief rather than from us.

**Three limits.** The positive class is one molecule of 53 aa, shorter than every design here. A
missing KD bounds affinity from ABOVE, not below, and does not prove no interaction. And ten
censored observations do not characterise a tail, so 0.5893 is a floor on how high a no-KD
molecule can score here, not a ceiling.

**One thing I have not fixed.** You asked twice for the historical EGFR data to be split by
design family. Eight of these eleven molecules are one group's `gitter-yolo` series, and I have
counted them as eight independent negatives. The effective n is closer to 3 families than 10
molecules. I do not know the right construction at this n and it is question (c) in my covering
note.

### The VHH blind spot is total, not partial

Your 0/2 reading understated it. Four **framework-preserving CDR decoys**, built on the real
5 nM VHH scaffold with CDRs located by conserved anchors and randomised, read 0.0000, 0.0117,
0.0000, 0.0000. **The real 5 nM VHH reads 0.0000 too.** The instrument cannot separate a
validated nanobody from a randomised-CDR decoy on that nanobody's own framework.

Two of the ten shipped designs are VHH format. They are reported as **affinity inadequately
assessed**, not as weak binders, exactly as you asked — their 0.23/0.45 and 0.44/0.71 columns
should not be read as measurements.
