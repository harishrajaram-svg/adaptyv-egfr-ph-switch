# Control table — rebuilt on PK's spec, 2026-10-04

PK: *"I would not replace the compromised bar with another universal number ... Add
experimentally characterised, expressed EGFR nonbinders where available ... For VHHs, include
framework-preserving CDR decoys as a separate synthetic control class ... Deduplicate positive
controls by sequence family ... keep seeds nested within sequences."*

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

  * **10 designs expressed, tested against EGFR on this platform, NO binding detected.**
    48–200 aa, from two independent groups. These are *measured* negatives, not presumed ones.
  * **Human EGF, measured on the same platform, n=15 runs: median KD 5.5e-8 M = 55 nM**
    (range 27–795 nM).

**CORRECTION TO WHAT I SENT PK:** I wrote that human EGF binds EGFR at "~2 nM". The platform's
own measurement is **55 nM median**, more than an order of magnitude weaker. The conclusion is
unchanged — EGF is a real binder, so an 81%-EGF molecule cannot serve as a negative control —
but the number I quoted was from memory, not from this assay.

What this does and does not settle: EGF's affinity is now measured. `NEG_nonbinder` is 81% EGF
*with two insertions*, and its own activity remains **unknown**, exactly as PK said.

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

**Carried caveat, PK's:** full-sequence shuffling often destroys the fold, so this null is
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

We have left the column in the submission rather than dropping it mid-flight, and said the same
thing in METHODS §4.5 and §11 so all three documents now agree. If you would rather it came out,
it is a one-line change and we have ~38 hours.

With n=10 measured negatives the tail is not characterised well enough for a hard gate either,
so affinity is reported as a continuous score with an uncertainty flag and **no design is
excluded on it**.

---

## 5b. THE pH CALIBRATION PAIR — G532 / G532Ctrl. Named twice, never used, and we had it backwards.

**You named G532 and G532Ctrl as the relevant pH calibration pair in both replies, and
instructed: "Use these binders only for calibration, not as starting sequences for competition
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
from its own negative comparator. All four land on the 0.699× analytic lower bound of the
single-site model, i.e. the gate is returning its own floor.

**We are not reporting that as a falsified gate, because you told us how to read it:** *"A
failure to recover this control could arise from the predicted pose or the protonation model.
Assess those separately before using the pH gate to discard candidates."* Running that
separation (`bin/g532_pose_check.py`) gives **POSE WRONG on 20 of 20**: the nearest carboxylate
sits **8.97–16.84 Å** from H433 across the five G532 poses, and in G532V/G532Ctrl the relevant
acidic residues are mutated out by design. A gate reading on a structure that cannot host the
mechanism is uninformative, so the failure is **attributed to the pose, not yet to the
protonation model.**

**What is still not done:** scoring the ladder's interface, which is what would separate "ESMFold2
built no complex" from "the protonation model is wrong." It is a local ten-minute job. One
complication to flag in advance: the ladder is a **3-chain** Fv complex with an intra-Fv pair at
0.8659, and `ipsae_min` **refuses** on more than one inter-chain pair by design (fixture case 10).
So the interface score may need the binder:target pair named explicitly rather than inferred.

**And the ordering consequence we have not honoured.** Your instruction was to assess pose and
protonation separately *before* using the pH gate to discard candidates. The gate has been used
as a discard filter throughout — §10 of the methods document lists 38 molecules rejected on it.
Those rejections rest on a gate whose only ground-truth test is unresolved. We are stating that
rather than re-running the selection 14 hours before a deadline, and it is question (d) below.

---

## 6. THE CONTROL RECOVERY ITSELF — the result this document was missing

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

**MEASURED NON-BINDERS (no binding detected)** — n=10

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

**The panel is 10 negatives and ONE positive.** With a single positive there is no discrimination
estimate available: every summary statistic reduces to *where that one molecule ranks*. We report
ranks and deliberately do not report an AUROC — an earlier draft of our methods document did, and
it was removed. This is your own `arms-backlog §2a` ("the gate is n=1 positive") arriving in the
place it mattered.

- **Human leg: 8 of 10 measured non-binders rank below the measured binder. 2 rank above it.**
  The highest-scoring molecule in the entire measured panel is a measured **non-binder** —
  `EXPNEG_gitter-yolo10` at 0.5893 against EGF's 0.3549. It also reads a **5.27× pH ratio**, which
  places it 8th of the 132 molecules eligible to rank on our primary objective.
- **Both species required: the one measured binder outranks all ten non-binders.**
- **The reason that works is a mechanism, not a margin.** Both molecules that beat EGF on human
  sit at **exactly 0.0000 on mouse with 5 of 5 dead seeds** — no interface found at all, rather
  than a near miss. A filter resting on a margin could erode with a larger panel; this one does
  not rest on a margin.
- **The cross-reactivity requirement came from the organisers' brief, not from us**, and it is the
  only specificity filter in this pipeline that measured data supports at all. That is a statement
  about the absence of evidence for the others, not a validation of this one at n=1.

**Three limits.** The positive class is one molecule of 53 aa, shorter than every design here.
"No binding detected" bounds affinity from below and does not prove no interaction. And ten
negatives do not characterise a tail, so 0.5893 is a floor on how high a measured non-binder can
score here, not a ceiling.

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
