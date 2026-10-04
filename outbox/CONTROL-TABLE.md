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

A provisional screening flag only: *above the 95th percentile of a matched calibration null*.
That means "above this computational null", **not** "binds EGFR", and it does not establish a
5% biological false-positive rate. With n=10 measured negatives the tail is not characterised
well enough for a hard gate, so affinity is reported as a continuous score with an uncertainty
flag and **no design is excluded on it**.
