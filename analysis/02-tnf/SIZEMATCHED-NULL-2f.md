# Item 2f: the designs scored against a null at their own length

2026-10-08. Two folds, same configuration as the 35 (ESMFold2-Full checkpoint, single sequence,
loops=10, steps=68, seed 42, L40S). Apps ap-snCV8ziLzEWYsls08s9Vh2 and ap-ihYnoYl5TTQ2UhM8vJA4Kp.

## Why this was needed

The control band committed in `p2_control_band_2026-10-08.tsv` shuffles the 164-residue TNFR2
ectodomain, so all five controls score at 302,382 cross-chain residue pairs. The designs are 76 or
84 residues and score at 219,486 / 227,022. `pae_interface_mean` is an average over those pairs,
so the band was not a yardstick for the designs. That gap was mine: I recommended that run without
checking size-matching first.

## The nulls

Chains A/B/C are the same 157-mer TNF-alpha protomers, read from the gated fixture. Chain D is a
composition-matched shuffle (seed 20261006, the convention `FIXTURES.json` records) of a design
named before scoring -- the first L76 and first L84 rows of the committed candidate TSV, picked by
file order so the choice cannot be read off the metric. `bin/build_sizematched_null.py` asserts the
resulting pair counts equal the designs' own.

| null | binder | pae_if_min | pae_if_mean | pairs |
|---|---|---|---|---|
| neg_shuffled_L76_1 | 76 | 0.854 | 14.248 | 219,486 |
| neg_shuffled_L84_1 | 84 | 0.855 | 14.996 | 227,022 |

## The designs against them

| | n | design mean (min-max) | median | null mean | beat the null |
|---|---|---|---|---|---|
| L76 | 21 | 11.218 - 16.586 | 14.047 | 14.248 | 15/21 |
| L84 | 14 | 13.945 - 15.434 | 14.648 | 14.996 | 10/14 |

25 of 35 beat their size-matched null on the mean (one-sided binomial p = 0.0083). The same 25 beat
it on the min, so the two metrics agree design-by-design at this configuration.

## What that is worth

Statistically above chance; biologically almost nothing. Scale the margins against what a REAL
receptor achieves over its own null -- TNFR2 scores 7.713 where its shuffles score ~18.99, a margin
of 11.28:

- best L76 design: 3.03 below its null = **27%** of the receptor's margin
- best L84 design: 1.05 below its null = **9%**
- median design: 0.20 (L76) / 0.35 (L84) = **2-3%**

So the median design is indistinguishable from a shuffled sequence of its own length and
composition, and the best one reaches about a quarter of the way to a real binder. This is the same
verdict Section 48 reached from the measured-design calibration, arrived at independently from
controls folded on our own configuration.

## Caveats, stated

1. **One draw per length.** The four TNFR2 shuffles had SD 0.08 on the mean, so a single draw is a
   defensible point estimate, but it is n=1 and carries no interval.
2. **The 27% figure assumes the null-relative gap is comparable across pair counts.** The size
   offset cancels to first order because each margin is measured against a null at its own size;
   it is not proven to cancel exactly.
3. This does not resolve item 2d. `pae_interface_min` was selected on AUC 0.901 from the MSA-fed
   arm, and our configuration has no MSA path. What 2f adds is that on OUR configuration the min
   and the mean rank the designs the same way relative to a size-matched null, which 2d could not
   determine.
