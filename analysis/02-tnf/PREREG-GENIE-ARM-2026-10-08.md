# Pre-registration: the Genie 3 arm's selection rule and its stop/go bar

Written **2026-10-08, 1:00 PM EDT, before a single SolubleMPNN sequence exists** and before any
complex from this arm has been co-folded. Committed before the numbers it governs.

## What is already known, so nothing here is blind

Honesty about non-blindness, in the manner of `analysis/02-tnf/calibrate_external.py`: the
**marginal** distributions of the three free inputs below are already measured and were in view
when this rule was written —

| input | measured, 2026-10-08 |
|---|---|
| monomer pLDDT, 100 seed sequences | min 0.656, median 0.824, max 0.897; 98/100 clear 0.70 |
| cation-reachable positions, 50 backbones | min 1, median 5.5, max 9 |
| interface contact pairs, 50 backbones | min 27, median 42, max 80 |

What is **not** known, and cannot be, is any co-fold outcome for this arm: no complex has been
folded, so no interface score, no margin over a null, and no linkage ratio exists for any Genie
design. The bar below is set against **previously published numbers from other work** and
**this project's own earlier designs**, never against this arm's output.

## The co-fold budget problem this rule exists to solve

The 0.70 monomer-foldability floor is the organisers' own pre-scoring filter
(`reference/anthropic-binder-design-protocol.md:100`). **0 of 35** of this project's Mosaic
designs clear it, so it eliminates all of them, against **2 of 100** Genie-derived sequences. At 98% survival it is a **floor, not a selector**: it cannot bound the co-fold spend,
which at ~90 s and ~$0.05 per 548-residue complex is the only expensive stage in the pipeline.

So the budget is **capped by rule**, and the rule is written here rather than chosen once the
scores are visible.

## Stage A — the 50 backbones already generated

**No ranking at all.** Take the **better-scoring monomer per backbone** by pLDDT, giving 50
candidates from 50 distinct backbones, and co-fold every one.

The reason is Genie3's own record: its 8 measured binders on this target are variants of **one**
backbone found in 23 samples, with pairwise identity 0.58-0.81 inside a set spanning 0.02-0.90.
The value sits in the **tail of backbone diversity**, not in the average of a family. Spending
A's budget on 50 distinct backbones rather than the top-scoring members of a few maximises the
chance of touching that tail, and 50 complexes is ~$2.50, below the point where selection pays
for itself.

A design is excluded from A only if it fails the 0.70 floor.

## Stage C — the rule if 1000 backbones are generated

1000 backbones x 2 seeds = 2000 sequences. Take the better monomer per backbone (1000), then
co-fold the **top 200** by this composite, computed on per-target z-scores:

    score = z(monomer_pLDDT) + z(cation_reachable_positions) + z(interface_contact_pairs)

Equal weights, and the reason is that **none of the three is a validated predictor of binding.**
The external calibration measured `plddt_binder` at AUC 0.568 (Boltz-2) to 0.766 (OpenFold3) --
a filter, not a ranker. The other two have no measured AUC at all. Weighting them differently
would be inventing a precision we have not earned, and the protocol's own 4:1 ipSAE:sc_DockQ
weighting is not transferable to different quantities.

**What this composite is, stated plainly:** a *necessary-conditions* budget allocator, not a
prediction. A design that is unfoldable, has no position within histidine reach of a target
cation, or has a vestigial interface cannot succeed even if it would otherwise bind. The
composite spends a fixed budget on designs that clear all three, and claims nothing more.

**Diversity constraint:** at most one sequence per backbone, as in A.

## The bar: what result justifies stage C

Each co-folded design is scored against a **size-matched shuffled null** -- a composition-matched
shuffle of its own 83-residue sequence against the same gated trimer, built by
`bin/build_sizematched_null.py`, so the comparison is at the design's own cross-chain pair count.
Primary metric `pae_interface_mean`, with `pae_interface_min` reported alongside; item 2f
established the two agree design-by-design at this configuration.

**Margin = null_mean - design_mean.** Positive means better than its own scramble.

The scale is fixed by numbers that already exist:

| reference | margin |
|---|---|
| **TNFR2, the real receptor**, over its own shuffles | **11.28** |
| this project's **best** Mosaic design | **3.03** (27% of the receptor's) |
| this project's **median** Mosaic design | 0.20 (L76) / 0.35 (L84) -- 2-3% |

**The bar is on the BEST design in A, not the median**, because the tail is what stage C buys:

- **ACTS, stage C is justified:** best margin **>= 5.64** (half the receptor's 11.28).
- **NULL, stop:** best margin **< 3.03** -- no better than what Mosaic already produced, so
  twenty times the sampling has no reason to help.
- **AMBIGUOUS, 3.03 to 5.64:** treated as **stop**. With under 67 hours to the retry-preserving
  upload and four METHODS items unwritten, an ambiguous result does not earn 16 GPU-hours and
  the attention that goes with them.

## What stage A cannot settle

A is a screen on the generator's **typical** output at n=50. It cannot prove anything about the
tail that stage C would buy. If A's median margin is poor while one design clears 5.64, that is
an argument FOR C, and the bar above is written to allow exactly that. The converse does not
hold: a best-of-50 below 3.03 is weak evidence about a best-of-1000, and the stop is a decision
about time and attention, not a proof that the arm cannot work. That distinction goes in METHODS.

## Downstream gates, in order, and what each may and may not claim

| # | gate | question | may claim |
|--:|---|---|---|
| 1 | 0.70 monomer floor | is it a protein? | nothing about binding |
| 2 | margin over size-matched null | does it bind at all? | the only binding evidence we have |
| 3 | `ph_gate_multisite.py --problem 2` | does the pH chemistry point the right way? | direction, not magnitude |
| 4 | PyRosetta `pH_mode` | how large with side chains repacked? | **only if it clears the control below** |

Gate 4 is a second opinion, never the pH screen. Its predecessor `bin/dddg_elec.py` passes a
**mechanism-absent far control 16/25 (64%)** against real anchors at 60% and 100% -- it cannot
separate mechanism from background. PyRosetta repacks where that file is rigid, which is the one
reason to expect better. **It must be run on the same far control first**, and if it also passes
at the anchors' rate it is abandoned and nothing downstream cites it.

## Standing caveat on every number this arm will produce

No structure predictor represents the pH-6.0 state. Gates 3 and 4 score protonation chemistry on
poses generated without any notion of pH. The pH objective remains a design-time geometric bet,
and nothing in this pipeline is a measurement of binding at either pH.
