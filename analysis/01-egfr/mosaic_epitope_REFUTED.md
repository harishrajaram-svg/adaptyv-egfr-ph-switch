# The 5/5 epitope result is explained by the crop, not by the objective

**Refuted 2026-10-03, same night, by the other window's conditional analysis.**
Commit 37bb380 said an optimized `BinderTargetContact` "looks like a different regime"
from BoltzGen's declared `binding_types` hint. **That claim is dead.** The correct
comparison is conditional on binding domain III at all — because a domain-III crop
removes the domain-I escape route by construction — and once conditioned, the two
methods are indistinguishable.

## The numbers that kill it

BoltzGen, across 8 arms, among designs that bound domain III:

| arm | n | domIII | patch given III | H433 given III |
|---|---|---|---|---|
| g532mimic | 20 | 12 | 8/12 | 6/12 |
| g532mimic_short | 20 | 13 | 12/13 | 10/13 |
| g532mimic_tiny | 20 | 7 | 6/7 | 6/7 |
| tinyHis | 20 | 7 | 7/7 | 7/7 |
| tight_h433_DD | 20 | 8 | 6/8 | 5/8 |
| tight_h433_EE | 20 | 7 | 7/7 | 7/7 |
| tight_h433_D1 | 20 | 5 | 4/5 | 3/5 |
| tight_h433 | 30 | 6 | 4/6 | 3/6 |
| **total** | | **65** | **54/65 = 83%** | **47/65 = 72%** |

Mosaic, cropped to domain III: 5/5 patch, 5/5 H433.
**Fisher exact two-sided p = 1.000.** P(5/5 | true rate 0.83) = **0.39** — a 5/5 run is
unremarkable at BoltzGen's own conditional rate.

## What this costs and what it buys

Costs: the headline. An optimized contact objective is **not** demonstrably better at
epitope targeting than a declared hint.

Buys something cheaper and more useful: **the crop is the fix, and it is free.** The
domain-I escape route that took roughly half of BoltzGen's designs disappears when the
target is cropped — a one-line spec change (`include: chain A` → a residue range), no
new tooling, no differentiable-geometry problem. That is the single most actionable
result of the night and it does not require Mosaic at all.

## What survives

The interface-spread diagnosis, for **both** methods. Mosaic: 22–49 target residues
contacted, 8–27% on patch. BoltzGen: ~24 target residues. Neither approach fixes it.

## Two corrections to my own framing

1. **I claimed geometry and interface confidence "pull apart" inside one objective**
   (s1 confident but missing by 0.96 Å; s2 clearing the bar with a soft interface).
   The other window's `rimA01_r15` breaks it: credible binding (ipSAE_min 0.5765, v2
   median human) **and** a 4.57× switch, dPKa +1.82, Ala knockout worth +2.72 pKa
   units. So the anti-correlation is strong but **not a law** — one objective can buy
   both, at least once. Do not conclude it cannot.
2. **The directional-term fix gains independent support.** `rimA01_r15` has ASP67
   carboxylate-O at **2.95 Å** from the H433 ring N with GLU92 at 4.76 Å — bidentate,
   and in the one real case the O is closer than the CA, matching the
   `acid_O < acid_CA` split that separated my near-misses from my failures. Prioritise
   the directional term over weight tuning.

## The experiment neither of us has run

`CROP=none` vs `CROP=284:453` at **matched seeds**. That isolates the objective from the
crop, and it is the only way to settle what the objective actually contributes.
