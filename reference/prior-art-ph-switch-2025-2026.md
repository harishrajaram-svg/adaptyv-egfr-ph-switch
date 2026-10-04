# Prior art on designed pH-switch binders — read 2026-10-04, ~21 h before the deadline

Two papers define the state of the art. Both matter for what we claim, and one of them
reframes the ceiling argument this project spent the day fighting.

## 1. Baker lab — "Computational design of pH-sensitive binders" (bioRxiv, 2025-09-29)
https://www.biorxiv.org/content/10.1101/2025.09.29.678932v1 | PMC12621740

**Direction: binders that RELEASE in acid.** pH 7.4 vs 5.4.
    EphA2    up to **1000x weaker** at 5.4
    TNFR2    122x weaker      TNFalpha  79x weaker
    IL-6     6x weaker        PCSK9     3.5x weaker
    Neo2     >2x weaker at pH 6.0

**Mechanism:** histidines on the BINDER, placed next to cationic residues (arginines) on the
target, so protonation at low pH creates electrostatic REPULSION. The 1000x EphA2 design put
His15/His19/His22 at 3.7 / 7.6 / 4.9 A from target arginines and showed **up to 11
histidine-cationic interactions**. A second strategy buries charged His networks in the
binder core to destabilise the fold in acid.

**They report NO design that binds MORE strongly at acidic pH.**

**Tools:** RFdiffusion -> ProteinMPNN with a histidine bias -> AF2 filtering -> Rosetta
sampling His in both protonation states. **No pKa predictor.**

**Hit rates:** interface design **4 pH-sensitive designs out of 12,000**. Buried-network
designs ~10-20% (6/43, 8/36, 7/40).

CAUTION: a search summary attributed a "10^8-fold" ILVBP result to this paper. That number
is NOT in the text. Do not cite it.

## 2. Jacobsen / Ovchinnikov et al — "pH-sensitive binder design with Proton-PottsMPNN"
(bioRxiv, 2026-09-30) https://www.biorxiv.org/content/10.64898/2026.09.30.755438v1

Represents protonated and deprotonated His, Asp and Glu as **distinct sequence tokens**, so
the design model can shape the microenvironment to favour one protonation state instead of
leaving pH sensitivity to be screened for afterwards.

8,407 de novo PD-L1 binders, yeast display across a pH series, **237 and 288 unique
pH-dependent designs recovered** (~3% each). **Transition pH 4.0-5.8.** Switches centred on
protonated **Asp and Glu as well as His**. No post-hoc structure-prediction filtering.

---

# WHAT THIS MEANS FOR OUR SUBMISSION

### 1. Our window is above everything published. Nobody has demonstrated a switch at 6.5/7.4.
Proton-PottsMPNN's transitions span **4.0-5.8**. Baker's headline comparison is 7.4 vs 5.4,
and the only result near our window is a **>2x at pH 6.0**. The competition asks for 6.5 vs
7.4 -- a **0.9 pH-unit** window sitting entirely above the published range. Our predicted
4.58x at 0.9 units is, if real, better than anything in either paper at a comparable gap.
That is a reason for humility about the prediction, not confidence about the design.

### 2. Our DIRECTION is the one neither paper achieved.
Both design release-in-acid. The Baker paper states plainly it has no design binding more
strongly at low pH. Our mechanism B -- a binder carboxylate raising a TARGET histidine's pKa
so the salt bridge only forms once the His protonates -- has the correct sign for this
challenge. That is a genuine point of distinction and belongs in the methods document.

### 3. THE CEILING WE FOUGHT ALL DAY IS A PROPERTY OF OUR MECHANISM CHOICE, NOT OF PHYSICS.
Our 7.94x-per-site bound is right, and over 0.9 units: 1 site 7.94x, 2 sites 63x, 3 sites 500x.
Baker reached 1000x over 2.0 units with **three or more** interface histidines (up to 11
His-cationic contacts). They multiplied sites; we spent the day trying to reach a SECOND site
and failed, because we chose mechanism B and **the target's histidines are wherever evolution
put them**. H370 is 8.5 A from H433 and unreachable; H358 is 26 A away.
**Binder-side titratable residues are placeable and multipliable. Target-side ones are not.**
That is the single most important strategic lesson available from this literature, and it says
our mechanism choice -- made before any of this was read -- capped the project at ~5.5x.

### 4. Both papers put the titratable residue on the BINDER. We put it on the target.
Mechanism A (binder His) was measured here as the weaker mechanism: 21/540 above the 1.20
bar versus a richer yield for B. But we never DESIGNED for mechanism A with multiple
histidines -- we only screened for it incidentally. Baker's hit rate for deliberately designed
interface His was 4/12,000, so incidental screening of 540 was never going to find it.

### 5. Hit rates are a warning about our numbers.
Baker: **4 experimentally pH-sensitive designs out of 12,000** for the interface strategy.
Ours are computational predictions with zero experimental validation, from ~2,000 designs.
Whatever our ranking says, the honest prior on any single design working experimentally is low.

### 6. Carboxylate-centred switching exists and we never considered it.
Proton-PottsMPNN recovered switches centred on protonated **Asp and Glu**, not only His. Our
entire project assumed a histidine is the titratable site in both mechanisms. A carboxylate
with an elevated pKa is a third mechanism we never tested -- and it is probably why their
transitions sit at 4.0-5.8, since a normal carboxylate pKa is ~4.
