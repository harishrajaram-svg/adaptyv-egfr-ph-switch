# Conditional TNF-α binder — methods

**Anthropic × Adaptyv Protein Design Competition, Challenge 2. Track 3.**
Repository: <https://github.com/harishrajaram-svg/adaptyv-egfr-ph-switch> (problem 1; problem 2
tree to follow)

> **STATUS: SKELETON, 2026-10-06.** Sections 1–6 are complete and will not change — they report
> measurements already made. Sections 7–10 depend on designs that do not exist yet; each states
> explicitly what it is waiting on. Nothing in this file is a placeholder number: where a count
> belongs it is generated from an artifact or the section says "pending".

---

## 1. Summary

We designed against the receptor-binding groove of human TNF-α with binder-side histidines
intended to release the target at pH 6.0, and we report a **negative methodological result that
is stronger than our design result**: the two published computational filters for this exact
objective, implemented faithfully and tested against the largest labelled dataset that exists,
**do not work on the class of design this competition asks for.**

Our submission is therefore ranked on quantities we can validate — predicted interface quality
and cross-species binding — with the pH mechanism **designed into the generator and reported as
unscorable** rather than claimed.

**Three results we did not expect going in:**

1. `dddG_elec` (Ahn et al., Baker lab) scores **41/82 = 50%** on single-point histidine variants
   with measured labels, against an **83% trivial null**. §5.
2. `his_cation_gate` — Ahn's second, geometric criterion — fires on **0 of 14** measured switches
   while a slack control finds a cation within 8 Å in **66 of 82** cases. §5.
3. The mechanism both are built to detect has **no measured precedent**: across ~100 real
   molecules it has never been observed in a confirmed pH switch. §5.4.

---

## 2. Target

Human TNF-α, UniProt **P01375** residues 77–233, assayed as the soluble trimer.

🔴 **The deposited structure is wrong at one position and it is in the epitope.** PDB **1TNF**
carries **Leu** where canonical TNF-α has **Asp219** — a hydrophobic→anionic change in the
consensus receptor core, present in 10 of 10 receptor copies. Reported by Ken Osumi in the
organisers' channel; verified independently. We rebuilt the sidechain with PDBFixer
(`targets/tnf/tnf_canonical_trimer.pdb`, 0 mismatches at all 152 observed positions per chain).

**Three of our own guards had the opportunity to catch this and all three declined**, each
because it verified identity only where we were already looking. The guard we then built
(`analysis/02-tnf/target_identity.py`) checks **every position of every consumed structure**, and
widening it from generation targets to analysis fixtures immediately found three more:

| structure | divergence | consequence |
|---|---|---|
| 3WD5 (+5 derived) | R107D, recorded by the PDB as a *conflict* | measured 15.82 Å from the nearest scored histidine — outside the 5.5 Å truncation, contributes exactly zero |
| 3IT8 | D219L | ledger-only entry; no number moves |
| **3ALQ** | **six mutations: K87M K141S K166P K174R K188N K204P** | **a six-lysine mutein, and K166 — a primary anchor — is PROLINE in it** |

**3ALQ is four of the ten receptor copies our conservation claims rest on.** The epitope itself
survives, because every consensus position is also in the contact set of the six wild-type TNFR1
copies — but **"present in 10 of 10 receptor copies" conflates contact position with residue
identity** and is restated throughout as **6 of 10 wild-type copies plus 4 in a mutein**.

Mouse (**P06804**, 2TNF) was checked for the first time and is **clean: 148 of 148 observed
positions identical**. Its numbering is human-aligned with a deletion at 73, so no single offset
fits it and a naive offset mislabels every residue past 72.

---

## 3. Epitope, and the three-protein constraint

The consensus receptor epitope (21 residues, present in all TNFR2 and TNFR1 copies) straddles two
protomers; a binder cannot sit on one chain.

**Lymphotoxin-α (P01374) binds the same groove**, so the design must hold **two** proteins —
human and mouse TNF — and release a **third**. Classifying all 21 positions on both axes at once:

| class | positions | n |
|---|---|---|
| **safe + selective** (conserved in mouse, divergent in LT-α) | 97, **108**, 162, 163, 167, 191, 220, 221, 222, 225 | **10** |
| safe, not selective | 109, 151, 153, **166**, 189, 219 | 6 |
| hazard (mouse differs) | 96, 107, 149, 161, 173 | 5 |

🔴 **R108 is the only cationic anchor that is conserved in mouse AND divergent in LT-α.** K166 —
the anchor every reachability and separation measurement favoured — is **identical in LT-α** and
therefore buys no selectivity. The easy anchor is the non-selective one.

✅ **The fold is conserved.** Superposing the real mouse trimer on canonical human gives an
epitope **Cα RMSD of 1.10 Å over all 21 positions**, so cross-reactivity is a sidechain-chemistry
problem, not a geometry problem. Two positions require ≈ **+40 Å³** of headroom for mouse — **96
(P→H) and 149 (H→Y)** — not one.

---

## 4. The designed mechanism, and what it is priced at

Binder-side histidines placed against the conserved cations **R108** and **K166**: protonation at
acid creates like-charge repulsion and releases the target. This is Ahn et al.'s mechanism (i).

**Thermodynamic ceiling.** The 7.4 → 6.0 window gives a one-proton ceiling of **25.1×**. Two
sites at a 2.0-unit pKa drop give 51.6×; at Broo et al.'s **measured** 0.5 units per flanking
cation they give **4.5×**. We report the **band 4.5× – 51.6×** and never an endpoint. ⚠️ A PubMed
search for a measured ΔpKa of a histidine in direct ~3 Å contact with an Arg or Lys across an
interface **returns nothing** — the geometry has no experimental calibration in either direction.

🔴 **And no published molecule meets objective 1 as specified.** Against a 10 µM quantification
ceiling, "no detectable binding at pH 6.0" requires ratio > ceiling/KD(7.4). The closest molecule
in the literature — Schröter PSV#2, KD 77 pM with a 785× off-rate shift — **would need to be 165×
better**. Adalimumab would need to be 241,000× better. **The escape is to bind weakly on
purpose**: at 500 nM the spec needs 20×, at 1 µM it needs 10×, and both sit inside the published
5–30× band. Objective 1 is reachable only where objective 3 is deliberately sacrificed, and that
tension is in the spec's arithmetic rather than ours.

---

## 5. 🔴 The instrument result — the strongest finding here

We implemented both published filters for this objective and tested them against
**Wei & Sulea's benchmark** (PMID 39289783): **279 single-point histidine substitutions across
nine campaigns with measured hit/non-hit labels**, of which four have deposited complexes.

### 5.1 `dddG_elec`, repacked with PyRosetta — Ahn's actual protocol

| campaign | correct | trivial null | TP | FN | FP |
|---|---|---|---|---|---|
| Pertuzumab (HER2) | 5/8 | 5/8 | 2 | 1 | 2 |
| bH1 (HER2) | 5/8 | 5/8 | 0 | 3 | 0 |
| Ipilimumab (CTLA-4) | 13/18 | 16/18 | 0 | 2 | 3 |
| Tocilizumab (IL-6R) | 18/48 | 42/48 | 1 | 5 | 25 |
| **pooled** | **41/82 = 50%** | **68/82 = 83%** | **3** | **11** | **30** |

**Exactly chance, and well below the null of predicting "no switch" for everything.**

It does work on **multi-histidine** molecules: on Schröter's adalimumab series (3–5 installed
histidines) it is 5 of 5. **The boundary is demonstrated inside a controlled pair:** Adafre's
AF-M2631 shares its heavy chain exactly with Schröter's PSV#3 and differs only by the absence of
two light-chain histidines. PSV#3 (3 His) is called correctly at **+3.390**; AF-M2631 (1 His,
measured 30× switch) is called **wrong at −0.209**. Same heavy chain, same VH histidine.

### 5.2 `his_cation_gate` — Ahn's second, geometric criterion

Fired **2 times in 82 variants, both on non-hits. 0 of 14 confirmed switches.**
A deliberately slack control (8 Å, no lone-pair angle) finds a cation in **66 of 82**, so the
rejection is **geometric, not absence of cations**.

### 5.3 Is there any structural property that does separate them?

Six features computed on the parent structure at the position about to be mutated, permutation
test with 20,000 shuffles: **nothing below p = 0.088** against a Bonferroni bar of 0.0083.

🟠 The only feature with any signal is **cation proximity** — switches sit a median **2.2 Å
closer** to an antigen cation (AUC 0.36, p = 0.088). That is the mechanism above, and our
criterion demanded an H-bond at 3.5 Å along a lone pair, which fires on none of them. **The
intuition may be right while the criterion is far too strict.** Stated as a hypothesis for future
work, not used here.

### 5.4 What this establishes

Across ~100 real molecules — 12 natural TNF complexes and 82 labelled engineered variants — the
mechanism has **never been observed in a confirmed switch**, while cations sit within 8 Å in the
large majority of cases. **Measured pH switches do not work this way.**

⚠️ Stated narrowly: all 82 labelled variants are antibody CDR scans against HER2, CTLA-4 and
IL-6R. TNF-α's epitope is unusually cation-rich, which is why we chose it. This does **not** show
the mechanism cannot work on TNF-α; it shows **there is no measured precedent for it working
anywhere.**

---

## 6. What we ranked on instead

With both purpose-built instruments retired on measured data, ranking uses only quantities with a
passing validation gate on this target.

**Validation gate, pre-registered before any score existed** (`analysis/02-tnf/validation_gate_tnf.py`):
positive = TNFR2 ECD (non-antibody, natural receptor, the receptor whose epitope we design
against); negatives = 4 composition-matched shuffles.

| case | A | B | C |
|---|---|---|---|
| TNFR2 positive | 0.6245 | 0.6161 | 0.4967 |
| all four negatives | 0.0000 | 0.0000 | 0.0000 |

**PASS under all four aggregation rules**, and the separation holds at the **seed** level: the
positive's weakest single seed is 0.4805, the negatives' strongest anywhere is 0.0000.

⚠️ Mean pLDDT is **0.7 for the shuffles against 0.9 for the positives**, so this shows the
instrument separates a binder from an **unfolded** chain. A folded non-binding negative is
_pending_ (§7).

🔴 **The frozen scorer could not score this target.** `ipSAE_min` asserts exactly one inter-chain
pair; a trimeric target gives six. The metric is unchanged — min over two alignment directions of
one interface — but combining three binder:protomer interfaces into one number requires a stated
rule, and all four candidates are reported rather than one chosen silently.

---

## 7. Designs — PENDING

_Waiting on: the target decision (point-mutated crystal vs predicted canonical structure) and
generation. Three BoltzGen arms are built and repointed, histidines pinned at generation rather
than filtered for afterwards; no compute has been spent._

## 8. Novelty and prior art — PENDING the organisers' checker

Sequence-keyed ledger: **22 structures → 119 binder chains → 47 distinct sequences**
(`analysis/02-tnf/prior_art_ledger.py`). Carried explicitly as a **gap**: Ahn et al.'s 72 TNF-α
designs (bioRxiv 2025.09.29.678932), prior art on our exact target and mechanism, sequences not
released.

Added during this work and not in the standard survey: **Adafre's AF-M2637 / AF-M2631**
(PMID 35896334), pH-sensitive monovalent adalimumab variants, **patent pending**.

## 9. Limitations — PENDING, accumulating

Carried so far: the 10 µM ceiling is our assumption about the assay; three different ratio
quantities (KD, off-rate, EC50 shift) appear in the literature and are not interchangeable;
ΔpKa's error bar is 1.0–1.2 units against a claimed 2.0-unit effect, so it is a tiebreak and not
a ranker; and §5's conclusions are drawn on non-TNF antigens.

## 10. Reproducibility

Every figure in §2–§6 has a named artifact under `analysis/02-tnf/` and a self-tested script.
Gate sweep: 20+ gates, all green. Packer trajectories are seeded; a single unseeded repack
carries a measured 13% false-negative rate on the one series where that could be quantified.
