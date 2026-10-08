# Methods — TNF-α conditional binder design (Challenge 2, Track 3)

**Objective.** A de novo binder to human TNF-α that binds at pH 7.4, shows no detectable binding
at pH 6.0, and also binds mouse TNF-α.

**Status of this document.** Written 2026-10-06, before generation. Sections 1–7 and 11–13 are
complete and rest on work already done. Sections 8–10 and 14 describe methods that are fixed but
whose results do not exist yet; they are marked accordingly. Nothing below is written as a
prediction of our own results.

**What this document leads with.** The strongest results in this work are negative: four ranking
instruments were tested against measured data and all four were retired, one of them against a
threshold that was committed to version control before the number existed.

Two findings are reusable independently of whether our designs succeed. **§7** concerns the
scoring instrument: an approved antibody's interface score collapses 18-fold when the target is
modelled as one protomer instead of three, which makes such a score uninterpretable without its
protomer count. **§6.1** concerns multi-objective design generally: a conditional objective can be
weighted above the objective it depends on, and the resulting failure looks like an insufficient
compute budget while in fact being a weight ratio.

---

## 1. Target, and the five numbering schemes

**Construct.** The sequence-corrected human TNF-α crystal trimer, from 1TNF. All three protomers
are identical contiguous slices of **mature residues 6–157** — 152 of the 157 residues the
challenge specifies. The five omitted residues are the disordered N-terminal `VRSSS` (mature 1–5);
nothing is omitted at the C-terminus. Verified per protomer, not assumed:

| protomer | residues | mature range | omitted |
|---|---|---|---|
| A | 152 | 6–157 | `VRSSS` (N-term) |
| B | 152 | 6–157 | `VRSSS` (N-term) |
| C | 152 | 6–157 | `VRSSS` (N-term) |

**The sequence discrepancy, and why it does not apply to us.** 1TNF as deposited carries **Leu**
where the canonical sequence carries **Asp** at mature 143 (canonical D219). Another entrant
reported publicly on 2026-10-05 that binders designed against unmodified 1TNF failed
re-prediction, and that redesigning against the canonical sequence produced binders that bound by
BLI. Our target was point-mutated to **Asp at mature 143** before any design work, and the full
152-residue sequence matches the canonical challenge sequence exactly — checked by string equality
against the sequence published on the challenge page, 2026-10-06, not inferred from the deposition.

**Confirmed by the organisers after the fact, 2026-10-06.** Asked directly to resolve the
discrepancy, they stated the assayed sequence is the one carrying **Asp143**, referenced
**UniProt P01375**, and named the specific reagent under test — **AcroBiosystems TNA-H4211**. This
work's choice therefore agrees with the confirmed target, and the agreement is recorded as
*confirmation of a prior decision* rather than as its basis: the point mutation was made on
sequence evidence a day before the confirmation existed. The vendor construct is the authority on
residue range and tag, and nothing in this work depends on either.

**Independently corroborated by a second campaign, read 2026-10-07.** A published TNF-α
optimization campaign (ArcRefine; manuscript DOI `10.5281/zenodo.23115832`, example data CC BY 4.0)
designed against canonical mature residues **12–157**, carrying **Asp143**. This work's target is
**6–157**, also Asp143; the two constructs differ only by six disordered N-terminal residues, and
that campaign produced a BLI-confirmed binder at this challenge's own assay vendor. So the
construct here agrees with the only TNF-α target on this challenge that has measured binders behind
it, reached by a route independent of both the sequence evidence and the organisers' answer. One
difference is recorded rather than reconciled: that campaign designed against **predicted**
structures, where this work locked the sequence-corrected crystal trimer. Both carry Asp143, so the
gap is small, and no decision here was reopened on the strength of it.

**Numbering.** Five schemes are live in this project, and conflating them has cost real compute
here. The anchor, in all of them:

| scheme | R108 | K166 | R158 | D219 |
|---|---|---|---|---|
| canonical (full-length) | **108** | 166 | 158 | 219 |
| mature (1–157) | **32** | 90 | 82 | 143 |
| our positional (renumbered 1–152) | **27** | 85 | 77 | 138 |
| mouse positional (2TNF, 1–148) | **24** | — | — | — |

Residue identity was confirmed at every mapped position rather than derived arithmetically.
`analysis/02-tnf/epitope_conserved.json` holds the mapping; `bin/mosaic_selftest.py` asserts it
against residue identity on every run.

**Mouse target.** 2TNF, mouse TNF-α, trimer, **100% sequence identity** to the mouse assay
construct, renumbered positionally 1–148.

---

## 2. Epitope selection

The challenge recommends the receptor-binding site at the interface between adjacent protomers.
That site is where the approved antibodies bind, and a binder placed there must straddle the
subunit interface rather than sitting on one chain.

**The epitope was then restricted to positions conserved between human and mouse**, because
objective 2 requires one sequence to bind both. Pairwise alignment of the observed chains gives
79.6% identity overall; within the 13-position epitope, **9 are identical**:

| our positional | human | mouse | kept |
|---|---|---|---|
| 16, 27, 28, 70, 72, 81, 82, 85, 86 | — | — | ✅ |
| 15 | P | H | dropped |
| 26 | R | Q | dropped — charge lost |
| 68 | H | **Y** | dropped — substitution, see the correction below |
| 80 | V | I | dropped |

🔴 **CORRECTION, 2026-10-07. An earlier version of this section said positional 68 is deleted
in mouse. It is not.** Positional 68 (mature 73) is **H → Y**, an ordinary substitution. The
alignment contains **exactly one gap**, and it falls on human **mature 71 = positional 66 (Ser),
canonical 147** — two residues upstream, and not an epitope position. The error was an
*attribution*: the register shift between the two chains is real (−3 then −4 along the epitope) and
the mouse indices derived from it are correct, but the gap was ascribed to the wrong residue.
Re-derived two independent ways — a full-mature-sequence alignment and an observed-chain alignment —
and **the nine mouse epitope indices recompute unchanged**, so no shipped number depended on it.
Stable at gap penalties −4 through −12. Checked by `bin/check_species_map.py`.

**What the retained epitope therefore rests on.** Nine positions are identical by substitution, not
by luck of gap placement. Epitope Cα RMSD mouse-vs-human over them is **1.10 Å** — similar surfaces
on average, with local differences an average conceals. The four dropped positions are all
substitutions (P→H, R→Q, H→Y, V→I); **none is a deletion**, so the earlier claim that a
backbone difference made sidechain selection insufficient does not hold and is withdrawn. The
honest statement is narrower: conservation was obtained by *restriction* — dropping four divergent
positions — rather than by sidechain choice across all thirteen.

**Mouse epitope indices are not the human ones.** The deletion shifts the register by **−3 at the
anchor and −4 further along**: human `16,27,28,70,72,81,82,85,86` maps to mouse
`13,24,25,66,68,77,78,81,82`. Mouse 68, 81 and 82 collide numerically with human 68, 81 and 82
while denoting different residues.

**External corroboration, and exactly how far it goes.** The campaign of §1 released per-design
lineage, assay outcomes and archived coordinates, which lets this epitope choice be checked against a
**measured** binder instead of against our own predictions alone. Its confirmed binder
`TNFA_OPT_07` — apparent K_D **66.2 nM** by BLI — was optimized from an archived starting scaffold
whose footprint we computed ourselves from the released coordinates:

| contact cutoff | target residues contacted | protomers | overlap with our 9 conserved positions | reaches R108 | reaches K166 |
|---|---|---|---|---|---|
| 4.5 Å | 10 | A, B | 2 of 9 | no | no |
| 6.0 Å | 20 | A, B | **7 of 9** | yes | yes |
| 8.0 Å | 38 | A, B | **9 of 9** | yes | yes |

Nearest approach **R108 at 5.12 Å and K166 at 4.92 Å**, spanning the A/B protomer interface — the
same inter-protomer surface and the same two candidate anchors chosen here. The 4.5 Å row reads low
for a mechanical reason, not a structural one: the archived binder is backbone-only (`C, CA, N, O,
OXT`, 4.0 atoms per residue), so an all-atom 4.5 Å cutoff registers only backbone-to-sidechain
contacts and undercounts.

The link from that scaffold to the confirmed binder is **measured retention rather than assumption**.
Re-predicted with the carried optimization state removed, `TNFA_OPT_07` keeps 60–70% of the
scaffold's contacts under seven predictors that supplied no gradients to the design: Boltz-2 0.703,
Protenix 0.676, OpenFold3 0.676, ESMFold2 0.649, OpenDDE 0.622, AlphaFold2 0.595, RF3 0.595.

**What this does not establish, including a counterexample on the same surface.** Every pose above is
predicted. The archived file is explicitly an *input* scaffold whose binder side chains may be `UNK`
placeholders — the release says in terms not to read it as a prediction of the parent sequence — and
the author states the experiments *"do not … establish that binding occurs at the predicted site."*
And the same surface hosts the scaffold of `TNFA_UNOPT_04`, which **was tested and did not bind**.
One scaffold on this epitope led to a confirmed binder; another was itself a measured failure.

**So the claim, at its true strength: this epitope was chosen independently by a campaign that got a
binder on it, and not that it is a measured binding epitope.** A first reading of the same data
overstated this, counting a non-binder's scaffold as support; the correction is recorded in §11.

---

## 3. Mechanism, and its arithmetic

**Design.** A binder histidine positioned against a fixed cation on the target — **R108**. At
pH 7.4 the histidine is predominantly neutral and the interface is intact. At pH 6.0 it is
predominantly protonated, and electrostatic repulsion against the arginine opposes binding.

**Per-site ceiling.** For the pH 7.4 → 6.0 window, the maximum achievable ratio from a single
titratable site, as a function of its free pKa:

| free pKa | 6.0 | 6.5 | 7.0 | 7.4 | 8.0 |
|---|---|---|---|---|---|
| max ratio | 1.92× | 3.70× | 7.87× | 13.1× | 20.3× |

Derived from the linkage relation
`ratio = ((1+10^(pKa_f−6.0))/(1+10^(pKa_b−6.0))) · ((1+10^(pKa_b−7.4))/(1+10^(pKa_f−7.4)))`,
reproduced independently in `analysis/02-tnf/species_and_histidines.py`. The useful ranking is by
**high free pKa**, and single-site designs cannot exceed ~25× in this window regardless of
geometry.

**🔴 The per-site shift we obtain is 0.40 pKa units, and it is a COMPUTED quantity.** PROPKA 3.5.1
on our own poses, neutral-His versus His⁺, paired by pose. It is **not an experimental
measurement**, and it is reported as a computational estimate supporting a testable hypothesis
rather than as an established mechanism.

**What 0.40 units implies.** Two sites give **2.8–3.9×**. Our earlier arithmetic had assumed 2.0
units per site; the computed value is well below the pessimistic end of that band.

**On the target ratio.** The published objective asks for binding at pH 7.4 and **no detectable
binding at pH 6.0**, and names no ratio. The 10× figure used internally in this project is **our
own target, not a challenge requirement** — a point we corrected after it had propagated through
several documents as though it were specified. What "no detectable" means depends on the assay's
detection limits, which we have not been able to confirm for this challenge and which we therefore
do not assume.

**Known-answer control, 9/9.** PROPKA 3.5.1 against JustHISpKa on apo 1TNF, all three protomers,
chain against chain: mean |Δ| **0.49**, RMSE **0.85**, max 2.40; excluding a single outlier at
chain C H91, mean |Δ| **0.26**.

---

## 4. Four ranking instruments, tested against measured data and retired

This is the section this work is most confident in. Each instrument was tested against molecules
whose behaviour was measured in a laboratory, and each failed.

| instrument | test | result |
|---|---|---|
| `dddG_elec` (Ahn et al. fa_elec ΔΔG) | 82 labelled single-His variants | **41/82 — chance**, against an **83% trivial null** |
| `his_cation_gate` (Ahn's geometric criterion) | 14 measured switches | **0 of 14**, while a slack control finds cations in **66 of 82** |
| structural features generally | 82 variants | nothing separates switches below **p = 0.088** |
| cross-species ranking | TNFR2 + human vs TNFR2 + mouse | **unresolvable** — the two cannot be distinguished |
| **`ipSAE_min` as an affinity ranker** | **measured 4.6–112 pM ladder** | **ρ = 0.400 against a pre-registered bar of 0.8** |

**The `ipSAE_min` test was pre-registered.** The aggregation rules and the success criterion were
committed to version control at **16:24:41 on 2026-10-06**, before any number was computed; the
commit timestamp precedes the scoring run. Three rules were declared in advance and all three
reported:

| rule | ρ | p (exact, 24 permutations) | predicted order |
|---|---|---|---|
| `max_pair` | **0.400** | 0.375 | adalimumab > **PSV3** > PSV1 > PSV2 |
| `sum_pairs` | **0.400** | 0.375 | adalimumab > **PSV3** > PSV1 > PSV2 |
| `top_contact` | **−0.200** | 0.625 | PSV3 > adalimumab > PSV1 > PSV2 |

True order, tightest to weakest: adalimumab (4.6 pM) > PSV#1 (46.3) > PSV#2 (77.3) > PSV#3 (112),
from Schröter et al., *mAbs* 7(1):138–151, DOI 10.4161/19420862.2014.985993. **The instrument
places the weakest binder in the panel second and the middle of the ladder last.** The rule whose
pair selection was made on inter-chain contact counts — and therefore could not see the score it
was about to report — is the one that comes out anti-correlated.

**The second failure is dynamic range.** All four constructs occupy **0.42–0.58** across a **24×**
affinity range, a separation smaller than the seed spread on most designs. The consistent reading
is that `ipSAE_min` separates binder from non-binder and carries no information about affinity
*among* binders — which is the only regime a design campaign operates in, since every candidate it
ranks is already a predicted binder.

**Scope of the claim, stated narrowly.** n = 4, so ρ = 0.400 is also not distinguishable from
chance. The bar was set in advance precisely so that neither reading could be argued after the
fact. This result retires **our affinity-ranking use** of the method under a prespecified
criterion. It is not a general claim that the method is random.

**Consequence for this submission.** No instrument available to us ranks candidates by predicted
affinity. The submitted ordering is therefore constructed on a documented nomination priority
(§10) and is labelled provisional; no retired score is presented as a calibrated prediction of
affinity or of pH selectivity.

---

### 4.5 All of the above, recalibrated against 150 measured designs on this target

Added 2026-10-07. `arms-backlog.md` recorded a per-target positive control as a **blocking
dependency** for validating any instrument. It was public throughout. Anthropic's released campaign
data (`huggingface.co/datasets/Anthropic/claude-protein-binder-design`, CC BY 4.0, ungated) carries
**150 de novo TNF-α designs assayed against the same construct as this challenge, Acro TNA-H4211**:
12 measured binders, 138 measured non-binders, each with co-folding metrics from ten models.
Fetched by `bin/fetch-anthropic-campaign.sh`; calibrated by `analysis/02-tnf/calibrate_external.py`,
which carries a self-test and two mutation tests. Full output, all 160 model × stoichiometry ×
metric cells, in `analysis/02-tnf/external_calibration_results.txt`.

Design-level AUC, binder vs non-binder, at **1 binder : 3 protomers** — the construct we submit.
Seeds are collapsed per design by median before ranking, because they are replicates of one design
and raw-seed ranking would inflate n tenfold.

| model | `ipsae_min` | `pae_interface_min` | `plddt_binder` |
|---|---|---|---|
| rf3 | 0.951 | **0.960** | 0.701 |
| of3 | 0.908 | 0.943 | 0.766 |
| af3of3 | 0.877 | 0.944 | 0.762 |
| **ef2full — our arm** | 0.822 | **0.901** | 0.579 |
| **boltz2 — our design oracle** | **0.781** | 0.820 | 0.568 |
| **ef2fast — our other arm** | 0.692 | 0.745 | 0.471 |

The top cells reach **p = 0.00005** by permutation on the labels (n = 20000) against a Bonferroni
bar of 0.00031 for the 160 cells scanned.

**Three conclusions, two of which change what we ship.**

**(a) `ipSAE` is a ranker for binding, and `plddt_binder` is not.** The three most anti-correlated
cells of all 160 are `plddt_binder` (AUC 0.413–0.471). This **confirms** §4's reading rather than
reversing it: ipSAE separates binder from non-binder and carries no affinity information among
binders. The retirement in §4 was measured on an affinity ladder — a different task, still retired
for that task.

**(b) We rank on an ESMFold2-Full interface-PAE term**, with `ipSAE_min` reported alongside.
Implemented in `bin/pae_interface.py` as a **sibling** of `ipsae_min.py` rather than an edit to it,
because that module is pinned by 12 fixtures and a fail-closed suite. Validated on those fixtures:
barnase/barstar reads **0.332, best of 12**; the shuffled negative **19.137, worst of 12**; and the
renumbered and chain-swapped variants return **0.332 identically**, so the metric is invariant to
the chain relabelling this project has repeatedly been bitten by.

⚠️ **The 0.901 that selected `pae_interface_min` was measured on a configuration we never ran, and
we no longer rest the choice on it.** The calibration's `ef2full` arm uses **target-chain MSAs**
(`msa_max_seq=2048`, protocol line 111). Our run is **single-sequence on both chains**: `grep -ci
msa` returns **0** in both the wrapper and its tracked patch. The AUC separating 0.901 from
`ipsae_min`'s 0.822 is therefore not transferable to our numbers, and
`analysis/02-tnf/CONFIG-COMPARABILITY-2d.md` records that gap rather than papering over it.

**What replaced it.** All 35 designs were scored on our own configuration, and each was compared to
a **size-matched shuffled null** at its own cross-chain pair count — the earlier control band paired
a 164-residue receptor against the trimer (302,382 pairs) while the designs sit at 219,486/227,022,
so it was never a yardstick (`analysis/02-tnf/SIZEMATCHED-NULL-2f.md`):

| | n | design `pae_if_mean` | median | own-length null | beat it |
|---|--:|---|--:|--:|--:|
| L76 | 21 | 11.218–16.586 | 14.047 | 14.248 | 15/21 |
| L84 | 14 | 13.945–15.434 | 14.648 | 14.996 | 10/14 |

**25 of 35 beat their own-length null** (one-sided binomial *p* = 0.0083), and **the same 25 beat it
on the min**, so at our configuration the two metrics order the set identically and the choice
between them no longer carries weight. The margin is what matters: TNFR2 scores 7.713 where its own
shuffles score ~18.99, a margin of **11.28**. Our best design clears its null by **3.03 (27%)**, our
best L84 by 1.05 (9%), and the median by **0.20–0.35 (2–3%)**. Statistically above chance;
biologically almost nothing.

**(c) ESMFold2-Full and ESMFold2-Fast are not interchangeable.** `arms-backlog.md` called their
agreement to four decimal places *"a symptom rather than a virtue"*. On external labels they split
**0.822 versus 0.692**. The complaint was right and was measuring the wrong thing. **Fast is retired
for this target.**

### 4.6 The seed-noise floor: what a single-seed margin can and cannot mean

**This is a property of the SCORER, not of any generator.** It measures how far
`pae_interface_mean` moves when the same complex is folded again with a different random seed, so
it applies to every single-seed number in this submission regardless of which method produced the
design.

Their released data carries **exactly 5 seeds per design × model × stoichiometry** — 3,000 groups,
15,000 rows, with a `seed` column. Ours carries **one**. `bin/score-esmfold2.sh` defaults to the
protocol's FINAL tier of five; every problem-2 run in this work overrode it to a single seed for
cost. That was a defensible economy and it has a consequence that was not priced at the time.

On our own arm, at our own stoichiometry, over their 150 designs × 5 seeds
(`analysis/02-tnf/per_seed_spread_2026-10-08.tsv`):

| | `ef2full`, 1:3 |
|---|--:|
| within-design seed SD, `pae_interface_mean` | **0.607** |
| within-design seed SD, `pae_interface_min` | 0.375 |
| within-design seed SD, `ipsae_min` | 0.050 |

A margin in this work is `null_mean − design_mean`, a **difference of two single draws**, so its
noise is `0.607 × √2 = `**`0.858`**.

| margin | value | in SDs of that difference |
|---|--:|--:|
| our **median** design | 0.27 | **0.3** |
| our **best** design | 3.03 | 3.5 |
| TNFR2, a real receptor | 11.28 | 13.1 |

**And a best-of-N is selected, so noise alone yields a positive best margin.** Simulated over
20,000 trials at this SD:

```
best-of-35 under pure noise:  mean 1.81,  95th percentile 2.56
best-of-48 under pure noise:  mean 1.91,  95th percentile 2.63
```

**Three consequences, stated against our own result.**

The **median** margin of 0.27 is **0.3 SDs**. It is indistinguishable from seed noise. §4.5's
"25 of 35 beat their own null, binomial *p* = 0.0083" stands as a **sign** test, and the *size* of
the typical margin does not survive this floor.

Our **best** design at 3.03 is 3.5 SDs — but it is a best-of-35, and pure noise reaches 2.56 at the
95th percentile for that N. So 3.03 is **only modestly beyond what selecting the maximum of 35
noisy draws produces anyway**, and we do not present it as a clear signal.

**The floor is probably optimistic for us.** Their `ef2full` arm uses target-chain MSAs and ours is
single-sequence (§4.5b). More input information generally stabilises a prediction across seeds, so
0.607 is better read as a **lower bound** on our own noise than an estimate of it. Measuring our
actual figure needs a multi-seed re-run that was not made.

**What would fix it, and what it would cost.** Five seeds per complex reduces the noise on a
median by roughly `√5`, taking the margin's SE from ~0.86 to ~0.48. At ~94 s and ~$0.05 per
548-residue fold, five seeds on 48 designs plus 48 nulls is ~$24 and ~12 GPU-hours. That is
affordable and was not spent; the single-seed economy is recorded here rather than defended.

---

## 5. What the surviving instrument can do

`ipSAE_min` reads a labelled negative correctly. Adalimumab is TNF-α-specific and does not
neutralise lymphotoxin-α:

| complex | ipSAE_min | n |
|---|---|---|
| adalimumab + human TNF-α trimer | **0.5757** | 5 seeds |
| adalimumab + LT-α trimer | **0.0000** | 5 seeds |

No seed disagreement. The LT-α trimer's own internal packing reads **0.79** in the same runs, and
the Fab's VH/VL interface **0.82–0.83**, so the zero is a property of the interface rather than a
failed prediction — a distinction the first attempt at this control could not make, because its
construct had not assembled at all (§11).

**This supports using the method as a gate, not as a ranker.** It is reported as a
structural-confidence annotation. **No candidate is rejected on it automatically**, because five
seeds are repeated predictions of one pair rather than five independent biological controls, and
one labelled pair cannot establish sensitivity or specificity for de novo binders.

---

## 6. The generation method

**Stack.** Mosaic (gradient-based multi-objective design) over Boltz-2, with soluble ProteinMPNN
inverse-folding recovery, optimised by `simplex_APGM`. Cysteine is removed from the alphabet
entirely rather than penalised.

**Three loss legs, and one of them is unusual.** The loss is a sum over:

1. **human complex** — all three protomers, 456 target residues, 532 tokens
2. **species complex** — mouse TNF-α, binding and pH terms
3. **binder alone** — so the binder is required to fold without the target

Leg 2 means **objective 2 is optimised rather than checked afterwards.** Published work finds that
ipSAE-type scores cannot resolve species (§4), which is a statement about *selecting* among
finished molecules and not about whether cross-reactivity can be *designed for*. The species leg
carries the binding and pH terms but **not** the folding terms, since folding is a property of the
binder alone and is already constrained by legs 1 and 3; including it a third time would
outweigh the objectives the leg exists to add.

**The species leg is run as a two-protomer target**, not three. This is a measured hardware limit
rather than a modelling choice: three mouse protomers is 976 tokens and exhausts a 48 GB L40S with
a single 34.24 GiB allocation inside the gradient computation. Two protomers run at a measured
25.9–30.6 s/step. A separate control (§7) establishes that two protomers recover 86% of the
three-protomer interface signal.

**The pH term.** A binder histidine is pulled toward the target cation, measured CA-to-cation-N
because design-time features give the binder no sidechains. The reduction over per-position scores
is **top-2, not a sum**. This matters, and the reason is recorded in §11.

**Schedule.** 150 soft + 50 sharp steps. The pH term has an effective range of ~20 Å and is
numerically flat beyond ~30 Å, where trajectories begin — so the contact terms must bring the
binder into range before the pH objective contributes any gradient at all. A short schedule spends
its budget in the regime where the mechanism is invisible.

### 6.1 The weight balance, and a tuning pathology worth reporting

Multi-objective design has a failure mode that is rarely written down: an objective can be
weighted strongly enough to compete with the objective it depends on. We hit it, and the diagnosis
is reportable independently of whether our designs succeed.

**The shipped balance was:**

| term | weight | what it asks for |
|---|---|---|
| `BinderTargetIPTM` | 1.0 | interface confidence |
| `BinderTargetContact` | 1.0 | contact with the epitope |
| `BinderTargetPAE` | 0.05 | interface precision |
| **pH term** | **2.0** | **a binder histidine against the target cation** |

**The pH objective carried twice the weight of the two binding objectives combined with
themselves** — on designs that were not binding. Measured at 50 optimisation steps, four
trajectories:

| | iptm_repred | nearest His → anchor |
|---|---|---|
| L76 seed 0 | 0.153 | 18.0 Å |
| L76 seed 1 | 0.120 | 24.8 Å |
| L84 seed 0 | 0.111 | 47.9 Å |
| L84 seed 1 | 0.156 | 39.2 Å |

**The distances are not measurements of histidine placement.** They are distances measured on
poses whose interface confidence is **0.11–0.16 against a 0.45 threshold** — i.e. poses where the
binder is not bound to anything. Our own gate (§8) classifies all four as *not applicable* rather
than as placement failures, which is the distinction the gate exists to enforce.

**Two independent observations point at the weighting rather than at the step budget.**

1. **Step count did not move the interface.** Five-step runs on the same configuration gave
   iptm_repred 0.178; fifty-step runs gave a median of **0.135**. A tenfold increase in
   optimisation produced no improvement in the quantity that has to improve first.
2. **The pH term cannot act at these distances.** Its sigmoid, at `d0 = 6.5 Å` and
   `width = 1.5 Å`, reads 1.2 × 10⁻⁴ at 20 Å and 1.4 × 10⁻¹¹ at 44 Å. It is numerically flat
   where the trajectories sit. So the weight it carries is not buying histidine placement — it is
   competing with the terms that would bring the binder close enough for placement to be
   achievable at all.

**The mechanism of the pathology.** A gradient-based multi-objective design has no notion that one
of its objectives is a precondition for another. The pH term is only meaningful inside ~20 Å; the
contact terms are what get the binder there. Weighting the conditional objective above the
precondition asks the optimiser to satisfy a constraint in a regime where the constraint is
unmeasurable, at the expense of reaching that regime. **The fix is a weight ratio, not more
compute** — which is the opposite of the conclusion a step-count sweep would have reached.

**Why it went unnoticed.** Six of the eleven loss weights — including both binding terms — were
present on the remote design function but absent from the command-line interface. Retuning the
balance required editing the source, so the shipped ratio was never treated as a parameter. **A
weight that cannot be passed is a weight nobody tunes.** All eleven are now exposed with identical
defaults, and the test suite asserts that every weight on the design function is reachable from
the interface with a matching default — verified by deleting one and confirming the test names it.

### 6.2 The weight ratio is a lever, measured — and it did not reach the threshold

A controlled single-variable probe tested the hypothesis: identical configuration and step count,
**pH weight reduced from 2.0 to 0.5**, four trajectories against the four-trajectory baseline in
§6.1. Since only the ratio of weights affects a weighted sum, reducing the pH weight is equivalent
to raising the binding weights. The success criterion was **fixed before the probe was launched:
a median iptm_repred of 0.20 or above.**

| | iptm_repred | binder pLDDT | histidine content | nearest His → anchor |
|---|---|---|---|---|
| L76 seed 0 | 0.175 | 0.51 | 2.6% | 24.7 Å |
| L76 seed 1 | 0.115 | 0.69 | 3.9% | 27.9 Å |
| L84 seed 0 | 0.183 | 0.57 | 7.1% | 21.0 Å |
| **L84 seed 1** | **0.200** | 0.37 | 2.4% | **2.62 Å** |

| | baseline (pH weight 2.0) | probe (pH weight 0.5) |
|---|---|---|
| median iptm_repred | 0.136 | **0.179** |
| max iptm_repred | 0.156 | **0.200** |

**Against the prespecified criterion, this is a miss: the median is 0.179 against a bar of 0.20.**
Reported as registered — median, not maximum.

**The ratio is nonetheless a measured lever.** A single-variable change moved the median by
**+31%**, and the quantity it moved is the one that has to move first.

**A second observation, not prespecified and therefore reported separately.** One trajectory placed
its histidine **2.62 Å from the anchor nitrogen** — inside the 4.0 Å criterion. Across the
preceding four-trajectory baseline the same measurement read 18.0, 24.8, 39.2 and 47.9 Å. The
geometry verdict for this design is nonetheless recorded as **not applicable**, because its
interface confidence of 0.200 is far below the 0.45 threshold at which this work treats a geometry
measurement as meaningful (§8); the distance is real, the placement is not verified.

**This result also falsifies the obvious objection to the intervention.** Reducing the weight on
the pH objective might be expected to reduce histidine placement. The opposite was observed: the
closest placement this work has produced came from the run in which the pH objective was weighted
**four times lower**. That is consistent with the mechanism proposed in §6.1 — the term was not
limited by its own weight but by the binder never entering the ~20 Å range in which the term has
any gradient. Allowing the binder to arrive does more for the conditional objective than weighting
the conditional objective more heavily.

**Scope, stated narrowly.** n = 4 per condition. The median difference (0.136 → 0.179) rests on
four trajectories per condition and is not a significance claim. The 2.62 Å placement is a single observation. Neither
result establishes that this configuration can produce a bound complex: **the best interface
confidence obtained anywhere in this work is 0.200 against a 0.45 threshold**, and no design has
yet reached a value at which this work would report a geometry verdict at all.

**Next test, launched before these results were written up.** Binding weights raised directly —
`w_iptm` and `w_contact` from 1.0 to 3.0, pH weight held at 0.5 — which was not expressible from
the command line until the interface gap described in §6.1 was closed. Same step count, same four
trajectories, same criterion. Its result, and the step-budget test that followed it, are §6.3.

*Results from the production method do not yet exist. Nothing in this section is a claim about
submitted designs.*

### 6.3 Two more conditions, and the measurement that ended the search

Two further single-variable tests followed §6.2, both prespecified, both against the same bar of a
**median iptm_repred of 0.20**.

**Condition 3 — binding weights tripled.** `w_iptm` and `w_contact` raised from 1.0 to 3.0, pH
weight held at 0.5, step count unchanged. If the interface was under-weighted, this is the direct
intervention.

**Condition 4 — step budget quadrupled.** 50 steps to 200 (150 soft + 50 sharp), weights identical
to §6.2's probe, **same two lengths and same two seeds**, so every one of the four trajectories is
a matched pair with a §6.2 trajectory and the step count is the only difference.

| condition | pH weight | binding weights | steps | median | max | closest His |
|---|---|---|---|---|---|---|
| 1 — baseline | 2.0 | 1.0 | 50 | 0.136 | 0.156 | 18.0 Å |
| 2 — pH weight down | 0.5 | 1.0 | 50 | **0.179** | **0.200** | **2.62 Å** |
| 3 — binding weights up | 0.5 | 3.0 | 50 | 0.156 | 0.197 | 9.69 Å |
| 4 — step budget up | 0.5 | 1.0 | **200** | 0.163 | 0.189 | 10.4 Å |
| 5 — free footprint | 0.5 | 1.0 | 100 | 0.134 | 0.137 | 20.1 Å |

**All five miss the bar. Condition 2 remains the best, and it was the second thing tried.**

#### The matched pairs

Condition 4 is the only paired comparison in this work, so it is reported pair by pair:

| length / seed | 50 steps | 200 steps | difference |
|---|---|---|---|
| 76 / 0 | 0.175 | 0.161 | −0.015 |
| 76 / 1 | 0.115 | 0.189 | +0.074 |
| 84 / 0 | 0.183 | 0.165 | −0.018 |
| 84 / 1 | 0.200 | 0.150 | −0.050 |
| **mean** | **0.168** | **0.166** | **−0.002** |

**Four times the compute changed the mean by −0.002.** Three of the four pairs got worse. The one
that improved was the weakest trajectory of the 50-step set, which is what regression to the mean
looks like and is not evidence of a step effect. Against this, the spread produced by seed and
length alone across the eight runs is **0.115 to 0.200** — a range of 0.085, roughly **forty times**
the size of the step effect.

#### The trajectories, which are the actual result

The step-budget test was designed to read a curve shape, not only an endpoint. The prespecified
reading was: if the design-time interface score plateaus by step 60–80, more steps are useless; if
it is still climbing at step 200, more steps are the answer. **Neither was observed.** Averaged in
25-step blocks — necessary because each step re-predicts with a stochastic structure-prediction
pass, so single-step values bounce by more than the whole optimisation moves:

| run | 0–24 | 25–49 | 50–74 | 75–99 | 100–124 | 125–149 |
|---|---|---|---|---|---|---|
| 76 / 0 | 0.132 | 0.126 | 0.133 | 0.130 | 0.132 | 0.126 |
| 76 / 1 | 0.112 | 0.133 | 0.162 | 0.126 | 0.127 | 0.132 |
| 84 / 0 | 0.137 | 0.120 | 0.134 | 0.139 | 0.144 | 0.131 |
| 84 / 1 | 0.133 | 0.129 | 0.122 | 0.123 | 0.124 | 0.142 |
| **mean** | **0.129** | **0.127** | **0.138** | **0.130** | **0.132** | **0.133** |

**The interface term does not climb.** It does not plateau after rising; there is no rise. Across
150 gradient steps the four-run mean moves from 0.129 to 0.133, which is smaller than the
step-to-step noise and smaller than the movement between any two adjacent blocks.

**The same measurement was then run on the eight 50-step trajectories from conditions 1–3**, where
it had not been looked at before. All eight are flat on the same reading, spanning 0.115 to 0.152
with no trend in any of them. So the finding is not specific to the long runs:

> **In twelve gradient-descent trajectories across three weight conditions and two step budgets,
> the interface objective never improved.**

**This is reported as the central negative result of the generation method**, because it changes
what kind of problem this is. A median that misses a bar invites more tuning. A flat trajectory
says the tuning surface is the wrong place to look: the optimiser is descending — total loss moves,
the composition and histidine-content constraints are satisfied and held — but it is not descending
on the interface. Searching weight space and step space more finely cannot fix a term that is not
responding to either.

#### What it does not establish

**It does not identify the cause**, and three candidates remain open:

1. **The pinned epitope.** The binder's footprint is restricted to nine specified positions. If no
   gradient path exists from the initial pose to a bound pose inside that restriction, the term
   would be flat for a reason that has nothing to do with weights or steps. This is testable by
   removing the restriction, and the free-footprint family (§10) does exactly that — a design
   choice made for an independent reason that now also serves as this test.
2. **Gradient quality through the structure predictor.** The per-step interface value swings between
   0.11 and 0.25 under a fixed sequence-space neighbourhood. If the gradient is dominated by that
   sampling noise, no step count recovers signal.
3. **Competing terms.** Eleven weighted terms share one scalar. Condition 3 is evidence *against*
   the simplest version of this — tripling the interface weights made the median worse, not better,
   and pushed histidine content through its 8% cap — but it does not rule out a subtler interaction.

**It also does not establish that the designs are bad**, only that this loss was not optimising
the quantity the loss intended. The distinction matters for §8: every reported number comes from an
independent re-prediction, not from inside a trajectory.

**Cost of the four conditions: approximately $26 of compute, against a production wave budgeted at
$49.** The negative result was bought at roughly half the price of the thing it was protecting.

*Results from the production method do not yet exist. Nothing in this section is a claim about
submitted designs.*

#### Reproducing the trajectory reading

The twelve run logs are committed at `analysis/02-tnf/traj/`, the endpoint table at
`analysis/02-tnf/deep_vs_probe.tsv`, and the block-averaging tool at
`analysis/02-tnf/loss_traj.py`, which carries a self-test:

```
analysis/02-tnf/loss_traj.py --selftest
analysis/02-tnf/loss_traj.py analysis/02-tnf/traj/deep_*.log
analysis/02-tnf/loss_traj.py --block 13 analysis/02-tnf/traj/probe_*.log
```

The self-test covers the two ways this reading can be made to lie: the step counter restarts at
zero when the soft phase hands off to the sharp phase, so a naive parse concatenates two phases
into one apparent trajectory; and a term absent from a log line must come back as absent rather
than as zero, since a zero would manufacture a downward trend.

### 6.4 Removing the epitope restriction: the one testable cause, eliminated

§6.3 left three candidate causes for the flat interface term and only one was testable in the
remaining time: that the nine pinned epitope positions admit no gradient path to a bound pose. The
test removes the restriction — `--epitope none` on both legs, so contact is rewarded against any of
the 456 human and 296 species target residues and the optimiser selects its own footprint. **The
anchor stays pinned in both arms**, because the pH mechanism needs its one cation contact; this is
a free footprint with a pinned anchor, not an unconstrained affinity problem.

Everything else is held: the same weights as condition 2, the same two lengths, the same two seeds.
100 steps rather than 50, chosen because a flat result at 50 could not be distinguished from "the
binder never had time to find a site" — and because the soft phase does not anneal against its own
length, so the longer run contains the matched 50-step comparison as a prefix.

**The restriction was not the blocker.**

| length / seed | pinned, 50 steps | free, 100 steps | difference |
|---|---|---|---|
| 76 / 0 | 0.175 | 0.135 | −0.040 |
| 76 / 1 | 0.115 | 0.137 | +0.022 |
| 84 / 0 | 0.183 | 0.133 | −0.050 |
| 84 / 1 | 0.200 | 0.134 | −0.066 |
| **mean** | **0.168** | **0.135** | **−0.034 ± 0.039** |

The trajectory is flat on the same reading as the other twelve: a four-run rise of
**+0.004 ± 0.013** across 75 soft steps, with one run marginally climbing, one falling and two
flat. **Twice the step budget of the pinned arm, a footprint 50× larger to choose from, and the
interface objective still does not improve.**

The endpoint difference is reported with its interval because it is **negative**: freeing the
footprint made the final score worse by 0.034, and the 2 SE interval is [−0.072, +0.005]. The
honest statement is **worse or no different, and certainly not better** — this work does not claim
a significant degradation from n = 4.

#### Two observations that were not predicted, and one that matters mechanistically

**The histidine moved away from the anchor.** Median nearest-histidine distance by condition:
22.8 Å pinned at 50 steps, 23.6 Å pinned at 200, and **34.8 Å free** — the worst of the three, with
one run at 54.5 Å. This is coherent rather than surprising once stated: the anchor is a single
pinned point on a 456-residue surface, and the epitope restriction was the only term holding the
binder in the anchor's neighbourhood. **Removing it let the binder drift away from the one contact
the mechanism requires.** The pinned epitope was not an obstacle to the pH objective; it was its
scaffold.

**The free-footprint scores cluster far more tightly than the pinned ones.** Four runs spanning two
lengths and two seeds returned 0.1327, 0.1336, 0.1351 and 0.1372 — a range of **0.0045**, against
0.085 for the pinned condition. Seed and length, which dominated every previous comparison, stopped
mattering. The reading offered, and it is an interpretation rather than a measurement: with no
epitope term the optimiser converges on the same generic surface contact regardless of where it
starts, and that solution is reproducible and mediocre. A tight distribution at a low value is not
better than a wide one; it is a sign of a single attractor that is not the one wanted.

#### What this eliminates, and what survives

**Cause 1 is eliminated.** The footprint restriction does not explain the flat interface term,
because removing it does not change the term. Causes 2 and 3 of §6.3 — sampling noise in the
structure predictor swamping the gradient, and interaction among the eleven weighted terms —
survive. **Correction, 2026-10-07: cause 2 is testable, and the claim that it was not was wrong.**
Mosaic's own API at the pinned commit exposes `build_multisample_loss(..., num_samples=4)`, which
re-runs the structure and confidence modules several times from a single trunk output and averages
— four samples halve the gradient-noise SD for well under 4× the cost. Momentum, the other half of
the same question, is two hardcoded literals in our wrapper (0.9 soft, 0.5 sharp). **Every one of the
sixteen trajectories drew one sample per gradient evaluation and integrated it under momentum 0.9,
on a per-step interface SD of 0.0345 against a signal of ~0.13**, which is a mechanism for a flat
trajectory on its own. That configuration was never varied, and calling it untestable was a failure
to read the dependency's API rather than a fact about the budget. Cause 3 remains untestable here.

**A prediction that failed, recorded because it was made in writing first.** The five-step smoke
run that validated this code path returned an interface confidence of **0.2932, the highest figure
this project has produced**, against a previous best of 0.200. It was recorded at the time as not
evidence — n = 1, five steps optimising nothing, binder pLDDT 0.34, and a histidine fraction of
13.2% against an 8% cap. **The four real runs maxed at 0.1372.** The flagged number did not
reproduce, and the reason to record this is that it was the most interesting-looking result of the
day and acting on it would have been wrong.

**Total: 16 trajectories, 5 conditions, 2 step budgets, 2 footprint families. The interface
objective improved in none of them.** The generation method does not produce bound complexes
against this target, and the search space in which it might have — weights, steps, footprint — is
exhausted at this budget. Cumulative compute: approximately $36.

*Results from the production method do not exist, and on this evidence will not. Nothing in this
section is a claim about submitted designs.*

---

### 6.5 What four other campaigns measured on this target

The search above ended flat. Before reading that as a defect peculiar to this method, it is worth
recording what other groups measured on the same target. The comparison below was assembled on
2026-10-07 from the published TNF-α campaign of §1 and its cited references, and it is external to
this work in every row.

| method | TNF-α outcome |
|---|---|
| AlphaProteo | **0 binders of 54 tested** |
| BoltzGen | no binders; missed TNF-α in **both** the nanobody and the protein arm |
| PXDesign | no binders — *"our pipeline encountered difficulty with one challenging target, TNF-α"* |
| Complexa | no binders |
| the challenge organisers' own TNF-α campaign | **12 of 150 = 8.0%**, best apparent K_D 0.70 nM |
| ArcRefine | 6 of 10 optimized vs 1 of 10 unoptimized, n = 10, selection-confounded (below) |

⚠️ **8.0% is a ONE-OBJECTIVE ceiling, not our calibration.** We have quoted it as the benchmark in
earlier drafts of this section, of limitation 19 and of the public methodology box. Their campaign
required **binding only**. Their protocol states the goals as *"1) high-affinity binders zero-shot
and 2) high overall hit rate"*, contains **no pH requirement at any point**, and makes cross-species
*"a secondary objective pursued only without compromising affinity or hit rate."* We are asked for
**three things at once**, and the added one has no computational precedent at pH 6.0 (§6.8). So the
honest framing of our negative result is not "we fell short of 8%" but "this is the
three-objective version of a task whose one-objective version runs at 8% on the organisers' own
instrumentation." Their generator breakdown is the sharper comparison, and it is in §6.9.

**AlphaProteo's published reason for failing is a description of the epitope this work chose.** It
names *"a flat, highly polar binding site at an interface between 2 subunits in a homotrimer."* That
is the protomer-spanning surface of §2, and it is the most direct explanation available for an
interface term that did not improve in sixteen trajectories — an explanation that is neither a
tuning fault nor reachable by any weight, step budget or footprint this work could have varied. It is
offered as the leading external hypothesis, not as a measurement of these runs, and it was found
*after* the search closed rather than used to justify closing it.

**The honest calibration for this target is 8%, not 60%.** A competent pipeline here returns roughly
one hit in twelve. The 6-of-10 figure is a hit-rate delta between two pools selected independently
from a common top-60 cohort — the source states they *"do not constitute matched experimental
parent–child pairs"* — so it is not a per-design rescue rate and should anchor nothing.

**Two consequences, both recorded before the submission closed.** First, the fully-specified
BoltzGen arm in this repository (`targets/tnf/boltzgen_tnf_2his.yaml`, histidines pinned at
generation, repointed at the corrected target) was deliberately **not launched**: BoltzGen is on the
zero-binder list above, and none of that campaign's 60 starting designs came from BoltzGen either, so
launching it on the theory that it produces binders by itself is unsupported. Second, the one axis
this work never varied is the **optimizer**. All sixteen trajectories held the gradient estimator and
the initialization constant while weights, steps and footprint were varied, so "the search is
exhausted" is true of weight space, step space and footprint space, and not of optimizer space. That
is stated here because it bounds the claim, and the budget to test it was not spent.

**The method that reports rescuing failed designs is not the one examined here.** HalluDesign
(bioRxiv 2025.11.08.686881) carries *coordinates* rather than trunk representations, noising the
previous iteration's structure through a truncated AF3 diffusion trajectory, and it reports wet-lab
rescue of designs that had already failed: **12 of 16** on PD-L1 and **8 of 16** on IL7RA, each
candidate derived from a previously failed design, against RFdiffusion baselines of 12/95 and 32/96.
Preprint, n = 16, with baselines taken from another group's paper. It is AF3-based and was out of
reach in the days remaining. It is cited as the honest pointer to what this work would try with more
runway — nothing here used it.

---

### 6.6 The estimator probe: built, pre-registered, and reported whatever it returns

Acting on the correction above rather than only recording it. The configuration every one of the
sixteen trajectories shared is now a parameter, and an arm exists that varies it:

| held identical to condition 2 — the best of the five | varied |
|---|---|
| pinned 9-position epitope, anchor 27 on chain B | **`--grad-samples 4`** (was 1) |
| `w_acid 0.5` and the condition-2 weight set | **`--momentum-soft 0.0`** (was 0.9) |
| lengths 76 and 84, seeds 0 and 1 | — |
| `steps_soft 38 / steps_sharp 12`, `sampling_steps 10`, `recycling_steps 1` | — |
| species leg on, mouse dimer, `w_species 1.0` | — |

Four samples go through `build_multisample_loss` on the **human and species legs only**. The monomer
leg stays single-sample: it carries pLDDT, within-binder contact and globularity and **no interface
term**, so averaging it spends GPU memory on a leg that cannot answer the question. **Correction,
stated because an earlier draft of this section had it wrong:** the human trimer is *not* the memory
ceiling — 456 residues / 532 tokens runs at 21.7 s/step. What exceeds the card is the dual-species
graph with a mouse trimer (~976 tokens); the mouse dimer at 828 tokens runs at 25.9 s/step, which is
why the species leg uses two protomers (§14).

**The bar is pre-registered, before any number exists**, and is the same statistic and tool as the
free-footprint probe so the two are comparable: **four-run mean rise in `iptm_repred` across the
soft phase**, via `analysis/02-tnf/loss_traj.py --block 13`. The free-footprint probe returned
**+0.004 ± 0.013**. A rise worth acting on is **≥ +0.12**; detectable rise over four runs is ~0.017,
so sensitivity is not the constraint. **A mean rise under +0.02 eliminates cause 2**, leaves cause 3
alone, and ends the search — it is not grounds for a seventh condition. Endpoint medians are
secondary, for the §6.3 reason that an endpoint cannot distinguish a badly-tuned method from one
that is not optimising.

**What it does not test, stated in advance.** Cause 3 is untouched, and condition 3 is direct
evidence for it — tripling the interface weights made results *worse*, on 20 terms across three legs
where the external method carries 8 on one leg at an ipTM weight of 0.025. A null here does not
exonerate the loss, and a rise would not prove noise was the only problem. AlphaProteo's diagnosis
of this epitope (§6.5) is optimizer-independent and survives either result.

**Status: RAN 2026-10-07, and returned a null.** Four trajectories, matched to condition 2 with the
estimator as the only variable. Four-run mean rise in `iptm_repred` across the soft phase:
**+0.0078 ± 0.0124**, against a pre-registered bar of ≥ +0.12 to act and < +0.02 to call it null.
**Cause 2 — gradient noise — is eliminated.** Endpoint medians were the project's best anywhere
(median 0.2159, max 0.2468) **on a null trajectory**, which is precisely the trap §6.3 names; they
are reported and not acted on.

**The hard stop written before the run is honoured.** §6.6 said a null "ends the search — it is not
grounds for a seventh condition." No seventh condition was run on this estimator.

**Status of the implementation: built, selftested, mutation-tested.** The default
path is unchanged — `grad_samples 1` and momentum 0.9/0.5 reproduce every earlier run — so nothing
already reported is affected. Guards added with it, because the free-footprint path taught that an
unexecuted path is probably broken: the selftest asserts momentum is not a literal again, that the
multisample path stays gated on `grad_samples > 1`, that **every keyword passed to
`build_multisample_loss` is one upstream `b94b9d4` accepts** (a bad keyword would otherwise surface
20 minutes into a billed run), that the monomer leg is not routed through it, that the run banner
states the estimator, and that **every flag in every probe launcher exists on the entrypoint**. Each
guard was mutation-tested in the direction it is meant to catch.

---

### 6.7 The epitope probe: the last testable variable, also null

Five conditions and the free-footprint probe had all used the **same epitope**. Three independent
sources said the epitope was the likeliest remaining cause: AlphaProteo's own diagnosis of this
target (§6.5), a reviewer's reading that *"it probably comes down to the epitope"*, and Glogl *et
al.* measuring that this surface class is hard. So the epitope became the variable.

**Region I** — a concave hydrophobic site **15.42 Å** from the receptor site, structurally
unrelated to it, where a published campaign reported 0.55 nM. Positional 69/70/92 on one protomer
plus 109/110 on the adjacent one = mature 74/75/97 + 114/115, with a **2.97 Å seam** between 92.B
and 110.C, so it genuinely spans two protomers as the challenge recommends. Matched to condition 2
with the epitope as the only variable, including keeping the **old** estimator, because arm A had
already returned a null and changing two things at once wastes the run.

**Bar pre-registered before the run**, same statistic and tool: four-run mean rise in
`iptm_repred`, `loss_traj.py --block 13`, ≥ +0.12 to act, < +0.02 null. The R108 null distribution
at that point was seventeen trajectories: +0.004 ± 0.013 (free footprint), +0.008 ± 0.012 (arm A),
and twelve flat tuning trajectories.

**Result: −0.0092 ± 0.0126. Null.** The interface term did not improve, and it moved *down*.

**The hard stop is honoured.** The launcher's own header recorded *"this is the last arm."* It was.
**Twenty-one trajectories, six conditions, two structurally unrelated epitopes.** The interface
term improved in none of them.

🔶 **Independent corroboration arriving after the fact, 2026-10-08.** Genie 3's published
BinderBench ships a TNF-α problem whose hotspot is mature **113 + 73 across two protomers**, tagged
to AlphaProteo. Mapped into its numbering, **four of Region I's five residues land in or adjacent to
its own `common` interface set, one exactly**. So the epitope §6.7 tested is the surface a published
campaign independently selected for this target. That makes §6.7's null likelier to be a
**generator** failure than an **epitope** failure — which is the conclusion §48 of the challenge
notes reaches from a different direction entirely.

---

### 6.8 No computational precedent exists for the pH leg, and no predictor represents it

The three objectives are not equally supported by prior work, and this section states the gap
rather than leaving it implied.

**No structure predictor used anywhere in this submission represents the pH-6.0 state.** Boltz-2,
ESMFold2, OpenFold3 and Protenix take a sequence, and optionally a structure or alignments. **None
takes a pH or a protonation state.** Every interface number we report is therefore computed on a
pose generated with no notion of pH, and the pH objective is a **design-time geometric bet** rather
than an optimised quantity. No amount of co-folding changes this.

**What the two closest prior works actually give us.** Ahn *et al.* (bioRxiv 2025.09.29.678932)
supply the geometric criterion we implement — a histidine within hydrogen-bonding distance of a
cationic partner across the interface — and Schröter *et al.* 2014
(doi:10.4161/19420862.2014.985993) supply the antibody-side precedent for pH-dependent binding.
Neither releases code, and neither reports a de novo design at pH 6.0.

**Our own His-cation contact count sits at the published floor.** Two. Working designs in the
released data carried **8 and 11**. That is the single most direct statement available about
whether the installed mechanism is dense enough to act, and it is unflattering.

**A direction-aware linkage measurement, added 2026-10-08, and it is worse than the count.** With
the TNF-α target registered in `bin/ph_gate_multisite.py --problem 2`, all 35 designs were scored
for thermodynamic linkage, `ratio = K(6.0)/K(7.4)`. Problem 2 binds at 7.4 and must be silent at
6.0, so **acid must weaken binding: ratio < 1**:

| | |
|---|--:|
| moving sites across 35 poses | 219 |
| pointing the **right** way | 80 — **37%** |
| pointing the **wrong** way | 139 — **63%** |
| **no counter-charge within 6 Å** | **167 — 76%** |
| that are histidines | 89 — 41% |

All-site product: min 0.326, median **2.779**, max 1406.8. Only **11 of 35** point the right way as
whole molecules. **The median design binds ~2.8× tighter at pH 6.0 — a reverse switch, gripping
hardest where it is meant to let go.** And 76% of the movement comes from sites with nothing to push
against, the desolvation signature this project named on H370; without that guard, 219 moving sites
would have read as a working mechanism.

**A hypothesis, offered as one.** These designs are acid-rich, and both generators available to us
produce acid-rich sequences. An interface carrying carboxylates on both sides is mildly
self-repelling at pH 7.4; protonate the acids at 6.0, the repulsion eases, and the complex
*tightens*. If that is what is happening, the composition has been building an inverted switch by
accident, and histidine placement cannot fix it while the acids dominate. We have not tested this
and do not claim it.

**Caveat on the whole section.** These linkages are computed on poses whose interfaces clear a
size-matched shuffled null by 2–3% of a real receptor's margin (§4.5). The direction of the finding
is worth something; the magnitudes are not.

---

### 6.9 Germinal: an attempted arm that never ran

Recorded because an attempted-and-failed arm is part of the method. Germinal was wrapped
(`patches/modal_germinal.patch`) and never produced a design — the wrapper required a vendor
dependency that 404s, patched around, and the arm still did not execute. It contributed **nothing**
to this submission and is listed so the generator inventory is complete rather than flattering.

**What the released data says about generator choice, which is the comparison §6.5 should have
drawn.** On the organisers' own 150 measured TNF-α designs:

| generator | binders / tested | rate |
|---|--:|--:|
| **Genie3** | **8 / 23** | **34.8%** |
| PXDesign | 4 / 57 | 7.0% |
| RFdiffusion3 | 0 / 40 | 0% |
| RFdiffusion | 0 / 20 | 0% |
| BoltzGen | 0 / 8 | 0% |
| FreeBindCraft | 0 / 2 | 0% |

Fisher exact, two-sided, Genie3 against all others pooled: **p = 0.00003**. Genie3 is also the
**only** generator in that set that produced mouse cross-reactivity (3 of its 8). And **all 12
binders used SolubleMPNN**; ProteinMPNN and Caliby produced zero. Two caveats that matter: the 8
are variants of **one** backbone found in 23 samples (pairwise identity 0.58–0.81 inside a set
spanning 0.02–0.90), so this is one lucky discovery rather than eight independent ones; and the
0/2 for FreeBindCraft is **almost no evidence** — a true 7% rate yields 0/2 about 86% of the time,
whereas 0/40 and 0/20 are genuinely informative.

**We used neither Genie3 nor SolubleMPNN for the designs in this submission.** That is the single
largest identifiable gap between this method and the one that worked on this target, and it was
visible in data published before our search closed.

---

## 7. Protomer count is a measured confound for any inter-protomer epitope

The challenge recommends an epitope spanning two protomers. Any ipSAE-type score computed for such
a design depends on how many protomers were modelled, and the dependence is large.

**Control.** Adalimumab — approved, binding exactly this site — folded against one, two and three
protomers. Fab chains and the TNF protomer lifted verbatim from the three-protomer construct, so
protomer count is the only variable. Five seeds each.

| target | ipSAE_min | relative to trimer |
|---|---|---|
| **1 protomer** | **0.0322** | **−94%** |
| 2 protomers | 0.4972 | −14% |
| 3 protomers | 0.5757 | — |

Per-seed on the monomer, with the Fab's own VH/VL interface as an internal control:

| seed | best TNF:Fab | Fab VH:VL |
|---|---|---|
| 1 | 0.0459 | 0.8586 |
| 2 | 0.0131 | 0.8395 |
| 3 | 0.1407 | 0.7931 |
| 4 | 0.0133 | 0.8187 |
| 5 | **0.0000** | 0.8434 |

**The Fab is intact in every seed while its interface with the monomer reads 0.000–0.141.** This is
not a folding failure: the binder builds correctly and has nothing to bind.

**An 18× collapse on an approved antibody.** The transferable statement is that **an ipSAE-type
number for an inter-protomer design is uninterpretable without stating the protomer count it was
computed against** — and single-chain targets are the default in most binder pipelines. A team
modelling the recommended epitope against one chain would conclude its designs did not bind.

**A dimer recovers most of the signal** (0.4972, 86%), which is what makes the two-protomer species
leg in §6 defensible. One of five dimer seeds failed to assemble (internal packing 0.067 against
0.398/0.376/0.375/0.401), so dimer runs carry roughly 20% attrition to non-assembly, and an
assembled dimer packs about half as tightly as the trimer.

---

## 8. Filtering and selection

*Method fixed and exercised; no design has passed it.*

**The gate has now been applied to all 21 candidate sequences and rejected all 21** (§10). That is
a result about the filter, not only about the designs, and it cuts in an uncomfortable direction:
**a gate that has only ever rejected is half-validated.** Its specificity was set from the negative
side, where measurements exist, and every application so far has been a rejection — so nothing in
this work establishes that it correctly *admits* a true conditional binder. The false-negative rate
is unmeasured and unmeasurable here, because no de novo TNF-α binder exists to anchor the upper
side. A reader should treat `n/a` on all 21 as evidence the designs did not bind, **and** as a
reminder that the instrument has never been shown to say yes.

**Geometry reporting is gated on binding.** The distance from the binder histidine to the target
cation is reported only when the re-predicted complex reaches **iptm ≥ 0.45**; below that it is
recorded as **not applicable**, never as a failure. The bar is set from the negative side, where
measurements exist: a confirmed dead fold reads 0.408–0.417 across five seeds, and 5-step
non-binders read 0.11–0.31, while the measured G1 binders read 0.77–0.86 across twenty. It is a
**necessary condition only** — deliberately not set at 0.77, because those are Fabs and a correctly
folded non-binder reads 0.46 (§5, G4b). No measured de novo binder against TNF-α exists to anchor
the upper side.

Applied retrospectively, this reclassifies every geometry verdict produced during development as
not applicable, including one 3.15 Å "pass" that came from a sequence carrying 12 histidines, where
one lands near the cation by chance.

**Liability screening.** `bin/express_qc_p2.py`. The specific hazard for this challenge: installed
histidines raise the binder's isoelectric point, and the assay is run at **both** pH 7.4 and
pH 6.0. A design that aggregates at pH 6.0 reports as "no binding at pH 6.0", which is
indistinguishable from a successful switch. This is the one screen that can manufacture a false
positive on the primary objective.

**Declared developability tension.** Raising free pKa raises pI, and elevated pI is associated with
faster non-specific clearance — Igawa et al. measured +2 pI units giving roughly 2× shorter
half-life, with a later study finding a bell-shaped relationship where charge patches of either
sign accelerate clearance (DOI 10.1080/19420862.2021.1993769). **The property that maximises the
switch is the property that degrades developability.** Not scored by this challenge, which measures
binding only.

---

### 8.5 Two filters re-scored against measured outcomes

**The 0.45 interface bar is inert, not strict.** Against the 150 assayed designs it admits **150 of
150** on 8 of the 10 co-folding models — every measured non-binder included — for a precision of
**8.0%**, which is exactly the base rate. It retained 12 of 12 binders only because it rejects
nothing. **It is therefore not a filter and is not described as one.** It survives solely as the
threshold below which we withdraw a geometry verdict to `n/a` (D-P2-1), which is a reporting rule.

**The filter that does work is monomer foldability, and this work never applied it.**
`reference/anthropic-binder-design-protocol.md:100` specifies a default floor of **0.70 mean pLDDT
on the binder alone**, to be run *before any co-folding spend*. We applied the 0.45 bar twenty-one
times and this floor zero times.

| `plddt_binder`, Boltz-2, 1 binder : 3 protomers | 150 assayed designs | our 35 |
|---|---|---|
| minimum | **68.9** | 32.1 |
| median | 85.3 | 47.0 |
| maximum | 93.0 | **68.8** |
| clearing the 0.70 floor | **146 / 150** | **0 / 35** |

**Our best design sits one-tenth of a point below the worst of 150 designs that were synthesised and
assayed.** Two caveats, stated rather than buried: one measured binder reads 69.4, so 0.70 is not an
absolute wall; and their co-folding configuration is not byte-identical to our re-prediction, so this
is order-of-magnitude placement rather than a matched measurement — though the gap is roughly three
times any plausible configuration difference. The floor is now declared and reported by
`analysis/02-tnf/design_inventory.py`, which refuses to export if `plddt_binder_repred` mixes the
0–1 and 0–100 scales across runs.

**Pre-submission screen.** `adaptyvbio/binder-prescreen`, the sequence-based screen that replaced the
novelty gate for TNF-α on 2026-10-08, run on all 35 on 2026-10-07: **35/35 pass**, whole-sequence
similarity to the 4,964-entry known-anti-TNF corpus **0.000 for every design**, binding-region
identity at most **0.260** (median 0.229). The public prior-art arm was not run — it requires the
10.2 M-sequence patent database — so this is the TNF-specific arm, which is the arm the submission
gate uses.

## 9. Controls

*Design fixed; results pending generation — and their value is now conditional, see below.*

🔴 **A matched pair is only informative if the parent binds.** With 0 of 21 candidates clearing the
interface gate (§8, §10), a histidine-removal control would compare two sequences that both fail to
bind, and the paired difference would measure nothing. **The controls do not rescue a non-binding
set; they only interpret a binding one.** This is stated here rather than discovered at analysis
time, because the temptation to report a paired difference between two non-binders as a mechanism
result is exactly the failure the preceding challenge produced — a 1.432× ratio that dissolved
under matched arms.

Their design is unchanged and is kept for the case where generation succeeds. If it does not, the
four control slots are better spent on additional candidates, and the submission says so rather
than shipping controls for an effect that was never established.

**Matched histidine-removal pairs.** Four of the twenty submitted designs are derived from shipped
backbones by removing the engineered histidine and re-predicting. The backbone is otherwise
identical, so the pair differs in the engineered feature and nothing else.

**These controls carry two independent functions.**

**(a) They test whether the engineered histidine does anything.** In the preceding challenge, a
1.432× ratio appeared to be a pH switch until matched arms pinning a non-titratable residue
switched at the same rate — 4% versus 6%, p = 0.53 — and a carboxylate-free patch switched at the
same rate as the acidic sites, 5% versus 4%, p = 0.74. Without the matched arm, a ratio in that
range is not evidence of a mechanism.

**(b) They cancel a target-stability confound.** Published kinetics give TNF-α a trimer-to-monomer
rate of **k = 1.66 ± 0.25 × 10⁻³ s⁻¹** — a ~7 minute half-life, complete within 30 minutes — with
trimer formation setting in only above roughly **10 nM** free TNF-α, and the dissociation rate
reported as significantly affected by solution pH (*Sci Rep* 10:9477, PMID 32518229; SPR kinetics
PMID 8186365). Two consequences follow: an inter-protomer epitope may be sparsely present at low
analyte concentration, and if the trimer dissociates faster at pH 6.0 than at 7.4, a loss of signal
at 6.0 is a property of the target rather than of the design.

**A matched pair run in the same assay experiences the same target dissociation.** The effect is
therefore common-mode and cancels in the paired comparison. **This work reports the paired
difference as its mechanism evidence, not an absolute pH ratio** — the paired difference is the
quantity the confound cannot reach.

**What a zero-histidine design is not.** It is not necessarily pH-insensitive: other titratable
groups and the target itself still contribute. These are described as tests of the engineered
histidine's contribution, and each pair is additionally checked for whether the substitution has
simply disrupted binding, which would make the comparison uninformative.

🔴 **The surviving limitation is pairing rather than chemistry.** If only one member of a pair is
synthesised, the primary evidence is lost. The pairs are identified explicitly in the submission
with a request that they be considered together.

**A control that would settle the confound at the assay level**, offered for interpretation of the
published dataset rather than as a request: running the pH 6.0 arm at two analyte concentrations
spanning the ~10 nM onset. A design whose signal loss is real gives the same ratio at both; one
riding on target dissociation does not.

---

## 10. The submission

*Pending generation.* The allocation and the ordering rule are fixed.

| arm | designs | anchor |
|---|---|---|
| R108-oriented | 12 | canonical R108 = our positional 27 |
| K166-oriented | 4 | canonical K166 = our positional 85 |
| matched histidine-removal controls | 4 | derived from shipped backbones |

At most two candidate slots may come from the inverted-objective pilot (§12), and only if it
produces a defensible result.

**Backbone and footprint diversity is an explicit axis, not a by-product.** Candidates are spread
across two footprint families: one with the epitope pinned to the nine conserved positions, and one
with **no epitope restriction**, where the optimiser selects the footprint subject only to the
anchor contact. The second is motivated by a published de novo result in which the epitope was not
specified at all and the optimiser's choice outperformed hand-chosen functional sites. Varying
length or seed does not produce distinct families and is not treated as diversity.

**Ordering.** Eligibility first, then structural plausibility, then liabilities and robustness,
while preserving family coverage and keeping matched pairs adjacent. **The ordering is provisional,
for experimental testing.** No retired score is used as a calibrated prediction.

### What exists right now, counted rather than estimated

Candidate material is not hypothetical, and stating its size is more useful than describing the
intended allocation. **Twenty-five distinct sequences exist**, one per trajectory across every
condition run, at two lengths. That is more than the twenty slots. **None of them is a success by
this work's own instruments:**

| | count of 25 |
|---|---|
| sequences, all distinct | 25 |
| clearing the 0.45 interface gate | **0** |
| therefore carrying a geometry verdict other than `n/a` | **0** |
| histidine fraction above the 8% cap | 3 |
| nearest histidine within the 4.0 Å criterion | 2 |
| of those 2, defensible | **0** — see below |

**Both sub-4 Å placements fail for separate reasons, and neither is reported as a result.** One
measures 2.62 Å but sits at an interface confidence of 0.200, so its verdict is `n/a` and the
distance is an annotation, not a verified placement — a precise distance inside a complex the
predictor does not assert exists. The other measures 3.15 Å and is invalid twice over: the sequence
carries 12 histidines in 76 residues, which is the composition pathology of §6.1, and it predates
the protomer fix of §11, so the measurement itself came from the code that searched the wrong
protomer. It appears in §11 as a withdrawn pass and nowhere else.

**The consequence for the submission is stated plainly.** A submitted set drawn from this material
would consist of sequences this work declines to validate, accompanied by the reasons it declines.
That is a worse experimental bet than a set of verified switches and a better one than a set
presented as verified switches when it is not. The instruments that would have let us claim
otherwise are the four retired in §4 and the gate in §5, and they were built and tested before
these designs existed, which is why they can refuse them.

**Diversity is now real but does not help.** 21 of the 25 belong to the pinned footprint family and
4 to the free-footprint family (§6.4), so the two distinct footprint families an external
reviewer asked for do exist. They do not reduce
the risk §14 names, because the second family scores *worse* than the first on both the interface
term and histidine placement. Two families of non-binders is better coverage of a space that does
not contain a binder.

---

## 11. Errors found in this work, and what they cost

Recorded because each one changed a method, and because several were invisible to the tests
written to catch them.

**A geometry check measured the wrong protomer.** TNF-α is a homotrimer: after renumbering, all
three chains carry `ARG 27`. The check searched all target chains for the first arginine at the
anchor residue number, with no knowledge of which protomer the loss was aimed at. With the anchor
declared on chain B, the loss pulled toward protomer B while the check measured protomer A — a
histidine placed perfectly reads ~30 Å. **Every geometry number produced before the fix carries
this defect.** A prior comment in the code asserted this had already been fixed by searching all
chains rather than one; searching all chains is the defect, because on a homotrimer the matches are
indistinguishable.

**The test written as that fix was vacuous.** It planted decoy protomers 60 Å away but named them
alanine, so the arginine/lysine filter disambiguated a search that is ambiguous in production. It
asserted a property of a fixture easier than the real input, and passed while the production path
was wrong. **A fixture easier than production is worse than no test, because it is reported as
coverage.** The replacement makes all three protomers arginine and asserts that every *wrong*
protomer index reports a miss; reverting the fix now makes the test fail.

**A loss term bought histidine quantity instead of histidine placement.** The pH term summed its
per-position scores, which rewards total histidine mass near the cation rather than one histidine
positioned. The cheapest path down that gradient is histidine everywhere: a smoke run returned
**12 histidines in 76 residues — 15.8% against a ~2.3% natural frequency** — while the diagnostic
reporting the best single placement read exactly 0.00 at every step. The incentive is provable on
the sigmoid without compute: 12 histidines at 9 Å score −0.646 and one histidine at 3 Å scores
+0.091, so the summed form strictly prefers the spam. Reduction changed to **top-2** — two rather
than one because the computed 0.40 units per site requires two sites. After the fix, the same
configuration returned **2 histidines, 2.6%**, and the count equalling *k* is what identifies the
cause: the composition cap is a soft hinge at six residues, so a cap-bound sequence would sit near
six.

**Two guards were inert on this challenge's entire mechanism.** The pass-versus-distance check and
a bounds check both read the dict key belonging to the *other* mechanism, so on every design in
this challenge they found nothing and skipped. They passed while completely blind.

**A single histidine-free design would have destroyed a whole run's output.** The results table
took its columns from the first row and indexed rather than looked up, while the geometry function
returns three different key shapes. The first heterogeneous pair raises `KeyError` — before the
only write to durable storage. Both fixed; designs now persist individually as each finishes, which
has since been validated by a network failure that killed a run and left the completed design
intact.

**A control construct was void, and the error had a five-file blast radius.** The first
lymphotoxin-α control folded at pLDDT 0.27–0.39 with every inter-chain score at exactly 0.0000 —
including the trimer's own internal packing, where human TNF-α reads 0.71–0.74. The construct
carried 37 residues of signal peptide and was missing ~24 residues internally, having been built
from a sequence alignment rather than from the solved chain. Rebuilt from the observed 1TNR chain,
it folds at 0.86 and packs at 0.79. An audit of every construct against every locally available
solved structure then found the same broken sequence in **five files**, all now quarantined.
`analysis/02-tnf/construct_audit.py` enforces the check and currently reports no live construct
carrying a signal peptide.

**Instrument scope was overstated and corrected.** Several documents described computed pKa shifts
as "measured" and described a 10× ratio as a challenge requirement. Neither was correct: the pKa
values are PROPKA output on our own poses, and the published objective names no ratio.

---

**Two summary cells in this document reported a row value as a median.** The condition comparison
in §6.2 gave the baseline median as 0.153 and the probe median as 0.183. Recomputed from the design
tables, the medians are **0.136 and 0.179**. Both wrong figures are real numbers from the tables —
0.1527 is the baseline's first row, 0.183 is the probe's third — so a cell meant to hold a statistic
had been filled with the row above it, twice, in the same table. The derived claim moved too: the
effect of the intervention is **+31%**, not the +20% stated. Caught by recomputing every summary
statistic in the document from the source tables rather than by reading the prose, and the direction
of the conclusion is unchanged. **Numbers that are individually real are the hardest transcription
errors to see**, because every spot-check of a cell against the data succeeds.

**A flat optimisation curve was not looked at until four conditions had been run.** Three weight
conditions were tested and compared on their endpoints before anyone plotted the interface term
against step number. The plot (§6.3) shows the term never improving in any of them, which would
have redirected the search after the first condition instead of the fourth. **An endpoint comparison
cannot distinguish a method that is being tuned badly from a method that is not optimising**, and
only one of those is worth more tuning. Roughly $15 of the $26 spent on the search went to a
question the first run's own log could have answered.

**The same error class recurred within hours, and the guard built that morning caught it.** The
condition table gained a fifth row, and its closest-histidine cell was filled with the condition's
**median** (34.8 Å) in a column that holds the **minimum** (20.1 Å). This is the identical mistake
as the two median cells above — a real statistic of the right data placed in a cell that holds a
different statistic — committed by the same author on the same day, hours after writing the entry
describing it. It was caught in seconds by `bin/check_p2_stats.py`, which recomputes every cell in
that table from the design tables and had been written specifically because the first instance was
invisible to proofreading.

**Two things follow, and the second is the useful one.** First, knowing about an error class does
not prevent committing it; the prose entry above was written by someone who then made the mistake
again. Second, **the value of a mechanical check is not that it finds errors a careful reader would
miss — it is that it does not get tired or confident.** The check cost roughly twenty minutes to
write and has now caught the thing it was built for, once, in the same session.

**An external corroboration was overstated, and half of it was withdrawn an hour later.** Reading
the released data of the campaign in §1, this work first reported that *two* archived structures
docked on our chosen epitope were binders of that campaign, and cited both as support for the epitope
choice. Both halves were wrong. The release's own README states that the archived files are **input
scaffolds**, not predictions of any assayed molecule, and that binder side chains may be `UNK`
placeholders; and the lineage table shows the parent of the repository's flagship TNF-α example,
`TNFA_UNOPT_04`, **was tested and did not bind**. So one of the two footprints cited as support is a
measured non-binder's scaffold. The claim was withdrawn and §2 now states what survives — one
scaffold, on our epitope, whose optimized child bound at 66.2 nM while retaining 60–70% of its
contacts across seven predictors — together with the counterexample.

**Why this one is worth the space.** The error was not caught by a test or by re-reading the prose. It
was caught by downloading the primary data release and checking a claim this work had already written
down, which took about fifteen minutes and reversed a conclusion. **A corroboration is the easiest
kind of finding to accept without checking, because it agrees with a decision already made** — this
one arrived on a morning when the method had just failed, and it was welcome. The surviving claim is
weaker than the one first written, and the target corroboration of §1, which is clean, was never in
doubt either way.

**Two correct decisions left one design oracle, and it was the one we had excluded.** On 2026-09-18
Boltz-2 failed this project's own validation gate as a ranking arm: with a staged target MSA a
shuffled negative scored **0.7326 against the true binder's 0.6463**, and the recorded diagnosis was
that Boltz needs alignment context on both chains to judge an interface, *which de novo binder design
can never supply*. Verdict: exclude Boltz. On 2026-09-29 a frozen decision forbade designing on the
ranking instrument, to keep the ranking question answerable in December. Both decisions were correct.
But the third permitted predictor emits no PAE, so the two together left **exactly one permitted
design oracle: the excluded one.** Twenty-one trajectories descended its gradient, and the external
calibration in §4.5 puts it last among credible models at AUC 0.781. **The project had a rule for
this** — an arm that fails its gate is excluded for that target — **and applied it to ranking but not
to generation, because the rule governed what a design objective may not use and never what it must
pass.** This is recorded as the leading explanation for the flat interface term, displacing the
term-interaction account in §6.3, which is retained as a contributing factor but no longer the
headline; it never explained why the monomer leg, which carries no interface or pH term, also folded
nothing.

**A number corrected in one document was still wrong in another, for the second time.** §11 corrected
a +20% effect size to **+31%** on 2026-10-07 at 7:35 AM, and the public methodology box was updated
the same morning — in one of its two prose fields. The other field still read 20% at 11:30 PM that
night. The lesson recorded the first time was *"a number fixed in one document is not fixed"*; the
recurrence shows that lesson was recorded and not operationalised. There is still no check that greps
every published surface for a corrected figure.

## 12. The inverted-objective pilot

*Time-boxed; result pending, and reportable either way.*

Proton-PottsMPNN was released during this competition and recommended in the competition channel
as relevant to pH-sensitive design. As shipped, its selective term minimises `e_P − e_D`, which
favours the protonated state and therefore designs for binding at **low** pH — the opposite of this
challenge's direction.

Three observations from reading the source. The relative weight is an unclamped parameter, so a
negative value inverts the selectivity term. The acid-dominated behaviour reported for the method
follows from a default list of titratable centre types, not from the method, and restricting
centres to histidine makes the contrast protonated-versus-neutral histidine. And the method is
inverse folding, so it composes with gradient-based backbone design rather than competing with it.

**What is not established, and governs the design of the pilot.** The published validation concerns
acidic-pH binding to a different target and does not validate the inverse objective. A negative
weight changes both the selectivity sign *and* the stability coefficient, so it is not a
sign-only comparison. And the candidate-placement step ranks positions by a protonated-state
preference that is **separate from the redesign term** — so inverting the weight alone could leave
placement still selecting for low-pH binding.

**Protocol.** Matched baseline, original-direction and inverted-direction runs on the same de novo
backbones, with the direction verified independently of the objective being optimised. The pilot
stops at half a day if the comparison is not interpretable or reproducible. A negative result is
reported. Method novelty alone does not displace a better-supported candidate.

🔴 **This pilot inherits the same precondition as §9, and the inheritance is the point.** It
redesigns sequences onto de novo backbones and compares three directions on them. With 0 of 21
backbones clearing the interface gate, all three arms would be built on scaffolds that do not bind,
and the comparison would resolve a question about pH selectivity on complexes that have no measured
affinity to be selective about.

**Stated once, because it applies to three sections at the same time.** §9's controls, §10's
submission and this pilot all assume a binding design upstream of them, and §6.3 reports that the
generation method did not produce one. **The failure is upstream of every downstream instrument in
this work**, which is why those instruments read `n/a` rather than `False`, and why this document
reports the generation result as the finding rather than burying it as a caveat under each
dependent section. The instruments are sound and unexercised in the direction that matters; the
generator is the thing that did not work.

---

## 13. Reproducibility, measured in a fresh clone

Claims in this section were tested by cloning the published repository into an empty directory
with no local state and running the tests there. **Two of the four did not work for any reader
when first tested**, and both faults are described below rather than silently fixed.

**MEASURED 2026-10-07, fresh clone, system `python3`, no virtualenv, no GPU, no Modal account:**

| | |
|---|---|
| `bin/mosaic_selftest.py` | **PASS** |
| `bin/esmfold2_selftest.py` | **PASS** |
| `bin/plan_wave.py --selftest` | **PASS** |
| `analysis/02-tnf/construct_audit.py --selftest` | **PASS** |

One preparation step is required and is the reason the first two initially failed:

```
mkdir -p biomodals
git apply --directory=biomodals patches/modal_mosaic.patch
git apply --directory=biomodals patches/modal_esmfold2.patch
```

### The two faults, because they are the same fault problem 1 had

**The wrappers are not published.** `biomodals/` is gitignored — the Modal wrappers are large,
change independently of the analysis, and are carried as patches under `patches/`. Both wrapper
self-tests read the wrapper file directly, so in a clone they failed with `FileNotFoundError`
**while passing locally**. The tests that gate both launchers could not run for any reader.

**And one patch could not reconstruct its file.** `modal_mosaic.patch` was generated as a creation
patch (`new file mode`) and reconstructs cleanly. `modal_esmfold2.patch` was a *modification*
diff (`index e3ea6da..ddf8b5e`), which requires the original file to already exist — in a clone
where `biomodals/` never exists, it can never apply. It is now a full-file creation patch.

**A second unpublished dependency, in the same test.** The esmfold2 naming check compares output
filenames against output the pre-fix code really wrote, which lived in `runs/` — gitignored and
tens of gigabytes. The check needs only the **names**, so 25 of them are vendored at
`outbox/esmfold2-fixtures/reference-names.txt` and the check uses the local cache when present
and the vendored list otherwise. The one assertion that genuinely needs the files — that each
`*_ipsae.json` has its `.cif` beside it — is skipped without the cache, and says so in its output
rather than passing silently.

**This is the same class of fault problem 1 documented**, where a scorer looked for its reference
only at a gitignored path and therefore failed closed for every reader while passing locally. It
recurred in a different file. The lesson that transfers is narrow and mechanical: **a test that
reads a gitignored path is not a published test**, and the only way to find out is to clone and
run, not to reason about it.

### What the tests themselves guarantee

- Every instrument's self-test runs with no GPU, no Modal account and no network, and both
  launchers refuse to start if theirs fails.
- Self-tests include **mutation tests**, which assert the *discarded* implementation fails — so a
  test cannot pass if the fix it guards is reverted. Verified by reverting: the protomer fix, the
  top-2 reduction, and the CLI weight coverage each make their test fail by name when undone.
- Numbering is checked against **residue identity** on every run, not against arithmetic.
- The `ipSAE_min` pre-registration is a commit timestamped **2026-10-06 16:24:41**, preceding the
  scoring run, so the ordering of registration and measurement is checkable rather than asserted.

### Which code produced which number

Every number in §6 was produced by a different build of the same wrapper, because the wrapper was
being fixed on the same day the conditions were run. That makes provenance a claim this document
has to support rather than assume, so it was checked mechanically rather than from memory.

**The last change with any effect on a result landed before the first reported run.** Four
behavioural fixes — the protomer-selecting geometry check, the top-2 histidine reduction, the
cation-atom guard moved out of the jitted loss, and the timeout safety factor — were all committed
by 20:15. The earliest run contributing a reported number started about three hours later. Four
earlier runs exist and **none of them supplies a figure in this document**; they were smoke and
configuration tests, and they are excluded for that reason rather than by argument.

**One commit does fall inside the reporting window**, between the baseline condition and the three
that followed it. Since that boundary separates condition 1 from conditions 2–4, and the +31%
figure in §6.2 is a comparison across it, it was checked rather than assumed. It changed one
function, in one statement: the call forwarding six newly CLI-settable weights. Each new default
was then compared against the weights the baseline run actually recorded — all six identical. The
forwarded values are the values that were already in use, so the comparison spans no behavioural
change.

**The free-footprint fixes do not reach any reported number either.** They were made after every
figure above was recorded, and the method used to establish that is worth stating because reading
the diff is not sufficient: the pre-fix and post-fix sources were reconstructed from their
committed patches and compared function by function at the level of parsed syntax, not text. Of 18
functions, 17 are identical — including every scoring, geometry and gating function. The one that
differs reduces to **ten minimal differing syntax nodes**: two display strings, two error-message
strings on a path that never executed, one guard whose outcome is unchanged for the inputs used,
and the parsing statements, which were then tested for equivalence across every epitope string any
run has passed. None lies inside the loss, the optimiser call, the metric computation or the
output writer.

**This is a weaker claim than it looks, and the limit is worth naming.** It establishes that the
code paths differ only in display, not that the runs are numerically reproducible — they are not,
for the reasons below. A trajectory re-run today would not return the same sequence.

### A dependency was unpinned, and the claim was wider than the evidence

**Found 2026-10-07 while evaluating an external method, not by a test.** The Modal image installed
Mosaic as `uv_pip_install("git+https://github.com/escalante-bio/mosaic")` with **no revision**. A
fresh clone therefore resolved the design dependency to whatever that repository's HEAD was on the
day it ran, which makes "reproducible in a fresh clone" a weaker statement than it reads as — the
tests above pin nothing about the optimizer that produced the design numbers.

**What the exposure actually was, measured rather than assumed.** The commit the wrapper's API was
read against is `b94b9d4` (2026-09-24). Only two Mosaic commits landed between that and the end of
this project's run window, `59492d4` and `1eb4af6`, both on 2026-10-06, and the compare across them
is **+21 / −0 lines across `CITATION.cff` and `README.md`** — measured on the GitHub compare API,
with **no Python file touched**. So every run in this project executed the same optimizer code, and
none of the comparisons in §6.1–§6.5 is confounded by a dependency moving underneath it. That is luck
rather than method: nothing in the setup would have revealed it if a commit had changed the
optimizer mid-campaign.

**Fixed.** The install is now pinned to the full revision
`b94b9d4eb9907a700a6d78ed2d29d3704c5df46c`, chosen because it is both the commit the wrapper was read
against and the code every run in this project actually ran — so the pin reproduces the campaign
rather than freezing it at an arbitrary later point. The change is in `patches/modal_mosaic.patch`,
which is the tracked artifact a reader applies; `biomodals/` is gitignored, so a pin landing only in
the working file would not have reached anyone. **`bin/check-pins.sh` now asserts the revision in both
places**, and the guard was mutation-tested in each direction — stripping the pin from the working
file alone, and from the patch alone, each turns it red.

### What is not reproducible

- **Pose caches** under `runs/` are gitignored and run to tens of gigabytes; available on request.
- **The design runs themselves** require a GPU and a Modal account. The wrappers, parameters and
  weights are published; the compute is not reproducible without equivalent hardware.
- **The platform's novelty checker** changed during the competition and its implementation is not
  published, so its verdicts cannot be reproduced locally.

## 14. Limitations

Grouped by what each one threatens. Several are consequences of decisions made deliberately, and
those are marked.

### The mechanism

1. **The pKa shifts are computed, not measured.** PROPKA 3.5.1 on our own poses. The 2.8–3.9×
   two-site estimate inherits that status entirely, and no experimental pKa was obtained for any
   design in this work.
2. **PROPKA's own accuracy bounds the estimate.** Published RMSD for histidine is of the order of
   ±0.8 pKa units. A 0.40-unit computed shift is inside that envelope, so the sign of the effect is
   better supported than its magnitude.
3. **The per-site ceiling is thermodynamic and cannot be designed around.** Over the 7.4 → 6.0
   window a single titratable site cannot exceed ~25× regardless of geometry, and reaches 1.92× at
   a free pKa of 6.0. No amount of interface optimisation moves that bound.
4. **Two sites are assumed independent.** The linkage arithmetic treats the sites as
   non-interacting. Coupled sites would give less than the product, and we have not measured
   coupling.
5. **A zero-histidine control is not pH-inert.** Other titratable groups and the target itself
   contribute. The controls bound the *engineered* contribution, not the total.

### The instruments

6. **n = 4 on the affinity ladder.** Sufficient to retire a prespecified use of the metric,
   insufficient to characterise the metric. ρ = 0.400 is also not distinguishable from chance.
7. **The four constructs in that ladder are not independent.** They are one antibody and three of
   its engineered variants, so the effective sample size is below four.
8. **The labelled-negative control is one pair, repeated.** Five seeds of adalimumab against LT-α
   are five predictions of the same two molecules, not five biological controls. It cannot
   establish sensitivity or specificity for de novo binders.
9. **No measured de novo binder against TNF-α exists** to anchor a confidence threshold, so the
   iptm gate is a necessary condition only and its upper side is unanchored.
10. **Design-time interface scores are systematically optimistic on the best-looking runs.**
    Measured translation from design-time to re-predicted: +17%, −0%, −49% and −60% across four
    cases, with the two largest drops on the two highest design-time scores. Any figure quoted
    from inside an optimisation run should be treated as an upper bound.
11. **The instruments we retired were tested on histidine point-variants**, not on de novo
    minibinders. The retirement is valid for our use; transfer to other formats is not established.

### The target and the assay

12. **The target omits five residues.** Mature 1–5 (`VRSSS`) are disordered and absent from the
    crystal. Designs are made against 152 of 157 residues.
13. **One sidechain is rebuilt.** Mature 143 was point-mutated from Leu to Asp with PDBFixer to
    match the canonical sequence. The backbone is crystallographic; that sidechain is modelled.
14. **Target stability at assay concentrations is unconfirmed.** Published kinetics give a ~7 minute
    trimer half-life with trimer formation setting in only above ~10 nM. *Deliberate:* handled by
    the matched pairs rather than assumed away, with the residual exposure that both members of a
    pair must be synthesised.
15. **The assay's detection limits are unconfirmed for this challenge.** "No detectable binding"
    depends on them, and the limits referenced elsewhere in this project were carried over from a
    different target and are not assumed here.
16. **Avidity is uncharacterised.** A trimeric analyte over immobilised binders can engage more
    than one site, raising apparent affinity and compressing the pH ratio being ranked.
17. **Buffer composition at the acidic point is unconfirmed** for this target, and ionic strength
    affects an electrostatically driven mechanism directly.
18. **Mouse pH coverage is unspecified in the published objective**, which names pH 7.4 only for
    objective 2. A design that releases at 6.0 may be assessed differently depending on whether
    mouse is measured at both points.

19. **Four published methods measured zero binders on this target.** AlphaProteo (0 of 54 tested),
    BoltzGen, PXDesign and Complexa each tested TNF-α designs and obtained none, and AlphaProteo
    attributes its failure to *"a flat, highly polar binding site at an interface between 2 subunits
    in a homotrimer"* — a description of the epitope chosen here (§6.5). The realistic calibration for
    this target is the organisers' own **12 of 150 = 8%**, not any higher figure. This limitation is
    about the target, not this method, and it is the leading external explanation for §6.3–§6.4.
20. **The external epitope corroboration rests on predicted poses, and the same surface hosts a
    measured non-binder.** The scaffold footprint in §2 was computed from archived *input*
    coordinates, every pose involved is predicted, the source states its experiments do not establish
    that binding occurs at the predicted site, and `TNFA_UNOPT_04` — whose scaffold also sits on this
    epitope — was tested and did not bind. The epitope is corroborated as an independently chosen,
    binder-producing site and **not** as a measured binding epitope.
21. **The external target corroboration comes from a campaign that designed against predicted
    structures**, where this work locked a sequence-corrected crystal trimer. Both carry Asp143, so
    the gap is six disordered N-terminal residues plus a structure-source difference that was recorded
    rather than reconciled.

### The design method

22. **The species leg uses two protomers rather than three.** *Deliberate, forced:* three exhausts
    a 48 GB GPU. A control establishes two recover 86% of the three-protomer signal, with ~20% seed
    attrition to non-assembly.
23. **Interface confidence has not exceeded 0.200 in this work**, against 0.45 as the threshold at
    which a geometry verdict is reported at all. No design has reached a value at which placement
    is treated as verified.
24. **Weight tuning appears to have a ceiling.** Three weight conditions spanning a sixfold range
    in the pH-to-binding ratio give medians of 0.136, 0.179 and 0.156. The best is the middle
    condition, which indicates an optimum rather than a direction. Quadrupling the step budget at
    that optimum gives 0.163 — a fourth condition inside the same band.
25. **The weights interact.** Tripling the binding weights pushed histidine content to 9.2%
    against an 8% cap and moved the closest histidine placement from 2.62 Å to 9.69 Å. The loss
    terms are not independent knobs.
26. **The composition cap is a soft hinge, not a constraint.** It has been exceeded. A design can
    ship above its nominal cap.
27. **The pH term is numerically flat beyond ~30 Å**, so it supplies no gradient at the distances
    trajectories begin from. The mechanism depends on the contact terms succeeding first.
28. **The pH term is measured CA-to-cation**, because design-time features give the binder no
    sidechains. The all-atom distance is only available on re-prediction.
29. **Seed variance is large relative to the effects being measured.** Four trajectories at
    identical settings spanned 0.111–0.156, and across the eight runs at the best weight setting
    seed and length alone span 0.115–0.200 — a range about forty times the measured effect of
    quadrupling the step budget. Differences between conditions of that order are not resolvable
    at n = 4, and every condition difference reported in §6 is of that order.
30. **The epitope is restricted to nine positions conserved to mouse.** *Deliberate:* four
    positions that differ were dropped, including one deleted in mouse, which forgoes whatever
    affinity those contacts offered on human.
31. **A single epitope carries the whole submission.** All candidates target one site, so a wrong
    epitope choice fails the set at once rather than independently. This is the largest correlated
    risk in the submission.
32. **Cysteine is excluded from the alphabet**, so no design can use a disulfide for stability.
33. **The interface objective did not improve in any trajectory run.** Twelve trajectories, three
    weight conditions, two step budgets: the design-time interface term is flat from the first
    block of steps to the last (§6.3). The generation method as configured is not optimising the
    quantity it is written to optimise, and **no result in this document should be read as
    evidence that gradient descent on this loss produces interfaces.**
34. **The cause of that flatness is unidentified, and only one of three candidates was
    actually tested.** The pinned epitope admitting no gradient path was tested directly and
    **ruled out**: removing the restriction on both legs changed nothing (§6.4). Sampling noise
    swamping the gradient and interaction among the eleven weighted terms both survive. *An earlier
    version of this limitation called both untestable within the remaining time; for sampling noise
    that was wrong* — the dependency exposes four-sample gradient averaging directly and momentum is
    a literal in our own wrapper, so **the configuration that every one of the sixteen trajectories
    shared (one sample per gradient, momentum 0.9) was never varied, and could have been.** The
    method fails for a reason this work cannot name, and one of the three candidate reasons was left
    untested by an unchecked assumption rather than by the budget. An arm that varies it is now
    built and its bar pre-registered (§6.6); whether it was run is stated there, and the sixteen
    trajectories reported in §6.1–§6.4 are unaffected either way.
35. **The step budget is therefore unjustified by measurement.** 50 steps and 200 steps give
    indistinguishable results, so whichever is used for the submitted designs is chosen on cost,
    not on evidence that it is sufficient.
36. **The epitope restriction was load-bearing for the pH objective, not an obstacle to it.**
    Freeing the footprint moved the median nearest-histidine distance from 22.8 Å to **34.8 Å**,
    the worst of any condition, with one run at 54.5 Å. The anchor is one pinned point on a
    456-residue surface, and the epitope term was the only thing holding the binder near it. Any
    future free-footprint attempt needs a replacement for that term, not simply its removal.
37. **A tight score distribution was obtained and is not good news.** The four free-footprint runs
    span 0.0045 against 0.085 for the pinned condition, so seed and length stopped mattering. The
    interpretation offered — a single generic surface-contact attractor, reproducible and mediocre
    — is an interpretation, not a measurement.
38. **Design-time scores are optimistic by an amount that is not a constant.** Measured
    design-to-re-predicted translation across four paired runs: +17%, −0%, −49%, −60%, with the
    two largest drops on the two highest design-time scores. No design-time number in this work is
    comparable to a re-predicted threshold, and no correction factor exists.

39. **Optimizer space was never varied, so "the search is exhausted" is narrower than it sounds.**
    All sixteen trajectories held the gradient estimator and the initialization constant: one sample
    per gradient evaluation, momentum 0.9 throughout the soft phase, and a binder starting undocked at
    44–93 Å. Weights, step budget and footprint were varied; sampling, momentum and initialization
    were not. An external method differing on exactly those axes reports binders on this target, and
    the budget to test that axis here was deliberately not spent (§6.5).

### The submission

40. **The CSV ordering is not a calibrated prediction.** No instrument available to us ranks by
    predicted affinity, so the order is a documented nomination priority and is labelled
    provisional.
41. **The matched pairs are only interpretable together**, and Track 3 does not guarantee both
    members are synthesised.
42. **Novelty was assessed with the platform's checker**, which changed during the competition and
    whose implementation is not published, so we cannot reproduce its verdicts locally.
43. **The pH mechanism and developability are in tension.** Raising free pKa raises pI, and
    elevated pI is associated with faster clearance. This challenge does not measure that side.
44. **The inverted-objective pilot, if any design comes from it, rests on an unvalidated
    inversion** — the published validation concerns acidic-pH binding to a different target.

### This document

45. **Sections 8–10 and 12 describe methods whose results did not exist when written**, and are
    marked as such rather than filled with projections.
46. **Several numbers in this work were corrected after first being recorded.** The computed pKa
    shift was described as "measured"; a 10× ratio was described as a requirement when the
    published objective names none; a geometry check measured the wrong protomer of a homotrimer
    for a full day. Each is recorded in §11 with its consequence.
47. **The errors found in §11 are the ones we found.** Two of them were invisible to the tests
    written to catch them, and one was inside a fix that had already been reported as complete.
    The rate at which this work discovers its own faults is not evidence that it has run out of
    them.

48. **The design dependency was unpinned for the whole campaign.** Mosaic was installed from an
    unrevisioned git URL until 2026-10-07, so a fresh clone resolved it to that day's HEAD. The
    exposure was then measured and found empty — the only commits in the run window changed
    `CITATION.cff` and `README.md`, +21/−0, no Python — and the install is now pinned to
    `b94b9d4eb9907a700a6d78ed2d29d3704c5df46c`. **The absence of a confound here was established after
    the fact, not guaranteed by the setup**, and nothing in the test suite would have caught it.

49. **The generator was never validated on this target and the calibration that would have shown it
    existed throughout.** §4.5 measures our design oracle at AUC 0.781 against 0.908 and 0.951 for
    two alternatives, on 150 assayed designs released months before this challenge. We built four
    ranking instruments and validated none of our generation path.

50. **No design in this submission clears the monomer-foldability floor in the organisers' own
    protocol.** 0 of 35 at 0.70; **our best is 68.8 against their lowest *assayed* value of 68.9** — the two
    figures were printed the wrong way round in an earlier draft of this limitation. We report
    this rather than relaxing the floor, and it is the strongest single reason to expect these designs
    not to bind.

51. **The AUC that chose our ranker was measured on a configuration we never ran.** The
    calibration's `ef2full` arm uses target-chain MSAs; our run is single-sequence on both chains,
    with `grep -ci msa` returning 0 in the wrapper and its patch. So 0.901 does not transfer, and
    §4.5(b) no longer rests on it. What we have instead is measured on our own configuration: all
    35 designs scored, each against a shuffled null at its own length, with the min and the mean
    ordering the set identically. An earlier version of this limitation said no ESMFold2 prediction
    of our 35 existed; that was true when written and is no longer.

52. **The calibration in §4.5 is not blind.** Medians for several models were inspected before AUC was
    chosen as the statistic. The guard is completeness — all 160 cells are printed, so no subset can
    be selected quietly — and permutation p values are computed only for the top cells. Single cells
    should be read as descriptive.

53. **A length hypothesis was raised and retracted within one evening.** The argument that 76–84 aa is
    too small for this surface class was retracted when the released data showed all twelve measured
    binders at 80 or 83 aa against a 60–115 range tested. It is recorded because the retraction was
    driven by measured outcomes overruling an a priori geometric argument, which is the direction this
    project has had to correct in most often.

54. **Every interface number in this work is a single-seed draw, and the scorer's seed noise is
    large relative to our margins.** Their data uses five seeds per complex; we used one, for cost.
    The within-design seed SD on `pae_interface_mean` is 0.607, so a margin — a difference of two
    single draws — carries ~0.858 of noise. Our median design's margin is 0.3 SDs of that, and our
    best is 3.5 SDs where pure noise reaches 2.56 at the 95th percentile for a best-of-35. The
    figure is also measured on their MSA-fed arm while ours is single-sequence, so it is a lower
    bound on our noise rather than an estimate of it. §4.6.
