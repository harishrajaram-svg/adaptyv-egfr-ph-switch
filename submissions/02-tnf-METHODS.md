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
| 68 | H | **absent** | dropped — **deleted in mouse** |
| 80 | V | I | dropped |

**Position 68 is why conservation could not be obtained by sidechain selection alone.** Human
carries a histidine there; mouse has no corresponding residue. A deletion moves the local
backbone, and no choice of which sidechains to contact repairs a backbone difference. Epitope Cα
RMSD mouse-vs-human over the retained positions is **1.10 Å** — similar surfaces on average, with
at least one local difference that an average conceals.

**Mouse epitope indices are not the human ones.** The deletion shifts the register by **−3 at the
anchor and −4 further along**: human `16,27,28,70,72,81,82,85,86` maps to mouse
`13,24,25,66,68,77,78,81,82`. Mouse 68, 81 and 82 collide numerically with human 68, 81 and 82
while denoting different residues.

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

**All four miss the bar. Condition 2 remains the best, and it was the second thing tried.**

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
intended allocation. **Twenty-one distinct sequences exist**, one per trajectory across every
condition run, at two lengths. That is more than the twenty slots. **None of them is a success by
this work's own instruments:**

| | count of 21 |
|---|---|
| sequences, all distinct | 21 |
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

**Diversity is the binding constraint on this material, not count.** All 21 belong to one footprint
family — the pinned nine positions — and span two lengths. §14 records single-epitope dependence as
the largest correlated risk in the submission, and 21 sequences from one family do not reduce it.
The second family is the subject of the test in §6.3's open causes.

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

### The design method

19. **The species leg uses two protomers rather than three.** *Deliberate, forced:* three exhausts
    a 48 GB GPU. A control establishes two recover 86% of the three-protomer signal, with ~20% seed
    attrition to non-assembly.
20. **Interface confidence has not exceeded 0.200 in this work**, against 0.45 as the threshold at
    which a geometry verdict is reported at all. No design has reached a value at which placement
    is treated as verified.
21. **Weight tuning appears to have a ceiling.** Three weight conditions spanning a sixfold range
    in the pH-to-binding ratio give medians of 0.136, 0.179 and 0.156. The best is the middle
    condition, which indicates an optimum rather than a direction. Quadrupling the step budget at
    that optimum gives 0.163 — a fourth condition inside the same band.
22. **The weights interact.** Tripling the binding weights pushed histidine content to 9.2%
    against an 8% cap and moved the closest histidine placement from 2.62 Å to 9.69 Å. The loss
    terms are not independent knobs.
23. **The composition cap is a soft hinge, not a constraint.** It has been exceeded. A design can
    ship above its nominal cap.
24. **The pH term is numerically flat beyond ~30 Å**, so it supplies no gradient at the distances
    trajectories begin from. The mechanism depends on the contact terms succeeding first.
25. **The pH term is measured CA-to-cation**, because design-time features give the binder no
    sidechains. The all-atom distance is only available on re-prediction.
26. **Seed variance is large relative to the effects being measured.** Four trajectories at
    identical settings spanned 0.111–0.156, and across the eight runs at the best weight setting
    seed and length alone span 0.115–0.200 — a range about forty times the measured effect of
    quadrupling the step budget. Differences between conditions of that order are not resolvable
    at n = 4, and every condition difference reported in §6 is of that order.
27. **The epitope is restricted to nine positions conserved to mouse.** *Deliberate:* four
    positions that differ were dropped, including one deleted in mouse, which forgoes whatever
    affinity those contacts offered on human.
28. **A single epitope carries the whole submission.** All candidates target one site, so a wrong
    epitope choice fails the set at once rather than independently. This is the largest correlated
    risk in the submission.
29. **Cysteine is excluded from the alphabet**, so no design can use a disulfide for stability.
30. **The interface objective did not improve in any trajectory run.** Twelve trajectories, three
    weight conditions, two step budgets: the design-time interface term is flat from the first
    block of steps to the last (§6.3). The generation method as configured is not optimising the
    quantity it is written to optimise, and **no result in this document should be read as
    evidence that gradient descent on this loss produces interfaces.**
31. **The cause of that flatness is unidentified.** Three candidates remain open — the pinned
    epitope admitting no gradient path, sampling noise in the structure predictor swamping the
    gradient, or an interaction among the eleven weighted terms. Only the first is scheduled to be
    tested (§10's free-footprint family), and it is being run for an independent reason.
32. **The step budget is therefore unjustified by measurement.** 50 steps and 200 steps give
    indistinguishable results, so whichever is used for the submitted designs is chosen on cost,
    not on evidence that it is sufficient.
33. **Design-time scores are optimistic by an amount that is not a constant.** Measured
    design-to-re-predicted translation across four paired runs: +17%, −0%, −49%, −60%, with the
    two largest drops on the two highest design-time scores. No design-time number in this work is
    comparable to a re-predicted threshold, and no correction factor exists.

### The submission

34. **The CSV ordering is not a calibrated prediction.** No instrument available to us ranks by
    predicted affinity, so the order is a documented nomination priority and is labelled
    provisional.
35. **The matched pairs are only interpretable together**, and Track 3 does not guarantee both
    members are synthesised.
36. **Novelty was assessed with the platform's checker**, which changed during the competition and
    whose implementation is not published, so we cannot reproduce its verdicts locally.
37. **The pH mechanism and developability are in tension.** Raising free pKa raises pI, and
    elevated pI is associated with faster clearance. This challenge does not measure that side.
38. **The inverted-objective pilot, if any design comes from it, rests on an unvalidated
    inversion** — the published validation concerns acidic-pH binding to a different target.

### This document

39. **Sections 8–10 and 12 describe methods whose results did not exist when written**, and are
    marked as such rather than filled with projections.
40. **Several numbers in this work were corrected after first being recorded.** The computed pKa
    shift was described as "measured"; a 10× ratio was described as a requirement when the
    published objective names none; a geometry check measured the wrong protomer of a homotrimer
    for a full day. Each is recorded in §11 with its consequence.
41. **The errors found in §11 are the ones we found.** Two of them were invisible to the tests
    written to catch them, and one was inside a fix that had already been reported as complete.
    The rate at which this work discovers its own faults is not evidence that it has run out of
    them.
