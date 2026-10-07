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

**Status.** A controlled single-variable probe is in progress at the time of writing: identical
configuration and step count, pH weight reduced from 2.0 to 0.5, four trajectories against the
four-trajectory baseline above. Since only the ratio of weights affects a weighted sum, reducing
the pH weight is equivalent to raising the binding weights, and it was the change reachable
without modifying code mid-experiment. The success criterion was fixed before the probe was
launched: **a median iptm_repred of 0.20 or above** indicates the ratio is a usable lever.

*Results from this method do not yet exist. Nothing in this section is a claim about output.*

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

*Method fixed; results pending generation.*

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

*Design fixed; results pending generation.*

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

---

## 13. Reproducibility

- All analysis code is published. Every instrument carries a `--selftest` that runs without a GPU,
  a Modal account, or network access, and the launchers refuse to start if it fails.
- Self-tests include **mutation tests**: each asserts that the *discarded* implementation fails, so
  a test cannot pass if the fix it guards is reverted. The protomer fix and the top-2 reduction are
  both pinned this way, and both were verified by reverting the fix and confirming the failure.
- The `ipSAE_min` pre-registration is a commit timestamped **2026-10-06 16:24:41**, preceding the
  scoring run.
- Pose caches are gitignored and run to tens of gigabytes; they are available on request.
- Numbering is machine-checked against residue identity rather than asserted, on every run.

---

## 14. Limitations

1. **The pKa shifts are computed, not measured.** PROPKA on our own poses. The 2.8–3.9× two-site
   estimate inherits that status.
2. **No measured de novo binder against TNF-α exists** to anchor a confidence threshold, so the
   iptm gate in §8 is a necessary condition only.
3. **n = 4 on the affinity ladder.** Sufficient to retire a prespecified use, insufficient to
   characterise the method.
4. **Target stability at assay concentrations is unconfirmed** (§9). Handled by the matched pairs
   rather than assumed away; the surviving exposure is that both members of a pair must be
   synthesised.
5. **The assay's detection limits for this challenge are unconfirmed.** "No detectable binding"
   depends on them, and the limits used elsewhere in this project were carried over from a
   different target.
6. **The pH mechanism and developability are in tension** (§8), and this challenge does not measure
   the developability side.
7. **The species leg uses two protomers rather than three** because three exhaust GPU memory. A
   control establishes two recover 86% of the signal, with ~20% seed attrition to non-assembly.
8. **Avidity is uncharacterised.** A trimeric analyte over immobilised binders can engage more than
   one site, which raises apparent affinity and compresses the pH ratio being ranked.
