# HANDOFF — Anthropic × Adaptyv 2026, Track 3 (all challenges)
# Updated 2026-10-08 6:35 PM EDT. Covers BOTH challenges; problem 1 is closed and archived below.

**Read order on a cold start:** this file → `projects/anthropic-adaptyv-2026/ROADMAP.md` (the live
plan) → `challenges/02-tnf-alpha.md` §1–§48. Then nothing else unless a link below says so.

---

## 🟡 AS OF 2026-10-08 6:35 PM: the Genie arm ran, and the instrument that was missing all along says it failed. Nothing spending. The design search is CLOSED.

`modal app list` is empty. **Two pre-registered arms ran on 10-07 and both returned nulls**, and both
carried hard stops written before the run. Those stops are honoured: **no further design arms.**

| arm | what it varied | pre-registered bar | result |
|---|---|---|---|
| **A** (§45) | the gradient estimator — 4-sample averaging, momentum 0.9 → 0.0 | ≥ +0.12 acts, < +0.02 null | **+0.0078 ± 0.0124 → NULL** |
| **Region I** (§47) | the **epitope** — a second site 15.42 Å away | same bar | **−0.0092 ± 0.0126 → NULL** |

**Twenty-one trajectories now, six conditions, two structurally unrelated epitopes.** The interface
term improved in none of them. Best `iptm_repred` anywhere is **0.2468** (arm A) against the **0.45**
gate we set before any number existed. 0 of 35 designs clear it (the 10-07 export
of `design_inventory.py --export` regenerated 25 -> 35).

**What survives as the explanation:** term interaction across 20 weighted terms on three legs.
Condition 3 is direct evidence — tripling the interface weights made results *worse*. Testing it
means stripping terms, and the terms *are* objectives 1 and 2, so the suspect cannot be fixed
without abandoning the brief.

**Spend:** ~$36 before today, ~$9 arm A, ~$9 Region I, ~$3 lost to two launch mistakes of mine.
**The $49 wave money is unspent and now has nothing to buy.**

### Arm A's endpoints were the project's best, and that is a trap, not a lead
Median 0.2159, max 0.2468 — above every other condition. **On a null trajectory.** §6.3 is why
endpoints are secondary: an endpoint cannot tell a badly-tuned method from one that is not
optimising. Do not reopen arm A on the strength of its endpoints.

---

## 📌 THE SUBMISSION IS PLAN C, AND THREE THINGS CHANGED TODAY THAT MAKE IT EASIER

**1. 🟢 THE NOVELTY GATE STOPS BEING A SUBMISSION FILTER — portal switches 2026-10-08.**
Simon (Adaptyv), 2026-10-07 6:46 PM:
> "switching to a **sequence based presubmission screening, which will flag only submitted sequences
> with resemblance to known TNFa binders** … **Novelty will still be computed but only after
> submission** together with other metrics on Proteinbase."

Offline tool: **`github.com/adaptyvbio/binder-prescreen`** — *"so that you can check your submissions
(or have Claude check them) before you submit."* Sequence-based, no structure needed.
**Consequence:** problem 1's novelty roulette (opposite verdicts on byte-identical sequences, cost us
a backbone family and the only two-site mechanism) **cannot recur**. Submit the best 20 by our own
ranking. ⚠️ Novelty has not vanished — it moves post-submission and **feeds selection**, so keep
reporting our levels honestly. **`binder-prescreen` was run on all 35 on 10-07; re-run it on
the final 20 before uploading.**

**2. Fewer than 20 designs is explicitly fine.** Tudor, 2026-10-02:
> "It's ok if you want to submit less than 20 … **we will not pick just 1-2** (so a bit more to give
> some statistical power)."
So padding to 20 buys nothing. Relevant because **§9 already says a matched histidine-removal
control is uninformative when the parent does not bind** — 0 of 35 bind, so the 4 control slots are
better spent or left empty.

**3. Upload early, then "designate".** Théo: *"There is an option on Proteinbase to 'designate' the
submission that will be taken into account when you've made multiple ones."* Multiple submissions are
retained. An early upload is free insurance; it does not burn the 24 h cooldown.

**Deadline, reconciled — no runway was lost:** **Sun Oct 11 23:59 AoE = Mon Oct 12 7:59 AM EDT.**
AoE is UTC−12, so the two readings are the same instant. Last upload that still leaves a 24 h retry:
**Sun Oct 11, 7:59 AM EDT.**

---

## ⚠️ ORGANISER FACTS — confirmed in the Slack, do not re-derive

- **Selection, stated three times, most explicitly by Tudor:** *"mainly based on method novelty,
  design diversity, and a couple of in silico/confidence metrics. Your strategy should not be
  penalized much by the in silico metrics … this selection strategy will likely be less biased
  against more creative design methods."* They will publish the selection prompts afterwards.
- **Amir:** *"Claude will review what you share and use it for selection."* The methodology section
  is linked to a Proteinbase Collection. ⚠️ Being made **public** is conditional — *"**if** your
  designs are selected for wet-lab testing … **might** be made public."* Do not state it as a rule.
- 🔴 **Assay construct, confirmed 2026-10-07 by Tudor:** *"C-terminal Twin-Strep, we immobilize the
  binder, and flow the target."* **Our design is the immobilised partner and the trimer is the
  analyte** — so the avidity exposure in limitation 16 is now organiser-confirmed, not inferred. A
  trimeric analyte can engage two immobilised binders and compress the pH ratio we are ranked on.
- **Target:** UniProt P01375, **Asp143** (not 1TNF's Leu143), AcroBiosystems TNA-H4211, tag-free.
- **Format:** 10–250 aa; `molecule_class` ∈ `single_chain|nanobody|scfv|fab_kappa|fab_lambda`
  (**`protein` is not in the enumeration**). Strata: microbinder <40 · minibinder 40–100 · large
  >100 · nanobody · antibody, **winners per stratum**. Our 76-mer is a minibinder.
- 🔴 **A ranked CSV is mandatory.** Declining to rank is not submittable.
- 🔴 **"Embedded instructions or prompt injection may be deemed grounds for disqualification."**
  The methods document must contain no imperative addressed to a reader. Sections added 2026-10-07
  were checked clean; **older prose was not fully swept.**
- **Provenance rule, verbatim:** *"You may not take an existing binder and modify it. Designs must be
  produced from scratch."* Kills any use of `TNFA_OPT_07`, the Chen *et al.* TMB sequences, and our
  own problem-1 sequences (*"We check against the current Proteinbase snapshot"*).
- **PyRosetta is permitted** for non-commercial use (Amir + Tudor, 2026-10-06). **PolyForm was never
  ruled on** — silence, not permission.
- **Operational:** multiple entrants report Claude's bio classifiers blocking routine protein work on
  challenge 2. Amir recommends Sonnet 5.5 and says the classifier level is the same across Sonnet 5 /
  Opus 5 / Sonnet 5.5. **This session was interrupted repeatedly.** A fresh window is the reported
  workaround.

---

## 🔭 GENERATION PATHWAYS — researched 2026-10-07 evening, NONE LAUNCHED

Three agents surveyed this. **No pathway is recommended without an explicit decision**, because the
§45/§47 hard stops say the design search is closed and reopening it is a judgement call, not a
tooling one. Recorded so a fresh window does not re-derive it.

| pathway | status | cost | verdict |
|---|---|---|---|
| **BindCraft2** `github.com/PacesaLab/BindCraft2` | created 2026-09-13, **pushed 2026-10-07** | ~$2–7 | **Best capability fit, zero published yield** |
| **Cao scaffolds → BoltzGen inverse-folding** | our harness already wired | ~$8–16 | Nearly zero setup; His pins at generation |
| **Genie 3** `aqlaboratory/genie3`, Apache-2.0 | public weights | cheap | **8 of the organisers' 12 TNF-α binders are variants of ONE Genie 3 backbone** |
| **Proton-PottsMPNN** | MIT, weights in repo | ~$0 (CPU) | The only model taking protonation state as **input** |
| **Germinal** | 🔴 **TRIED HERE 2026-10-02 AND FAILED** | ~$80 for 20 | **Out** |
| BoltzGen `nanobody-anything` | — | — | Externally measured **zero on TNF-α, both arms** |
| RFantibody · RIFdock · Chai-2 · Latent-X · IgDesign | — | — | No wrapper / 100k CPU-h / not available / needs a parent |

**Germinal is out on measured grounds, not theory.** `runs/egfr-d3-scfv/`: four attempts on
2026-10-02, **zero designs**, abandoned in 52 min with a `STOPPED_germinal_boxed_out` tombstone.
Failures: `No module named 'yaml'` · return code 1 · `No module named 'pkg_resources'` ·
`ValueError: Unrecognized amino acid token: A` inside vendored `colabdesign/iglm/model.py`. Three
never reached the design stage. Also ~$80 for 20 designs (over budget), hard-requires PyRosetta, and
**§26 measured that our ranking instrument substantially fails on nanobody format** while a ranked
CSV is mandatory.

**BindCraft2, what is actually verified:** ships `scaffolds/{VHH,scFv,Fab,ARP}.cif`;
`examples/pdl1_vhh.json` does fixed-scaffold CDR-only design with variable CDR lengths,
`min_scaffold_sequence_retained_final: 1.0` and `max_off_paratope_contact_final: 0.0`;
`pdl1_crossreactive_detarget.json` carries **two orthologs as a native joint objective** plus
negative design; `pdl1_homotrimer.json` exists. License is source-available, free for any purpose
except hosted resale, **no PyRosetta**. 🔴 **Unverified: yield. No paper, no wet-lab numbers.** pH is
not an input. VRAM ≈ 2.0 × (3.4 GB + 38 kB·N²) → full trimer + 100 aa ≈ 30 GB → **one worker on an
L40S; truncate to a protomer pair or an epitope shell.**

🔴 **RULES TRAP:** BindCraft and BoltzGen both ship **adalimumab** and **golimumab** Fab scaffolds.
Both are anti-TNF-α drugs. Using either as a framework against TNF-α is *"taking an existing binder
and modifying it."* The `nanobody_scaffolds/` set is clean. **Check every scaffold against the
target.**

**Honest expectations if any of this ever runs.** The Baker pH paper's TNF-α row: **72 designs → 48
retained binding → 3 pH-sensitive**, best 79× weaker and that at **pH 5.4**, not 6.0. At a 20-design
cap that is **~0.8 expected switches**. The best computational pH-**6.0** result anywhere is **>2×**.

---

## 🔴 2026-10-08 EVENING — THE GENIE ARM, AND THE TWO MEASUREMENTS THAT SETTLE IT

Full write-up in METHODS §6.10a–f. Total spend today ~$35.

**The arm works mechanically.** Genie 3 at 39.9 s and $0.023 per backbone, 3.39 GB peak, 4 chains
of a 4-chain ceiling. **98 of 100** seed sequences clear the 0.70 monomer floor against **0 of 35**
for the Mosaic designs. 48 designs co-folded, best margin **+6.287** on five seeds that excluded
the selected one, which cleared the pre-registered 5.64 bar.

**Then two measurements that should have come first.**

**1. The margin is validated, and our designs fall below NON-binders (§6.10e).** 12 measured
binders vs 30 measured non-binders through our own configuration: **AUC 0.794, p = 0.0010**. On
that scale binders sit at **+7.814** and non-binders at **+6.347**, while our 48 sit at **+0.320**
with a best of +6.287. **1 of 48** reaches the weakest real binder; **0 of 48** reach a typical
non-binder. The pre-registered bar of 5.64 was anchored on TNFR2, a natural receptor, when the
relevant population is designed binders — so it was set *below* the non-binder median. It was met
and it was mis-set; both are in the record.

**2. Pose agreement is near zero (§6.10f).** The protocol ranks on ipSAE **and sc_DockQ**, 4:1, at
every promotion step. We never computed sc_DockQ. It gives median **0.025**, **0 of 40** passing
the protocol's 0.23, **Fnat ≈ 0**, ligand RMSD median **31 Å** — the predictor does not put the
binder at the designed epitope — and it is **uncorrelated with our margin** (r = −0.159). So every
margin in this work scored an interface we did not design.

**What the arm did establish, and it is worth keeping:**

| | |
|---|---|
| sequence variance **1.6×** backbone variance | SD 1.462 vs 0.931; which sequence matters more than which backbone |
| best-of-4 by MPNN score is worthless | 4 of 24 alternates beat the pick; the new best **+6.459** came from a backbone ranked 7th of 8 |
| their 8 binders are **ONE** family | 23 Genie3 designs in 10 families; one family of 8 gave 8/8 binders, the other nine gave zero |
| objective 3 is near-chance | mouse median **−1.264**, below its own scramble; human/mouse r = +0.122 |
| novelty | **16–21 of 48 at level 2**, zero at level 4, driven by fold similarity (median TM 0.732) not sequence (median identity 0.000) |
| prescreen | **48/48 pass**, `outbox/prescreen-genie48/` |

**One design is not bad at anything:** `tnfa_corrected_30_s1` — human +5.814, mouse +1.991, pH
linkage 0.149. Not good at anything either.

**Two dramatic results today were bugs the controls caught before they were reported**: a
design/null display collision in the Stage A margins, and a five-residue register shift in
sc_DockQ that the target-superposition control exposed. Both tools now refuse rather than report
when their own control fails.

**What this closes.** Stage C — 1000 backbones, 200 co-folds, ~$47 — is **not justified**. Not by
the pre-registration's ambiguity clause but by a clear negative on a validated instrument. And the
plan as written was structurally wrong anyway: 1000 backbones × 1 sequence each scales the axis
that matters least.

---

## 🔴 2026-10-08 — THE RANKER WAS NEVER TESTED AGAINST A NULL OF THE RIGHT SIZE, AND NOW IT IS

Two folds, ~$0.35, `analysis/02-tnf/SIZEMATCHED-NULL-2f.md`. Apps `ap-snCV8ziLzEWYsls08s9Vh2`,
`ap-ihYnoYl5TTQ2UhM8vJA4Kp`.

The control band folded on 10-07 (`p2_control_band_2026-10-08.tsv`) shuffles the 164-residue TNFR2
ectodomain, so all five controls score at 302,382 cross-chain residue pairs. The designs are 76 or
84 residues — 219,486 / 227,022 pairs. `pae_interface_mean` averages over those pairs, so that band
was **never a yardstick for the designs**. Mine: I recommended that run without checking
size-matching first.

| null | binder | pae_if_min | pae_if_mean | pairs |
|---|---|---|---|---|
| `neg_shuffled_L76_1` | 76 | 0.854 | **14.248** | 219,486 |
| `neg_shuffled_L84_1` | 84 | 0.855 | **14.996** | 227,022 |

| | n | design mean | median | null mean | beat it |
|---|---|---|---|---|---|
| L76 | 21 | 11.218–16.586 | 14.047 | 14.248 | 15/21 |
| L84 | 14 | 13.945–15.434 | 14.648 | 14.996 | 10/14 |

**25 of 35 beat their size-matched null on the mean** (one-sided binomial p = 0.0083), and **the
same 25 beat it on the min** — so at this configuration the two metrics rank the designs
identically, which is the part of item 2d that was answerable without an MSA path.

**Above chance, and almost nothing biologically.** TNFR2 scores 7.713 where its own shuffles score
~18.99 — a margin of 11.28. Best L76 design: 3.03 below its null, **27%** of that margin. Best L84:
1.05, **9%**. Median design: 0.20 / 0.35, **2–3%**. The median design is not distinguishable from a
shuffle of its own length and composition. §48 reached the same verdict from Anthropic's 150
measured designs; this reaches it from controls folded on our own configuration.

⚠️ Caveats, stated in the file: one draw per length (n=1, no interval), and the 27% assumes the
null-relative gap cancels the size offset to first order rather than exactly. **2d is not closed** —
`pae_interface_min` was selected at AUC 0.901 on the MSA-fed arm and our configuration has none.

### The liability screen fired, and it is near-binding on the submission

`bin/express_qc_p2.py` (thresholds declared before any problem-2 design existed) on all 35.
`analysis/02-tnf/p2_liability_screen_2026-10-08.txt` + `p2_liability_2026-10-08.tsv`.

| flag | hits |
|---|---|
| pI **inside** the assayed 6.0–7.4 | **13/35** |
| pI within 0.5 of an assay pH | 5/35 |
| near-neutral net charge at an assay pH | 15/35 |
| **any pH/charge flag** | **20/35** |
| no flag of any kind | 7/35 |
| composition / homopolymer / GRAVY / free Cys | 0/35 |

pI range 4.49–11.21, median 6.43. L76: 10/21 inside the range; L84: 3/14.

The hazard was predictable from the assay alone and is written into the tool's docstring: problem 2
is measured at **two** pH values and the mechanism is installed with **histidines**, whose pKa sits
between them by construction. The better the mechanism is installed, the likelier the binder's pI
lands inside its own measurement window.

**Consequence for nomination: 20 slots, and only 15 of 35 are clean on both pH and charge.**
Liability can no longer be a tiebreak — it has to be a term in the ranking, or five slots go to
designs sitting at their solubility minimum during the measurement that decides them. Problem 1 lost
a VHH to exactly this (pI 6.38, assay 6.5) and the QC in place then had no pI term at all.

### Timing fact worth keeping
The first fold took ~19 minutes, the second ~4. The difference is a **cold model fetch**: the
container downloaded 11 weight files (1:44) then loaded for 117 s before folding. `score-esmfold2.sh`
sizes `MODAL_TIMEOUT` from the **fold** only, so its 5-minute estimate is blind to the fetch. The
timeout is a ceiling and Modal bills actual use, so nothing was lost — but do not read the estimate
as wall-clock on a cold image.

---

## 🔬 THE CALIBRATION WAS TOO GENEROUS, AND THIS IS STILL OWED IN THE WRITE-UP

We have been quoting the organisers' own **12/150 = 8%** on TNF-α as our calibration, in METHODS
§6.5, limitation 19 and the public methodology box. **Their campaign required only binding.** Their
protocol (`reference/anthropic-binder-design-protocol.md`, TNF-α is target #10) states the goals as
*"1) high-affinity binders zero-shot and 2) high overall hit rate"*, contains **no pH requirement at
all**, and makes cross-species *"a secondary objective pursued only without compromising affinity or
hit rate."*

**So 8% is the one-objective ceiling, not our calibration.** We are asked for three things at once,
and the added one has no computational precedent at pH 6.0. All three places need qualifying —
it reframes our negative result as the three-objective version of a task whose one-objective version
runs at 8% on the organisers' own instrumentation.

---

## ✅ CORRECTIONS SHIPPED 2026-10-07 — all guarded, do not reintroduce

| what was wrong | now | guard |
|---|---|---|
| Mosaic installed from an **unpinned** git URL all campaign | pinned `b94b9d4eb99…`, in the wrapper **and** the tracked patch | `bin/check-pins.sh`, mutation-tested |
| METHODS §2: positional 68 *"deleted in mouse"*, with an argument built on it | it is **H→Y**; the single gap is at **mature 71 (Ser)**, not an epitope position | `bin/check_species_map.py`, 6 mutations |
| *"an L40S cannot hold three protomers"* | human trimer = 456 res / **532 tokens, 21.7 s/step, runs**. The **dual-species** graph with a mouse trimer (~976 tokens) is what OOMs | stated in §6.6 |
| methodology box: `+20%` effect size | **`+31%`** — §11 fixed it at 7:35 AM and it was never carried across | — |
| ROADMAP: methodology *"linked to a public Collection"* | conditional — *may* be public *if* selected | — |
| novelty gate = *"3 of 4 checks"* | a **4-level scale, bar at Level ≥ 3** | — |

**`bin/gate_sweep.py` now runs 24 gates** (`check_species_map` and `check_published_counts` added
since the 21 of 10-06). The methodology box no longer states a gate count at all, which is why
`check_published_counts.py` audits every surface that does -- it was extended on 10-08 to cover
gate counts after this file and ROADMAP.md were both found still saying 22.
⚠️ I briefly wrote "15 of 15" into that box from a bad regex — fixed, but re-read it before upload.

---

## 🔧 TOOLING ADDED 2026-10-07 — one command each

```
python3 bin/check_species_map.py                  mouse mapping + gap + Region I, from sequence
python3 analysis/02-tnf/region1_verify.py         Region I by motif, in all three numbering schemes
bash    bin/check-pins.sh                         Mosaic revision pin, wrapper AND patch
python3 analysis/02-tnf/loss_traj.py LOG --block 13   the trajectory read; the ONLY decisive statistic
bash    bin/probe-grad-noise.sh  <tag> <L> <seed>     arm A (GRAD_SAMPLES/MOMENTUM_SOFT env)
bash    bin/probe-region1.sh     <tag> <L> <seed>     the Region I arm
python3 bin/gate_sweep.py                         all 24 gates
python3 bin/express_qc_p2.py <designs.csv>        liability screen (pI at BOTH assay pH values)
python3 bin/build_sizematched_null.py <design> <out.faa>   a null at a design's own length
python3 bin/pae_interface.py --dir <run>          pae_if_min / pae_if_mean / n_pairs
python3 bin/check_published_counts.py             limitation AND gate counts, every surface
```

**Launcher guards, learned the hard way today:** a tag containing whitespace is refused, length and
seed must be digits, and every launcher echoes its resolved config **before** launching. This exists
because `set -- $a` in a **zsh** loop does not word-split, so four runs silently became four
identical `L=76 seed=0` runs. Also: **every TNF launcher must name `--epitope` on both legs** —
`DEFAULT_EPITOPE` is problem 1's **EGFR** site and inheriting it designs against the wrong protein.

---

## 📋 WHAT IS LEFT — all writing, except the upload

1. **Three documents now overstate the ranker**, listed in `analysis/02-tnf/CONFIG-COMPARABILITY-2d.md`:
   METHODS §4.5, limitation 51, and the box's FIELD 2 point 4 each present `pae_interface_min` as
   selected at AUC 0.901, which was measured on the MSA-fed arm. 2f adds the finding that at **our**
   configuration the min and the mean rank the 35 identically, so the claim to make is "either metric,
   same order, margin 2–3% of a real binder's" — not an AUC we did not measure here.
2. **METHODS:** fold in arm A (§6.6) and Region I (new §6.7) with both hard stops recorded as
   honoured; qualify §6.5's 8% per the section above; add the pH-6.0 ceiling with citations
   (Ahn *et al.* bioRxiv 2025.09.29.678932; Schröter 2014 doi:10.4161/19420862.2014.985993); state
   that **no structure predictor represents the pH-6.0 state**; note our **2** His-cation contacts sit
   at the published floor where working designs had **8 and 11**; record Germinal as an attempted,
   failed arm; add the size-matched null band and the liability screen.
3. **Sweep older prose for reader-addressed imperatives** — disqualification risk, not style.
4. **Three stale claims:** `submission-basis-p2.md:48` says `molecule_class: protein` (code is
   already correct); `setup-checklist.md:24` still says avoid PyRosetta; problem 1's METHODS §27 owes
   the PyRosetta correction with today's date.
5. **Nominate ≤20 + ranking rationale.** The screen changed what this is: **only 15 of 35 clear
   both the pH and charge flags**, so liability is a ranking term, not a tiebreak. Needs the item-1
   decision first.
6. **Re-run `binder-prescreen`** on the final 20 (10-07's pass was on all 35).
7. **Upload early, designate, attach methods, finalise.** ⬅ **Harish's**
8. **Send the reviewer their 4 artifacts.** ⬅ **Harish's** — critical path, see the note below.
9. Optional, $0, still unasked under his own name: **does a large K_D shift with binding at both pH
   values qualify, or must binding be undetectable at 6.0?** It decides whether any single-site
   design counts. Two other teams asked; never answered.
10. Open and undecided: **a Modal cap.** October billed **$590.71 with no credit balance**, while
    project documents still record $168.44 from 10-04 — stale by ~$420. Yesterday and today added
    ~$2. Nothing is running, but no cap exists.

**⚠️ ROUTE SCIENTIFIC JUDGEMENT TO THE REVIEWER, NOT TO CLAUDE.** Across 2026-10-07/08 more than six of my
responses were stopped mid-generation by a safety classifier, every one of them on explanatory prose
about this project, while every tool call, file write, commit and GPU run completed. Execution here
is reliable; reasoning delivered in prose is not. That makes item 8 critical path rather than
courtesy. Artifacts and tables have gone through consistently — prefer them.

## 🧭 THE PLAN LIVES IN `ROADMAP.md` IN THE OTHER REPO

`~/code/context-directory/projects/anthropic-adaptyv-2026/ROADMAP.md` (updated 2026-10-08 9:20 AM).
Everything below this line is problem 1, closed and archived.

---

# ══════════ PROBLEM 1 — CLOSED, ARCHIVE BELOW ══════════

**✅ SUBMITTED 2026-10-07. 16 designs, all 16 cleared novelty. Nothing owed.**
Methodology box pasted, METHODS attached, repo public and pushed, gates green.

**The one carry-forward lesson:** the resubmission swapped the set's *composition*, not its size.
Both designs the 10-05 check rejected were accepted on 10-06, and two it had passed at 3/4 were
rejected — **on byte-identical sequences.** Their novelty check is not stable, so every novelty
level we state is one run of their pipeline on one day. METHODS limitation 37; lessons 10 and 11
in `projects/anthropic-adaptyv-2026/lessons-problem-1.md`.

Everything below this line is the problem-1 record, kept because the lessons carry forward. **No
action in it is outstanding.**

---

## 1. STATE — everything below is verified, not remembered

**Submission: `submissions/01-egfr.csv`, 16 designs** (the resubmission set, as accepted). Track 3
allows 20 — see METHODS §11. *(This line read "18 designs" until 2026-10-07; the paragraph two
below has recorded the drop to 16 since 10-06, so the file contradicted itself for a day.)*
Ranked on the two-partner histidine-only pH product, not the superseded target-only basis.
`bin/check_discards.py` PASSES, exit 0. **`bin/gate_sweep.py`: 14 gates, all green at 16
designs** (d63c61f).

🔴 **Both counts in this paragraph were stale and are corrected 2026-10-06.** It read
"18 designs" and "12 gates: 11 green, `novelty_coverage` RED". What happened in between, from
the commit record rather than from this file: the 18-row set **was** pushed to Proteinbase on
10/5 ≈ 12:35 PM, their novelty filter scored 16 at 3/4 and **two at 2/4**, and it **blocked the
submission** until those two were removed — `bc_s360518_mpnn9_A22D` (Level 3 on our gate at
qTM 0.7924, clearing the cliff by 0.0076; limitation 33 named it in advance as the sharpest
eligibility exposure in the set) and `h370_020_vhh`. Removed as an INELIGIBLE map in
`bin/emit_submission_csv.py`, so a re-emit cannot reinstate them. `novelty_coverage` went green
on 10/5 when foldseek was reinstalled and the gate re-run — all four unlevelled chains cleared
Level 3. **`rimA02_d3_rimA_14_vhh` PASSED at 3/4**, which retires the ANARCI worry recorded
elsewhere: the one surviving antibody-format row is known-good on their gate.

✅ **ANSWERED 2026-10-06 11:50 AM, from the portal. IT WAS UPLOADED AND IT IS DESIGNATED.**
~~UNRECORDED: was the corrected 16-row set ever successfully uploaded?~~

Read first-hand while logged in as `harishrajaram-svg`:
- The profile carries a collection **"[Anthropic × Adaptyv] Submission 1", 16 proteins, marked
  DESIGNATED**, created 10/5. Its 16 AA lengths match `submissions/01-egfr.csv` exactly.
- `/rounds/1/submit` returns *"You've already submitted today … again after 6 Oct, 16:49 UTC."*
  A blocked upload does not start that cooldown, so the accepted submission went in at
  **5 Oct 16:49 UTC = 12:49 PM EDT**, fourteen minutes after `d63c61f` removed the two designs
  their novelty filter rejected.

**DESIGNATED means it is the one that counts, so no nomination step is outstanding.** Nothing is
owed on problem 1. The only thing still live is an opportunity, not a risk: challenge 1 stayed
open until Wed 10/7 7:59 AM and one further submission was available from 12:49 PM on 10/6 —
submissions are retained with the designated one counting, so an upload could not have lost what
was already banked.

🔴 **The 10-row table that stood here is replaced, 2026-10-06 — it was wrong twice over.**
It listed 10 of 16 shipped designs, it still carried `h370_020_vhh`, which Proteinbase rejected,
and its numbers came from the **`ph_ratio_target_only_SUPERSEDED`** column — contradicting the
sentence directly above it, which declares the ranking basis to be the two-partner his-only
product. Below is all 16, generated from `submissions/01-egfr.csv` on that declared basis
(`ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE`), never typed. The basis change reorders the set
heavily — the four S15D designs drop from ranks 1–4 to 11–14 — so the CSV's own
`ph_rank_range_across_bases` is carried in the last column. **If the intended ranking key is a
different column, say so; this table follows the sentence above it.**

     1  5.656  rimA01_r15_L133E                              hu 0.5978  mo 0.4338  150aa  rank-range 1-6   
     2  5.546  c5_cf_short__boltzgen_egfr_cropfree_short_48  hu 0.2415  mo 0.1812   70aa  rank-range 2-10  
     3  4.838  rimA02_d3_rimA_14_vhh                         hu 0.2186  mo 0.4468  129aa  rank-range 3-13    nanobody
     4  4.812  c5_cr_crop_patch__boltzgen_egfr_crop_patch_05 hu 0.1320  mo 0.2041   66aa  rank-range 4-11  
     5  4.256  rimA01_r15_boltzgen_egfr_d3_rimA_20           hu 0.5938  mo 0.5668  150aa  rank-range 3-14  
     6  3.545  ss_bc_s831683_mpnn6_S15D_S62H_routeA          hu 0.7654  mo 0.7443   65aa  rank-range 4-11  
     7  3.526  d2c_mpnn13_S88D_serasp                        hu 0.6031  mo 0.5283  147aa  rank-range 5-9   
     8  3.180  cons_gap_h370_only__boltzgen_egfr_h370_018    hu 0.4566  mo 0.2150   90aa  rank-range 4-15  
     9  2.914  bcr_d3acid3_l60_s647537_mpnn3                 hu 0.5736  mo 0.1681   60aa  rank-range 10-16 
    10  2.747  bcr_d3acid3_l60_s647537_mpnn11                hu 0.4429  mo 0.2342   60aa  rank-range 11-17 
    11  1.835  bc_s831683_mpnn6_S15D                         hu 0.7803  mo 0.7507   65aa  rank-range 2-13  
    12  1.774  bc_s831683_mpnn19_S15D                        hu 0.8077  mo 0.7859   65aa  rank-range 3-14  
    13  1.062  bc_s831683_mpnn9_S15D                         hu 0.8040  mo 0.8037   65aa  rank-range 7-16  
    14  1.023  bc_s831683_mpnn8_S15D                         hu 0.7760  mo 0.7637   65aa  rank-range 8-17  
    15  0.737  bc_d3acid_l65_s831683_mpnn11                  hu 0.7963  mo 0.7838   65aa  rank-range 5-17  
    16  0.627  bc_s831683_mpnn9_WT                           hu 0.7857  mo 0.7840   65aa  rank-range 12-18 

**Methods document: `submissions/01-egfr-METHODS.md`, 12 sections, DONE.** Nothing is owed to
it. Leads with §4.4.

**Git remote: EXISTS, public, pushed.**
<https://github.com/harishrajaram-svg/anthropic-adaptyv-2026> — 200 on the README, the methods
doc, the CSV. *(Renamed from `adaptyv-egfr-ph-switch` on 2026-10-07; the old URL 301-redirects.)* The README is a submission front door, not the old setup log (that moved to
`docs/setup-notes-2026-09-18.md`).

Every pH number is the MEDIAN over every ESMFold2 refold pose of that exact binder sequence,
pooled across runs, joined to affinity BY SEQUENCE in `analysis/01-egfr/master_rank.json`.
**That file is the only source for both columns. Do not re-derive either from a run name.**

---

## 2. [2026-10-04] The result that changed that day — control recovery

`expctrl`/`ctrl2` had written raw PAE matrices and **no `*_10_10.txt`**, so the ipSAE step had
never run on them and the only measured-outcome test in the project was unavailable. Backfilled
`bin/ipsae_min.py --dir` over all 8 shards: **270 poses, 0 failures.** `truenull` was scored
after (200 poses; it is the superseded 1-shuffle×5-seed design, kept for n only).

**The panel is 10 negatives and ONE positive. Do not quote an AUC off it** — with one positive
every summary statistic is just "where does that one molecule rank", and an earlier draft of
METHODS §4.4 led with AUC 0.80 / 1.00 before this was corrected. Report ranks.

**Human leg: 8 of 10 RIGHT-CENSORED molecules rank below the measured binder, 2 rank ABOVE it.**
They are molecules with no reported KD, not molecules measured not to bind (the reviewer, 2026-10-04;
METHODS §4.1) — this line said "MEASURED non-binders" after that correction had landed.
The top-scoring molecule in the entire measured panel is one with no reported KD (right-censored, not a measured zero) —
`EXPNEG_gitter-yolo10` at 0.5893 against human EGF's 0.3549 — and it reads a **5.27× pH ratio,
rank 18 of the 246 rankable at n ≥ 5**, the top 7.3%. "Rank 8 of 2,009" appeared here and in README, and "8th of 132" replaced it without reproducing either; METHODS §7
records it as a ~15× overstatement, because 1,877 of that denominator were never rankable.

**Requiring both species: the one quantified binder outranks all ten no-KD molecules.** Both
molecules that beat it on human are at **exactly 0.0000 on mouse, 5 of 5 dead seeds** — no
interface at all, so the filter does not rest on a margin that a bigger panel could erode.
That mechanism is the defensible part; the ordering statistic is not. Cross-reactivity came
from the brief and is the only specificity filter here that measured data supports at all.

**The matched null is degenerate:** 12/12 shuffles at exactly 0.0000, 5/5 dead seeds. So
`affinity_above_null` means "above zero" and carries almost nothing. **CDR decoys on the real
5 nM VHH framework also read 0.0000 — same as the real VHH** — so the nanobody blind spot is
total, not the partial 0/2 of the old panel.

Artifacts: `analysis/01-egfr/control_recovery.{json,tsv}`. Written up as METHODS §4.4 + §4.5.

---

## 3. [2026-10-04] Decisions made that day — do not relitigate

- **Submission cut 20 → 10** (Harish, 3:30 PM). Ranks 11–20 did not stand on a measurement:
  8 read below 1.0× on the 0.702× analytic acid-weakening extreme. (Corrected 2026-10-05:
  that is NOT a "no switch" reading — no linkage gives 1.0 — and a pile-up on an analytic
  extreme indicates protonation-model saturation, so those 8 are uninterpretable on this
  gate, which is still a reason not to submit them but a different reason.) `bg04_r03`
  was 1.94× with 0.0000/0.0000 and an unpaired cysteine. `LIMIT = 10` in
  `bin/emit_submission_csv.py`, reasoning recorded in the docstring.
- **Six of ten slots on one backbone** (`d3acid_l65_s831683`) — four S15D sequences, mpnn11, and
  the matched WT. Deliberate: diversity traded for wet-lab replication of the only causal result.
- **Tier-2 sort key: MOOT and left unchanged.** At LIMIT = 10 no tier-2 design ships. Not edited
  blind at submission time.
- **Repo public** (Harish, 3:30 PM), full history including the error-correction commits.
- Standing: one submission made LATE; affinity gate retired (`MIN_AFFINITY = None`); mouse
  dropped as a tier-2 sort key; Mosaic gated for problems 2–5; S88D replaces its parent.

**New fact worth carrying:** inverting the linkage equation, the four S15D ratios imply
pKa_bound of 8.89 / 9.11 / 8.98 / 8.98 — the whole 5.40–5.46× spread is **0.22 pKa units wide**
and the four are NOT distinguishable from each other. Do not claim rank 1 beats rank 4.

---

## 4. WHAT WAS LEFT — one item survives, and it is still open

1. ~~**UPLOAD, and upload early enough to resubmit.**~~ ✅ **DONE** — uploaded 10-05 12:49 PM EDT,
   designated, 16 designs, all cleared novelty on the 10-06 re-check. The VHH worry recorded here
   resolved: `rimA02_d3_rimA_14_vhh` passed at 3/4, and `h370_020_vhh` was one of the two the
   10-05 check rejected and is not in the shipped set.
2. 🔴 **STILL OPEN — send the reviewer their artifacts.** All of their asks are done and **still
   unsent**: `outbox/ipsae-fixtures/` (8 cases, `run_fixtures.py --check`, VHH-zero trace),
   `outbox/CONTROL-TABLE.md`, `outbox/PREREGISTRATION.md`, and the construct audit. §4.4 answers
   their control-recovery ask and they should see it. **This is ROADMAP item 19 and it is Harish's
   to send** — the Gmail connector drafts only, it never sends.
3. ~~Optional: better designs hiding in the `sd`/`sd2` arms of `master_rank.json`.~~ **Moot** —
   problem 1 is submitted and designated, and iterating on a submitted design is disallowed.

---

## 5. TOOLING

    bin/ph_gate_refolds.py      pH gate on ESMFold2 refolds -- run on ALL 131 run dirs
    bin/ph_pool_by_sequence.py  ONE ratio per MOLECULE, pooled by binder sequence
    bin/master_rank.py          pH + affinity joined BY SEQUENCE -- the single ranking table
    bin/instrument_v2.py        ipSAE_min, v1 max / v2 median over seeds   (--self-test)
    bin/ipsae_min.py --dir D    *** writes the _10_10.txt every other reader parses ***
    bin/novelty_gate.py         FOLDSEEK=tools/foldseek/bin/foldseek ... --db tools/fsdb/pdb
    bin/antibody_novelty.py     Adaptyv's ANTIBODY branch (CDRH3)          (--selftest)
    bin/express_qc.py           expression liabilities (OVERWRITES express_qc.tsv)
    bin/emit_submission_csv.py  builds the CSV                             (--selftest)
    bin/check_discards.py       hard gate: nothing unmeasured may outrank a submitted design
    bin/bindcraft-sync.sh       pull BindCraft output off the Modal volume. RUN AFTER EVERY JOB.

**Proton-PottsMPNN** (third pH predictor, local CPU):
    HBPLUS_PATH=/tmp/HBPLUS/hbplus/hbplus DYLD_LIBRARY_PATH=/opt/homebrew/opt/libomp/lib \
      /tmp/ProtonPottsMPNN/.venv/bin/python /tmp/ppm_batch.py <list.txt> <out.csv>
  `/tmp` IS NOT PERSISTENT. Rebuild: clone github.com/christian-creator/ProtonPottsMPNN,
  `./install.sh`, build HBPLUS from github.com/RomanLas/HBPLUS (`make` in `hbplus/`).

**Standing traps, all of which have bitten:**
  - **a fold job is not a scored job.** ESMFold2 writes `*_ipsae.json` (raw PAE); the ipSAE
    step that writes `*_10_10.txt` is separate, and EVERY reader in `bin/` parses the latter.
    8 control shards sat "done" and unscored for hours because of this. After any fold job:
    `bin/ipsae_min.py --dir <rundir>`.
  - **names are not identities.** The same molecule appears under up to 3 run names, so a
    name-keyed `n` read 1 where 11 poses existed. Join on the binder SEQUENCE. Broke 5
    analyses; `master_rank.py` exists so it cannot happen a 6th time. A 6th instance did
    surface today and the CHECK caught it — one of `check_discards`' 45 "rejections" is an
    alias of a submitted design (METHODS §10).
  - novelty: **high structural similarity ALONE (TM >= 0.80) is Level 2**, no sequence
    condition. Missing that clause produced the bogus "0/120 pass novelty". 114 of 238 clear.
  - zsh does NOT word-split unquoted variables (broke 3 scripts)
  - chain order is NOT consistent: 29 of 69 run dirs have the TARGET in chain A. Use
    `ph_gate_all.orient()`. Broke 3 analyses.
  - construct numbering differs: crop H433=123, d3 H433=99, ECD H433=409. UniProt = PDB + 24.
  - `python3` has gemmi+propka; `.venv` has numpy but NOT gemmi
  - ProcessPoolExecutor cannot fork from a heredoc -- write a real .py file
  - score-esmfold2.sh cd's into biomodals/, so OUT_DIR must start `../`
