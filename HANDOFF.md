# HANDOFF — Anthropic × Adaptyv 2026, Track 3 (all challenges)
# Updated 2026-10-07 9:25 AM EDT. Covers BOTH challenges; problem 1 is closed and archived below.

**Problem 2 is the live work. Due Mon Oct 12, 7:59 AM EDT.** 20 slots, Track 3.

**Problem 1 is CLOSED** — submitted, 16 designs, deadline Wed Oct 7 07:59 EDT now passed. Nothing
in this file is a problem-1 action any more; its history is kept below because the lessons carry
forward (`lessons-problem-1.md`), not because anything is pending. **The repo is still named
`adaptyv-egfr-ph-switch` after problem 1**, which is worth renaming or noting before submission,
since the methods document links it.

## 🟢 NOTHING RUNNING as of 2026-10-07 9:15 AM — the search is finished

All 16 trajectories complete, `modal app list` empty, ~$36 spent, **the $49 wave unspent.**

**Result: the free footprint changed nothing.** Four-run rise +0.004 ± 0.013 over 75 soft steps —
flat, like the other twelve. Endpoints worse: mean 0.168 pinned → 0.135 free. **Cause 1 of
METHODS §6.3 is eliminated**; causes 2 and 3 survive and neither is testable in the remaining time.

```
cd ~/code/adaptyv-2026
python3 analysis/02-tnf/score_free_probe.py        # the full reading, verdict included
python3 analysis/02-tnf/design_inventory.py        # what is shippable: 25 seqs, 0 clear the gate
python3 bin/check_p2_stats.py                      # every methods figure, recomputed
```

**Five conditions tried. Best anywhere 0.200 against a 0.45 bar, from the second thing tried.**

| # | condition | steps | median | max |
|---|---|---|---|---|
| 1 | baseline | 50 | 0.136 | 0.156 |
| 2 | **pH weight down** | 50 | **0.179** | **0.200** |
| 3 | binding weights up | 50 | 0.156 | 0.197 |
| 4 | step budget up | 200 | 0.163 | 0.189 |
| 5 | free footprint | 100 | 0.134 | 0.137 |

**Keep this finding:** the pinned epitope was *load-bearing for the pH objective*, not an obstacle.
Freeing it pushed the median nearest-histidine distance from 22.8 Å to 34.8 Å. The anchor is one
pinned point on a 456-residue surface and the epitope term was the only thing holding the binder
near it.

### 🎯 Blocked on one decision: ROADMAP item 10

**Recommended: plan C**, the methods-first submission — 20 of 25 real sequences in two footprint
families, plus the measured negative result. The organisers confirmed on Oct 6 that *"Claude will
review what you share and use it for selection"* and that the methodology text is linked to a
public Proteinbase Collection.

**One unevaluated alternative:** Ken Osumi's structural carryover
(`github.com/ken-osumi/ArcRefine`), 6/10 TNF-α designs bound vs 1/10 unoptimized, BLI-measured at
Adaptyv, built on Boltz-2. But ~$677 for their campaign, PolyForm Noncommercial, the author's own
caveat that the comparison does not isolate the mechanism or establish the binding site, and no
statement that it rescues non-binders. Read before deciding item 10.

### Organiser facts confirmed 2026-10-06, worth not re-deriving

- The assayed target is **Asp143**, UniProt **P01375**, reagent **AcroBiosystems TNA-H4211**. Our
  target already carried Asp143 — confirmed, not corrected.
- The novelty gate is **"at least 3 of 4 checks"**, and the organisers are reviewing whether to
  relax it to 2 of 4. Watch before item 13.
- Ingmar's seven TNF-α questions (Oct 6, 11:09 AM) still have **zero replies**, including whether
  selection models against the full trimer or a single chain, and whether the trimer stays intact
  at pH 6.0. **Our §7 and §9 contingencies stand.**

### Recovery, if ever needed

```
mkdir -p /tmp/rec && cd /tmp/rec
modal volume get mosaic-weights runs/p2free-p     # destination must be the CWD; an argument raises Errno 21
```

## 🧭 THE PLAN LIVES IN `ROADMAP.md` IN THE OTHER REPO

`~/code/context-directory/projects/anthropic-adaptyv-2026/ROADMAP.md` (updated 2026-10-07 9:15 AM)
is the live plan: the result, the one outstanding decision, the 20-slot allocation, and what blocks
what. Read it before anything in this file. Then `challenges/02-tnf-alpha.md` §1–§44 for the
decision record.

**What a reader of THIS repo needs, because the code is here:**

* **`biomodals/` is gitignored.** `modal_mosaic.py` survives *only* via
  `patches/modal_mosaic.patch` (1,636 lines) and `modal_esmfold2.patch` (353 lines).
  **Regenerate after every edit** — it has gone stale once and nearly lost a day. Verified current
  as of 9:20 AM: regenerating produces a byte-identical file.
* **Restore it in a clone with:**
  `git apply --directory=biomodals patches/modal_mosaic.patch` and the same for `modal_esmfold2`.
  Both self-tests then pass with bare `python3`, no venv, no GPU, no Modal account.
* **`runs/` is gitignored too.** The 25 candidate sequences are committed separately as
  `outbox/02-tnf-candidates.faa` and `.tsv`, with each sequence's condition, score and histidine
  distance in its own FASTA header. Nothing else in git holds them.
* **Every tool carries a `--selftest`** and all eleven pass. The launchers refuse to start if the
  wrapper self-test fails.

## 🔧 TOOLING ADDED 2026-10-07 — one command each

```
python3 analysis/02-tnf/score_free_probe.py    the free-footprint reading, verdict included
python3 analysis/02-tnf/design_inventory.py    what is shippable, counted from the run tables
python3 analysis/02-tnf/loss_traj.py <logs>    block-averaged trajectory, the decisive reading
python3 bin/check_p2_stats.py                  recompute every METHODS figure from designs.tsv
python3 bin/check_no_verbatim.py               reviewer anonymity across the published tree
bin/probe-free-footprint.sh <tag> <len> <seed> the free-footprint launcher, reasoning in its header
```

`check_p2_stats.py` exists because two cells in METHODS §6.2 reported a row value where a median
belonged. It has since caught the same error class twice more, including one committed hours after
the entry describing it was written. **Run it before any edit to §6 or §10 ships.**

## ⚠️ ORGANISER FACTS — confirmed in the Slack, do not re-derive

Harish **is** a member of `#anthropic_adaptyv_competition` (`C0C4VEG57HU`) and the other four
Proteinbase channels. Organisers: Tudor-Stefan Cotet, Simon Dürr, Théo Jalabert, and Amir
Shanehsazzadeh (Anthropic).

**Problem 2, confirmed 2026-10-06/07:**
* **The assayed target is Asp143**, UniProt **P01375**, reagent **AcroBiosystems TNA-H4211**. Our
  target already carried Asp143 — organiser confirmation of a prior decision, not a correction.
* **Selection reads the methodology text.** Amir, 10-06: *"Claude will review what you share and
  use it for selection."* Tudor added that the section is **linked to a public Proteinbase
  Collection** — so it is published, which is why the reviewer-anonymity check matters and why the
  document carries no reader-directed imperatives.
* **The novelty gate is "at least 3 of 4 checks"**, and the organisers are reviewing whether to
  relax it to 2 of 4. Watch before the novelty step.
* **Model refusals are hitting other entrants** on protein-design tasks. Amir: Sonnet 5, Opus 5 and
  Sonnet 5.5 carry the same classifiers, and he recommends Sonnet 5.5. We have not been blocked.

**Still unanswered, and our contingencies depend on it.** A competitor posted seven TNF-α questions
on 10-06 at 11:09 AM with **zero replies**, two of which are ours: whether the selection step models
against the full trimer or a single chain, and whether the trimer stays intact at pH 6.0.
**METHODS §7 and §9 carry the contingencies; they stand.** Also open: the exact mouse TNF-α
sequence and vendor, and whether the expression system and C-terminal tag match challenge 1.

## 🔭 ONE UNEVALUATED ALTERNATIVE

Mosaic's author pointed at **Ken Osumi's "structural carryover"**
(`github.com/ken-osumi/ArcRefine`): **6 of 10 optimized TNF-α designs bound vs 1 of 10
unoptimized**, measured by BLI at Adaptyv, built on Boltz-2 which we already run. Relevant because
1-of-10 means their starting designs mostly did not bind either — which is our position.

**Not a recommendation yet.** The author's own caveat: candidates "were selected independently for
each group… They do not isolate the effect of structural carryover alone or establish that binding
occurs at the predicted site." Their campaign cost **~$677**, roughly 14× our unspent budget.
PolyForm Noncommercial licence, and a tool-licence question in `#design-methods` is unanswered.
Nothing states whether it rescues non-binders. **Read before deciding ROADMAP item 10.**

## 📛 THE REPO WAS RENAMED 2026-10-07 — and problem 1's artifacts deliberately were not updated

`adaptyv-egfr-ph-switch` → **`anthropic-adaptyv-2026`**, because it now covers both challenges and
the methods documents link it publicly. GitHub serves a **301 redirect** from the old URL
(verified: resolves 200), so nothing that references the old name is broken.

🔴 **Four files still carry the old name on purpose. Do not "fix" them.**
`submissions/01-egfr.csv`, `submissions/01-egfr-METHODS.md`, `bin/emit_submission_csv.py` and
`outbox/PREREGISTRATION.md`. The first two are **submitted artifacts** and the third's job is to
**reproduce the submitted CSV byte-identically** — a claim published in problem 1's methods and in
problem 2's methodology box. Changing the URL in the emitter would break that claim, because the
submitted CSV contains the old string in all 16 assessment fields. Re-verified after the rename:
the emitter still reproduces the submitted file byte-identically. The pre-registration is
immutable by definition.

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
