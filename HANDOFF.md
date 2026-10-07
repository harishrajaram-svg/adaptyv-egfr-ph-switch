# HANDOFF — Adaptyv challenge 1 (EGFR pH-switch)
# Written 2026-10-04 ~3:40 PM EDT. Supersedes the 2:00 PM / 3:10 PM revisions entirely.

**Deadline: Tue Oct 6, 23:59 AoE = Wed Oct 7, 07:59 EDT.** *(Weekdays corrected 2026-10-05: Oct 6 2026 is a TUESDAY and Oct 7 a WEDNESDAY. The dates and times were always right; both weekday labels were wrong, which is the kind of error that makes someone plan a final upload a full day late.)* Confirmed from an organiser
message, not from our own notes. Adaptyv extended it; the original Oct 4 date is wrong
everywhere it still appears.

## 🧭 PROBLEM 2's PLAN IS `ROADMAP.md` IN THE OTHER REPO

`~/code/context-directory/projects/anthropic-adaptyv-2026/ROADMAP.md` (2026-10-06 11:20 PM) is the
live plan: the 20-slot allocation, the wave, what the dry run established, the order of work, and
what is decided vs open. Read it before anything in this file.

**What a reader of THIS repo needs:** `bin/plan_wave.py` prints the wave split and launches only
under `--launch`. `biomodals/` is gitignored — the wrappers survive only via
`patches/modal_mosaic.patch` and `patches/modal_esmfold2.patch`, so regenerate after every edit.
`bin/mosaic_selftest.py` and `bin/esmfold2_selftest.py` run with bare `python3`, no GPU, and the
launchers refuse to start if they fail.

## 🔵 THIS FILE IS PROBLEM 1 ONLY. Problem 2 is live and lives elsewhere.

**Problem 1 is SUBMITTED** — 16 designs, 10/5 12:49 PM, DESIGNATED. **RESUBMISSION IN
FLIGHT 10/6:** all 18 re-uploaded after the organisers recalibrated nanobody novelty; the
final count is whatever their check clears (16, 17 or 18). See METHODS limitation 37.

**Problem 2 (conditional TNF-α binder: bind at pH 7.4, release at pH 6.0, plus mouse
cross-reactivity) closes Mon 2026-10-12 07:59 EDT**, retry-preserving upload Sun 10/11. Its
decision record is **not in this repo** — it is
`~/code/context-directory/projects/anthropic-adaptyv-2026/challenges/02-tnf-alpha.md` §1–§28,
with the current-state block in `SESSION-HANDOFF-2026-10-06.md`. Read those, not this file.

What a reader of **this repo** needs to know, because the code is here:

* **`biomodals/modal_mosaic.py` is now problem-2 capable** and `biomodals/` is **gitignored** —
  the wrapper survives only via **`patches/modal_mosaic.patch`** (1122 lines). Regenerate it
  after every edit; it has gone stale once already and nearly lost a day's work.
* It takes `--mechanism {acid_near_his,his_near_cation}`, a comma-separated `--target-chain`, and
  `--anchor-chain`. `p2trimer02` passed step 4 against **all three TNF protomers** (532 tokens).
  Cost: **21.7 s/step at 532 tokens vs 3.16 s at 228** — quadratic, so a light wave is ~$20.
* 🔴 **`HisNearCation` reduced with `score.sum()` and was gameable** — it rewards histidine mass
  near the cation, not histidine placed, and returned **15.8% histidine** while `his_best` read
  **exactly 0.00 every step**. Fixed to **top-2** via module-level `his_reduce()`, plus
  `cap_his=0.08`. `bin/mosaic_selftest.py` pins it with a **mutation test** that fails if the
  reduction is reverted. **`AcidNearHis` (problem 1) sums too and is deliberately NOT fixed** —
  nothing shipped came from a Mosaic trajectory, so no problem-1 result is affected.
* **`bin/mosaic_selftest.py` runs with bare `python3`, no numpy, no GPU, no Modal account.**
  `bin/design-mosaic.sh` calls it first. Keep it dependency-free — `bin/gate_sweep.py:32`
  records what hardcoding `.venv/bin/python` cost: two gates dead in every fresh clone.
* ⚠️ **Do not launch a Mosaic production wave.** The target and anchor are an open question for the
  reviewer; R108 is in the smoke command only because the smoke needed something.

## ⚠️ READ THIS FIRST — THE RULES LIVE IN A CHANNEL WE ARE NOT IN

**`#anthropic_adaptyv_competition` (`C0C4VEG57HU`) on the Proteinbase Slack.** Harish is **not a
member**; it is public and searchable. **Every organiser clarification since 2026-09-28 is
there**, and neither repo referenced it until now — which is why the facts below were missing
from this file for a week. Organisers: Tudor-Stefan Cotet, Simon Dürr, Amir Shanehsazzadeh.
`#design-methods` is dead for this challenge; `#feedback` carries a few answers.

**Submission mechanics, from organiser messages:**
  * **One submission per 24 hours.** Submissions are **retained, not replaced** — you nominate
    which one counts, or the most recent is designated by default. So an early upload costs
    nothing. **The last upload that still permits a second attempt is Tue Oct 6, 07:59 EDT.**
  * Novelty runs **at upload**, in 3–5 minutes, and a **self-service novelty pipeline ships
    Oct 5** so designs can be checked before committing. A competitor has reproduced our exact
    VHH failure mode on the live platform and Adaptyv have said the antibody threshold is being
    re-tested, so **the bar our two VHH rows sit near may move before the deadline.**
  * Track 3 needs nothing beyond the CSV on Proteinbase; track is assigned by account email.
    The cap is 20 designs per account per collection. **We ship 16** (updated 2026-10-06; this line said 10 and the set has been 16 since d63c61f).
  * **Iterating on any previously submitted design is explicitly disallowed** — stricter than
    the challenge page's "existing binder" wording. Our submission is clean on it, verified.

**Measurement spec, which retires a decision recorded below:**
  * human EGFR at **pH 6.5 and 7.4**; mouse EGFR at **pH 6.5 only**. Mouse is NOT 7.4-only, and
    the cross-reactivity/pH tension that premise created does not exist.
  * Targets are the **full ectodomains, tethered**: human **Met1–Ser645**, mouse **Met1–Ser647**,
    Sino Biological 10001-H08H and 51091-M08H. Met1–Ser645 minus the 24-residue signal peptide
    is **621 residues** — which is the construct we folded against.
  * You do **not** prepend the initiator Met; they add it when building constructs.
  * The assayed target carries a **His tag**, and the organisers have said a binder engaging it
    "might look pH-selective but would bind to anything with a His tag", and that in-silico
    evidence will be weighted more heavily to catch it. **All 10 of our designs switch on the
    target's native H433** (48 of 50 pool-wide; zero tag). Say so.

**The two blockers from the last handoff are CLOSED.** The methods document is finished and
the git remote exists. What is left is upload.

---

## 1. STATE — everything below is verified, not remembered

**Submission: `submissions/01-egfr.csv`, 18 designs** (the resubmission set). Track 3
allows 20 — see METHODS §11.
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
<https://github.com/harishrajaram-svg/adaptyv-egfr-ph-switch> — 200 on the README, the methods
doc, the CSV. The README is a submission front door, not the old setup log (that moved to
`docs/setup-notes-2026-09-18.md`).

Every pH number is the MEDIAN over every ESMFold2 refold pose of that exact binder sequence,
pooled across runs, joined to affinity BY SEQUENCE in `analysis/01-egfr/master_rank.json`.
**That file is the only source for both columns. Do not re-derive either from a run name.**

---

## 2. THE RESULT THAT CHANGED TODAY — control recovery

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

## 3. DECISIONS MADE TODAY — do not relitigate

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

## 4. WHAT IS LEFT — in order

1. **UPLOAD, and upload EARLY enough to resubmit.** Two of ten rows are VHH format and clear
   only on Adaptyv's antibody novelty branch. `h370_020_vhh` is antibody-rule Level 4 but
   general-rule Level 2 (TM 0.914); `rimA02_d3_rimA_14_vhh` is general-rule Level 1. **If their
   ANARCI does not call them antibodies, those two rows fail the gate at upload.** METHODS §9
   states this. Do not upload at the wire.
2. **Send the outbox.** All three of the reviewer's asks are done and still unsent:
   `outbox/ipsae-fixtures/` (8 cases, `run_fixtures.py --check`, VHH-zero trace),
   `outbox/CONTROL-TABLE.md`, `outbox/PREREGISTRATION.md`. §4.4 is the answer to his
   control-recovery ask and he should see it.
3. Optional: `master_rank.json` now has molecules with `ratio_n >= 5` and affinity that were
   never considered for the submission. Nothing in the top 10 is at risk, but the `sd`/`sd2`
   arms are the place a better design would hide.

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
