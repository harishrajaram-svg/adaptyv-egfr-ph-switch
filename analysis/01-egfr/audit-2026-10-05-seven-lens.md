# Seven-lens deliverable audit — 2026-10-05

14 agents (7 finders, 7 adversarial verifiers), 0 errors, 33 min, 2.71M tokens.
118 findings raised, **107 survived** verification (40 HIGH), 11 refuted.

Verified by hand before publishing this file: clusters A, B, C(rank 4), D, G.

| cluster | findings | HIGH |
|---|---|---|
| C. GRADED CSV assessment strings carry withdrawn/superseded/unsourced numbers | 12 | 9 |
| I. other | 35 | 7 |
| A. README/HANDOFF still describe a 10-design submission on the old basis | 14 | 6 |
| B. §12 attestations list a removed design, omit shipped ones | 9 | 5 |
| E. Asn420 glycan: 3 designs contact it, documents say 1 | 7 | 4 |
| G. CONTROL-TABLE still calls censored molecules "measured negatives" | 7 | 4 |
| D. §11.7 perturbation section contradicts its own artifact | 16 | 4 |
| H. §1 calls a 5-row table a full census of 17 histidines | 3 | 1 |
| F. §10b/§11 computed on a 17-design set, presented as the submission | 4 | 0 |

## HIGH

### [pk-coverage] README's four headline claims re-assert four claims PK ordered withdrawn and METHODS retracts
- `/Users/harish/code/adaptyv-2026/README.md`:21
- verdict: CONFIRMED
- quoted: 2. **A molecule with no reported KD reads a 5.27× switch.** Of 11 molecules Adaptyv measured on this / platform, the highest-scoring one on our own ranking metric is a design already measured / **not to bind** — and it ranks 8th of 2,009 on the competition's primary objective.
3. ... it turns out to be the only specificity filter here that / measured data supports ... this is a statement / about r
- artifact: All four are withdrawn elsewhere. (a) PK 2026-10-04 §1 + METHODS:367-371 — the ten are right-censored 'no KD reported', not 'measured not to bind'; and METHODS:701 lists 'a rank among 132 rankable molecules quoted against a denominator of 2,009 | overstated by ~15×'. (b) PK 2026-10-04 and METHODS:432-433 — 'no specificity filter in this pipeline is supported by measured data — including this one';
- matters: README is the repository front door ('Start here') and these four lines are the submission's summary of itself. A grader or PK reading only them gets four claims the methods document spends sections retracting, including the two PK objected to most directly (censoring and the dual-species rescue). Ironically README:111-117 prints the censoring correction 90 lines below claim 2.
- gate: Add README.md to check_claims.py's citation-window scope with the existing banned-phrase rules ('measured not to bind', 'of 2,009', 'route past it is closed', 'no cost in predicted affinity'), and generate the 'four lines' block and the design count from the emitted CSV the way METHODS §11 blocks al

### [pk-coverage] Three submitted designs contact the Asn420 glycan sequon; METHODS, README and limitation 14 report only one
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1011
- verdict: CONFIRMED
- quoted: | **glycan** | **1 of 17** touches an N-glycosylation sequon — and it is the top-ranked design |
... **The glycan flag is on rank 1.** `c5_cf_short__boltzgen_egfr_cropfree_short_48` contacts **Asn420**
- artifact: analysis/01-egfr/finalist_footprints.json (committed; regenerates byte-identically) reports 18 designs and three with glycan_sequon_hits = [420]: c5_cf_short__boltzgen_egfr_cropfree_short_48 (rank 1), bcr_d3acid3_l60_s647537_mpnn3 (rank 8) and bcr_d3acid3_l60_s647537_mpnn11 (rank 9). bin/finalist_footprints.py prints 'across 18 designs: ... 3 touch a glycosylation sequon'. The same section's other
- matters: This is an undisclosed liability on two shipped rows. The two bcr_d3acid3_l60_s647537 designs were added on 2026-10-05 as the diversity replacements, and their CSV assessment strings list novelty, pose spread and 'Expression QC NOT RUN' but say nothing about a glycan contact — while the organisers' own answer (reference/organizer-answers-slack.md:49-51) is that screening uses 'glycosylated, tether
- gate: Make §10b a GENERATED block in gen_methods_submission.py, sourced from finalist_footprints.json, emitting the per-design glycan hit list by name; and add a glycan_sequon_hits column (or an explicit sentence) to each CSV assessment string so a sequon contact cannot ship undeclared.

### [pk-coverage] The graded CSV's rimA01_r15_L133E assessment carries five values METHODS §11.6 withdraws by name, contradicting its own numeric columns
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr.csv`:15
- verdict: CONFIRMED
- quoted: all-site median 5.659x, above the 5.55x thermodynamic ceiling ... Two of five poses exceed the 7.94x one-proton bound ... Affinity HELD: human 0.594 -> 0.616 ... the pose spread is 4.24-11.64, i.e. 1.31x the median ... L133D at the same position reads 0.723x with human affinity 0.000 -- Asp spans ~2.5A from CB and LEU133 sits 6.28A from H370
- artifact: The same CSV row reads ph_poses_n=20, ph_pose_spread_over_median=4.38, ipsae_min_human=0.5978, ph_ratio_allsite_SENSITIVITY=52.181, his-only=5.656. analysis/01-egfr/ph_sensitivity.json: n_poses 20, hisonly_median 5.6558, hisonly_min/max 3.8379/28.6292, allsite_median 52.1807. METHODS:1365-1377 withdraws each item explicitly: 'At n = 20 the human leg is **0.598** against the parent's 0.594 — **flat
- matters: The CSV is THE graded upload and the assessment column is the only prose a grader sees. One row simultaneously reports spread 4.38 in its numeric column and '1.31x the median' in its text, and asserts an affinity gain the methods document calls five-pose noise. All nine gates pass (bin/gate_sweep.py: 9/9 green), so nothing catches it. It also mislabels the ranked basis: 5.659 is quoted as 'all-sit
- gate: Extend emit_submission_csv.py --selftest to assert that every number appearing in an assessment string also appears in that row's own numeric columns or in the row's ph_sensitivity.json entry (tolerance 0.01), failing on any orphan figure; and add the retired values (5.659, 1.31x, 0.616, 0.723x) to 

### [pk-coverage] METHODS §5 uses the EGF-derived controls' Potts ranks as validation and asserts their agonist activity — the exact inference PK ruled out
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:652
- verdict: CONFIRMED
- quoted: The controls validate it independently of us: all six positive controls fall in the lower
half of the design distribution, and both EGF-derived "nonbinders" fall in the **top 14%** —
a model that never saw this target says protonating H433 is maximally bad for an EGF-like
complex, which is exactly right for a neutral-pH agonist.
- artifact: PK 2026-10-04 §5, final line: 'The EGF-derived sequences' Potts ranks cannot independently establish their binding or agonism.' And PK §1: 'Relabel both EGF-derived controls as activity-unknown, document their provenance gap, and withdraw claims based on their supposed negative status.' METHODS §4.1 (lines 310-316) complies: 'Both EGF-derived controls are therefore relabelled **activity-unknown**.
- matters: This is the one paragraph that still treats the EGF-derived molecules as a known agonist and uses their Potts rank as corroboration of the protonation model — the single inference PK named and forbade, in a section whose whole purpose is the third pH predictor's credibility. It also contradicts §4.1 eleven sections earlier and describes the panel's scare-quoted 'nonbinders' as validating 'independ
- gate: Add a check_claims.py rule binding the tokens NEG_nonbinder / NEGd3_nonbinder / 'EGF-derived' to a required 'activity-unknown' qualifier within the same paragraph, and ban the words 'validate'/'agonist' in any sentence citing their scores.

### [pk-coverage] §12 Declarations attest 17 sequences, name one that was removed, and omit two that shipped
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1702
- verdict: CONFIRMED
- quoted: **Five sequences were added on 2026-10-05 and their review is recorded separately, because it is
a human attestation and must not be inflated by restating a count.** They are
`c5_cf_short__boltzgen_egfr_cropfree_short_48`, `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05`,
`sd_d2c_101_l147_s144898_m_T65D`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA` and
`cons_gap_h370_only__boltzgen_egfr_h370_018`.
...
*
- artifact: submissions/01-egfr.csv holds 18 rows. `sd_d2c_101_l147_s144898_m_T65D` is not one of them — METHODS §10 and §11.3 record that it was removed as a near-duplicate and replaced by `bcr_d3acid3_l60_s647537_mpnn3` and `_mpnn11`, neither of which appears in this attestation. 12 original + 5 − 2 removed + 2 replacements + 1 restored = 18, so the 12/17 arithmetic in §12 predates the swap.
- matters: §12 is the declarations section the organisers require, and it states in its own words that a human attestation 'must not be inflated by restating a count'. As written, two shipped sequences are covered by no human review and no coded provenance check, and the document attests to a sequence that is not in the upload. The 'checked by code over all seventeen' provenance claim (no G532 ancestry, max 
- gate: Make the §12 attestation list a GENERATED block keyed on the emitted CSV, with gen_methods_submission.py --check failing if any CSV name is absent from the review lists or any listed name is absent from the CSV.

### [pk-coverage] CONTROL-TABLE still calls the ten control molecules 'measured negatives' four lines after printing the correction that they are not
- `/Users/harish/code/adaptyv-2026/outbox/CONTROL-TABLE.md`:363
- verdict: CONFIRMED
- quoted: (Corrected 2026-10-05: a missing KD is right-censoring, not a measured zero, so this
  separates one quantified binder from ten censored observations.)
  The highest-scoring molecule in the entire measured panel is a measured **non-binder** —
  `EXPNEG_gitter-yolo10` at 0.5893 against EGF's 0.3549.
[and §2, line 29:] **10 designs expressed, tested against EGFR on this platform, NO binding detected
- artifact: PK 2026-10-04 §1 required the opposite framing, and the same file's §6 header already applies it: 'NO KD REPORTED (right-censored; affinity weaker than the quantifiable limit, or an expression/QC failure)'. METHODS:366-371 states it fully. data/proteinbase/egfr_round1_second.csv reports no KD for these ten; it does not report non-binding.
- matters: PK named this file as one of the two prerequisites for the next decision ('Please send the scorer/fixtures and the new control table first'). It is unsent, and in its current state it answers his central correction by restating the error twice — once in the section that introduces the panel and once in the sentence that carries the panel's headline result. §2's 'These are *measured* negatives, not
- gate: Run check_claims.py over outbox/ as well as submissions/, with 'measured negative(s)', 'measured non-binder' and 'NO binding detected' banned outside a quoted-correction block.

### [arithmetic] §11.7 and limitation 17 both assert no design holds top-3 in >50% of draws; the artifact says 60%
- `submissions/01-egfr-METHODS.md`:1492
- verdict: CONFIRMED
- quoted: **No design holds a top-three position in more than 50% of draws.** The conclusion does not depend on σ: it already holds at the optimistic 0.4.
- artifact: analysis/01-egfr/ph_pka_perturbation.json (sigma 0.8, draws 400, 18 designs): rimA01_r15_L133E has top3_fraction = 0.6. Next highest c5_cf_short = 0.487, rimA02 = 0.427.
- matters: This is the load-bearing sentence of the submission's headline epistemic claim ("the pH ratio cannot order this submission"), repeated verbatim in limitation 17 at line 1817. The shipped artifact contradicts it: the top-ranked design does hold a top-three slot in a clear majority of draws. A grader who opens ph_pka_perturbation.json finds the paper's own file refuting the paper's own conclusion, w
- gate: Add a check_claims rule that reads top3_fraction from ph_pka_perturbation.json and asserts max(top3_fraction) < 0.5 whenever the >50% sentence is present; better, generate the sentence (and limitation 17) from the artifact's max top3_fraction and n.

### [arithmetic] "Designs move by up to 8 ranks" understates the measured maximum by a third, and the exemplar's two ranks are both wrong
- `submissions/01-egfr-METHODS.md`:1429
- verdict: CONFIRMED
- quoted: Designs move by up to **8 ranks**
(`rimA02_d3_rimA_14_vhh`: 2nd on his-only, 10th on partnered).
- artifact: Recomputed from ph_sensitivity.json over the 18 shipped designs, and confirmed by the CSV's own ph_rank_range_across_bases column: rimA02_d3_rimA_14_vhh is 3rd on his-only and 13th on partnered (CSV reads "3-13"), a 10-rank move; the largest move is bc_d3acid_l65_s831683_mpnn11, rank 5 to rank 17 (CSV reads "5-17"), 12 ranks. The 2nd/10th/8-rank figures are the 12-design values from before the 202
- matters: The sentence is the quantitative statement of how unstable the ranking is, and the graded CSV carries a column that contradicts it in the same submission package. The error runs in the submission's favour (the table looks a third more stable than it is), and the same "up to 8 ranks" is repeated in README.md line 97.
- gate: Compute max rank span from the emitted ph_rank_range_across_bases column and emit the "moves by up to N ranks" clause plus its worst-case exemplar as a generated block; the data for it is already in the CSV.

### [arithmetic] §11.7's "top set" of five names a design that is not in the submission and omits one that qualifies
- `submissions/01-egfr-METHODS.md`:1499
- verdict: CONFIRMED
- quoted: five designs (`rimA01_r15_L133E`, `c5_cf_short…_48`, `sd_d2c…_T65D`,
`rimA02_d3_rimA_14_vhh`, `rimA01_r15`) hold a top-three slot in 25–50% of draws and the rest
hold it in 0–18%.
- artifact: `sd_d2c_101_l147_s144898_m_T65D` appears in ph_sensitivity.json but is absent from submissions/01-egfr.csv and absent from ph_pka_perturbation.json entirely — §10 line 977 records that it was removed as a near-duplicate. In the artifact the 25–50% band holds four designs: c5_cf_short 0.487, rimA02 0.427, d2c_mpnn13_S88D_serasp 0.275, rimA01_r15 0.263; rimA01_r15_L133E is 0.600, i.e. outside the st
- matters: This sentence is the only positive claim §11.7 makes ("a top set exists"), and the named set is wrong in three ways at once: it includes a molecule that was not submitted, excludes d2c_mpnn13_S88D_serasp which does fall in the band, and includes L133E whose fraction sits above the band's upper edge. A grader told which five designs are the defensible top set would be looking for one that is not in
- gate: Derive the top-set membership list from ph_pka_perturbation.json filtered to names present in the emitted CSV, and fail the gate on any design name cited in METHODS §11 that is not a row of submissions/01-egfr.csv.

### [arithmetic] The generated binder-histidine block credits d2c_mpnn13_S88D_serasp with 14 of its own histidines; the sequence has 2
- `submissions/01-egfr-METHODS.md`:1144
- verdict: CONFIRMED
- quoted: **nine of the eighteen submitted designs carry at least one histidine of their own**: `d2c_mpnn13_S88D_serasp` (14); `ss_bc_s831683_mpnn6_S15D_S62H_routeA` (4); ...
- artifact: The shipped sequence for d2c_mpnn13_S88D_serasp contains 2 histidines (positions 42 and 96). ph_sensitivity.json gives n_his = 19 for this complex because it is the one design whose target is the full 621 aa ECD, which itself carries 17 histidines; every other design's target is the 170 aa crop with 5. The generator evidently reports n_his − 5. METHODS line 1569 and the CSV's own d2c assessment bo
- matters: This is inside a GENERATED block that gen_methods_submission.py --check passes, so the gate re-derives the same wrong number and certifies it. The block exists precisely because this count was wrong five times by hand (line 1148 says so), and it is still wrong — by 7× on the one design the error applies to. It also inflates the stated severity of the binder-histidine desolvation drag, which is §11
- gate: Count histidines in the CSV `sequence` field directly (seq.count('H')) instead of subtracting a hard-coded 5 target histidines from ph_sensitivity.json's n_his; add a self-test asserting the two agree for every row.

### [arithmetic] Three of eighteen designs contact the Asn420 glycan sequon, not one of seventeen
- `submissions/01-egfr-METHODS.md`:1011
- verdict: CONFIRMED
- quoted: | **glycan** | **1 of 17** touches an N-glycosylation sequon — and it is the top-ranked design |
- artifact: finalist_footprints.json: glycan_sequon_hits = [420] for c5_cf_short__boltzgen_egfr_cropfree_short_48 AND for bcr_d3acid3_l60_s647537_mpnn3 AND for bcr_d3acid3_l60_s647537_mpnn11 — 3 of 18 designs, all on the same sequon (mature N420, contacted by 3/18 designs).
- matters: The glycan check is one of the four footprint checks PK asked for by name, and the section's follow-up ("The glycan flag is on rank 1", line 1017), limitation 14 (line 1801) and README line 172 all present it as a single-design liability. Two further shipped designs — ranks 8 and 9, both added on 2026-10-05 — carry the same liability and are nowhere flagged. The designs were scored against glycan-
- gate: Emit the glycan row of the §10b table as a generated block that counts designs with non-empty glycan_sequon_hits and names every one of them, and add a CSV column or assessment sentence for any row with a sequon hit.

### [arithmetic] The two bcr replacement designs are 86.7% identical to each other, while both CSV rows and §11.3 claim "<0.47 to anything else submitted"
- `submissions/01-egfr.csv`:9
- verdict: CONFIRMED
- quoted: Highest sequence identity to any other submitted design: <0.47.
- artifact: Recomputed pairwise over the 18 shipped sequences: bcr_d3acid3_l60_s647537_mpnn3 and bcr_d3acid3_l60_s647537_mpnn11 are both 60 aa and differ at 8 positions — identity 0.867. The same "<0.47" claim is in the mpnn11 row (CSV line 10), and METHODS repeats it as "at most 0.467 identical to anything else submitted" at lines 979 and 1275. METHODS' own generated FAMILY-LIST (line 1290) already shows the
- matters: This number is the entire stated justification for the swap: two additions were pulled for being 0.985/0.986 near-duplicates of shipped designs and replaced by this pair, declared non-duplicative at <0.47. The replacement pair is itself an 86.7% pair — the fourth-closest pair in the submission — so the fix reproduced the problem PK asked to avoid ("avoid filling available slots with nearly identic
- gate: Compute the full pairwise identity matrix over the emitted CSV sequences in the emitter and generate both the "three pairs exceed 90%" table and each row's max-identity sentence from it; fail the gate if any assessment string's stated identity bound is below the computed value.

### [arithmetic] The L133E assessment string in the graded CSV carries four figures §11.6 explicitly retracted, each contradicting its own row's columns
- `submissions/01-egfr.csv`:15
- verdict: CONFIRMED
- quoted: all-site median 5.659x, above the 5.55x thermodynamic ceiling for H433 alone ... Two of five poses exceed the 7.94x one-proton bound ... Affinity HELD: human 0.594 -> 0.616, mouse 0.567 -> 0.435 ... the pose spread is 4.24-11.64, i.e. 1.31x the median
- artifact: Same CSV row: ph_ratio_allsite_SENSITIVITY = 52.181 (not 5.659), ipsae_min_human = 0.5978 (not 0.616), ph_poses_n = 20 (not 5), ph_pose_spread_over_median = 4.38 (not 1.31). METHODS §11.6 (lines 1368–1381) retracts all four by name: "At n = 20 the human leg is **0.598** against the parent's 0.594 — **flat, not up**. The apparent gain was five-pose noise"; "On twenty it is **4.38×**, range **3.838 
- matters: The CSV is the graded upload. §11.6 was rewritten the same day specifically to withdraw these five-pose numbers, and the methods document even names the old values — but the assessment string shipped unchanged, so the graded file asserts a held/improved affinity, a 5.659× all-site median under the 5.55× ceiling, and a 1.31× spread that the row's own four numeric columns each contradict. A grader r
- gate: Have check_claims.py parse every number in each assessment string that is also an emitted column for that row (pH ratio, all-site, poses_n, spread, ipsae human/mouse) and fail on mismatch; and add a rule forbidding any value named in the §11.6 "what did not survive" list from appearing anywhere in t

### [arithmetic] Four S15D rows quote the superseded target-only swing as the design's switch, next to graded columns reading ~1.0×
- `submissions/01-egfr.csv`:19
- verdict: CONFIRMED
- quoted: Matched wild-type, same batch, same gate: 3.15x -> 5.46x ... 5.46x is 98% of H433's single-site thermodynamic ceiling (5.55x), so the four S15D designs span 0.22 pKa units and are NOT distinguishable from one another.
- artifact: That row's own graded column ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE = 1.023; the 3.15/5.46 pair is the ph_ratio_target_only_SUPERSEDED basis (5.461). Same pattern in bc_s831683_mpnn9_S15D (line 16: "3.52x -> 5.43x" against columns 0.627 → 1.062), bc_s831683_mpnn6_S15D ("3.85x -> 5.40x" against 1.835) and bc_s831683_mpnn19_S15D ("3.68x -> 5.43x" against 1.774). Only bc_s831683_mpnn9_WT (line 1
- matters: METHODS §11.1 line 1194 says four designs "are no longer switches" on the shipped basis — these are those four. Their assessment strings nonetheless present a 5.4× switch at "98% of the thermodynamic ceiling" with no basis label, while the column being graded says 1.0×. One row out of five labels the basis, which shows the fix was known and applied inconsistently. A grader comparing narrative to c
- gate: Require any ratio in an assessment string to carry its basis token (his_only / target_only_SUPERSEDED / all_site), and have check_claims.py reject an unlabelled ratio that matches the superseded column but not the graded one.

### [arithmetic] §11.2 says the VHH "leads the submission" with the "second-highest" pH ratio at "rank 6"; it is third-highest at rank 12
- `submissions/01-egfr-METHODS.md`:1234
- verdict: CONFIRMED
- quoted: On a pure pH ordering
`rimA02_d3_rimA_14_vhh` leads the submission at 4.838× ... **The cost, stated:** rimA02 carries the second-highest honest pH ratio in the submission and
sits at rank 6, below a design at 1.774×.
- artifact: The GENERATED rank table 25 lines above (lines 1211–1229) lists rimA01_r15_L133E at 5.656 and c5_cf_short at 5.546 ahead of rimA02's 4.838, and places rimA02 at rank 12; the design at 1.774 (bc_s831683_mpnn19_S15D) is rank 11. The 2nd-place / rank-6 figures are the pre-2026-10-05 twelve-design values.
- matters: This paragraph is the stated justification for the one deliberate departure from ranking on the primary objective — PK's "apply eligibility and credible-interface checks first". The justification is built on a pH standing (leads / second-highest) that the section's own generated table refutes, and on a rank that is six places off, so the "cost, stated" is not the cost actually paid.
- gate: Forbid bare rank references in METHODS §§10–11 prose (the document's own §11.4 note already says ranks drift) and generate the few that are load-bearing; or add a check that any "rank N" citation resolves to the same design in the generated RANK-TABLE.

### [unsupported] The graded CSV's rank-14 assessment is the pre-triad version: five numbers §11.6 explicitly retracted, three contradicting columns in the same row
- `submissions/01-egfr.csv`:15
- verdict: CONFIRMED
- quoted: Affinity HELD: human 0.594 -> 0.616, mouse 0.567 -> 0.435. WHY IT IS RANKED LAST: the pose spread is 4.24-11.64, i.e. 1.31x the median [...] L133D at the same position reads 0.723x with human affinity 0.000 -- Asp spans ~2.5A from CB and LEU133 sits 6.28A from H370 [...] while Glu (~3.9A) reaches. [...] it needs 10-15 seeds and an isosteric L133Q control to be a result.
- artifact: METHODS.md §11.6 lines 1365-1377 withdraws each of these after the 15-seed triad landed: (1) "'affinity held: human went up, 0.594 -> 0.616.' At n = 20 the human leg is 0.598 against the parent's 0.594 -- flat, not up. The apparent gain was five-pose noise" — and the row's own ipsae_min_human column reads 0.5978; (2) "L133D his-only '0.723x'. At n = 15 it is 2.221x"; (3) "The Glu/Asp reach argumen
- matters: Every one of these is a claim about what was measured, in the one file the organisers grade, and each is the version the project itself overturned by deeper sampling. Three of them are contradicted by numeric columns in the same CSV row, so a grader comparing the prose to the columns sees the submission disagreeing with itself about its own design — and a grader who reads only the prose is told af
- gate: Generate the numeric spans inside assessment strings rather than hand-writing them: add assessment fields to the GENERATED blocks so every figure quoted in prose (`human X -> Y`, `pose spread Z`, per-site ranges, variant ratios) is emitted from master_rank.json / ph_sensitivity.json. Minimal version

### [unsupported] §5 calls the Potts controls a validation, on reasoning the document's own §4.1 withdrew — and no control molecule appears in either published Potts artifact
- `submissions/01-egfr-METHODS.md`:652
- verdict: CONFIRMED
- quoted: The controls validate it independently of us: all six positive controls fall in the lower
half of the design distribution, and both EGF-derived "nonbinders" fall in the **top 14%** —
a model that never saw this target says protonating H433 is maximally bad for an EGF-like
complex, which is exactly right for a neutral-pH agonist.
- artifact: The two published Potts artifacts — analysis/01-egfr/potts_ddg_pool.json (1,693 records) and potts_ddg_extra.json (399) — contain no control molecule. A scan of all 2,092 labels for POS_/NEG_/expctrl/ctrl2/nonbinder/cetuximab/gitter/yolo/deepsat/rac1/4uip returns zero hits; the only 'g532' labels are this project's own boltzgen_egfr_g532mimic_* designs, not the antibody. No analysis/ file carries 
- matters: Three lines above, the document states the correct position — "Agreement between two methods is evidence, not validation. Neither has experimental ground truth on this target" (line 649). Line 652 then calls the same kind of agreement a validation, and the inference that makes it work ("which is exactly right for a neutral-pH agonist") treats the EGF-derived control's agonist activity as establish
- gate: A check_claims rule keyed to the word "validate"/"validates"/"validated": any sentence using it about a predictor must cite a file under analysis/ that contains a measured outcome for the molecules named. Operationally cheaper: require every percentile/rank claim about a named control to resolve to 

### [contradiction] README calls the submission 10 designs; the graded CSV has 18
- `README.md`:13
- verdict: CONFIRMED
- quoted: | **[submissions/01-egfr.csv](submissions/01-egfr.csv)** | The submission: 10 designs, ranked. |
- artifact: submissions/01-egfr.csv has 18 data rows (18 designs, 15 columns). METHODS §11 line 1131: "**18 designs, ranked on the two-partner histidine-only pH product.**" README repeats the 10-design figure at line 34 ("We submitted **10 designs of the 20 allowed**"), while its own line 158 says "(12 → 17 of the 20 permitted)". Three different submission sizes in one file, none of them 18.
- matters: README is the submission front door (HANDOFF line 66-68: "The README is a submission front door"). A grader comparing the front door against the uploaded CSV sees the file describing a different submission than the one graded, and cannot tell which 10 of the 18 rows the claims apply to.
- gate: Generate the design count, the slot count and the per-rank list in README from submissions/01-egfr.csv inside GENERATED markers, the same way bin/gen_methods_submission.py already does for METHODS §11; add README to its --check set.

### [contradiction] The §12 human-review declaration lists a sequence that is not in the submission and omits two that are
- `submissions/01-egfr-METHODS.md`:1702
- verdict: CONFIRMED
- quoted: **Five sequences were added on 2026-10-05 and their review is recorded separately, because it is
a human attestation and must not be inflated by restating a count.** They are
`c5_cf_short__boltzgen_egfr_cropfree_short_48`, `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05`,
`sd_d2c_101_l147_s144898_m_T65D`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA` and
`cons_gap_h370_only__boltzgen_egfr_h370_018`.
- artifact: `sd_d2c_101_l147_s144898_m_T65D` does not appear in submissions/01-egfr.csv. METHODS §10 line 977-981 and §11.3 line 1271-1276 both record that it was removed and replaced by `bcr_d3acid3_l60_s647537_mpnn3` and `_mpnn11`, which ARE in the CSV (ranks 8 and 9) and are named in neither declaration. 12 reviewed + 5 listed = 17, against 18 CSV rows.
- matters: This is the attestation block the organisers ask for. It certifies review of a sequence that was not uploaded, leaves two uploaded sequences uncovered by any review statement, and its own arithmetic accounts for 17 of 18 rows.
- gate: Generate the §12 added-design list from the CSV minus the frozen 2026-10-04 twelve, and make the generator raise when the union of the "twelve reviewed" plus the listed additions is not exactly the CSV's name set.

### [contradiction] §12 declares structure coverage as 10 of 12 with two designs missing; 8 of 18 are missing
- `submissions/01-egfr-METHODS.md`:1731
- verdict: CONFIRMED
- quoted: **Coverage is 10 of 12**: the two designs added latest, `bc_s360518_mpnn9_A22D` and
`rimA01_r15_L133E`, have no published structure yet.
- artifact: submissions/structures/ holds 10 .cif files for an 18-design submission. Eight submitted designs have no structure: c5_cf_short__boltzgen_egfr_cropfree_short_48, c5_cr_crop_patch__boltzgen_egfr_crop_patch_05, bc_s360518_mpnn9_A22D, ss_bc_s831683_mpnn6_S15D_S62H_routeA, cons_gap_h370_only__boltzgen_egfr_h370_018, bcr_d3acid3_l60_s647537_mpnn3, bcr_d3acid3_l60_s647537_mpnn11, rimA01_r15_L133E.
- matters: The line names exactly which designs lack a published pose, and it names two when the real answer is eight — including rank 1 and rank 2. A reviewer who wants the pose for the top-ranked design is told coverage is 83% when it is 56%. The sentence even flags that an earlier version of itself was wrong about this same count.
- gate: Generate the coverage sentence by diffing the CSV name column against `ls submissions/structures/*.cif`; fail the gate sweep if the stated count or the named-missing list differs.

### [contradiction] CONTROL-TABLE §2 still presents the ten Adaptyv molecules as measured non-binders — the exact relabel PK demanded and METHODS performed
- `outbox/CONTROL-TABLE.md`:28
- verdict: CONFIRMED
- quoted: * **10 designs expressed, tested against EGFR on this platform, NO binding detected.**
    48–200 aa, from two independent groups. These are *measured* negatives, not presumed ones.
- artifact: METHODS §4.4 line 364-373: "Earlier versions of this section called them 'measured non-binders'. That is not what the data says ... Throughout this section they are therefore **"no KD reported"**". CONTROL-TABLE's own §6 header (line 254) reads "NO KD REPORTED (right-censored; affinity weaker than the quantifiable limit, or an expression/QC failure)" and line 361-362 carries the "Corrected 2026-10
- matters: CONTROL-TABLE is one of the two artifacts PK named as prerequisites for the next decision (2026-10-04 reply, closing line). Its introduction of the control class asserts the exact claim he asked to be withdrawn, and the document contradicts itself four sections later. A reviewer reading §2 concludes the project has assay-confirmed negatives; it does not.
- gate: Add a check_claims.py rule banning "measured negative(s)", "measured non-binder" and "NO binding detected" applied to the EXPNEG class anywhere in submissions/ or outbox/, with the right-censored phrasing as the only permitted form.

### [contradiction] CONTROL-TABLE tells the reviewer the affinity_above_null column is still in the submission and that METHODS §4.5 agrees; §4.5 says it was removed and the CSV has no such column
- `outbox/CONTROL-TABLE.md`:80
- verdict: CONFIRMED
- quoted: We have left the column in the submission rather than dropping it mid-flight, and said the same
thing in METHODS §4.5 and §11 so all three documents now agree. If you would rather it came out,
it is a one-line change and we have ~38 hours.
- artifact: The CSV header is `name,sequence,molecule_class,ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE,ipsae_min_human,ipsae_min_mouse,ph_poses_n,ph_pose_spread_over_median,ph_ratio_allsite_SENSITIVITY,ph_ratio_partnered_SENSITIVITY,ph_rank_range_across_bases,ph_tier_provisional,ph_ratio_target_only_SUPERSEDED,affinity_assessable,assessment` — no affinity_above_null. METHODS §4.5 line 558-559: "**The `affini
- matters: This is the outgoing artifact to the external reviewer. It asks him to decide whether to remove a column that is already gone, and it cites METHODS §4.5 as agreeing with the opposite of what §4.5 says. The §4.5 pointer is a live cross-reference to content that states the reverse.
- gate: Have check_claims.py assert that any document claiming a CSV column exists names a column present in the emitted CSV header, and that any "all three documents now agree" sentence is backed by a matching-string check across the three files.

### [contradiction] README and HANDOFF still quote "rank 8 of 2,009", a figure METHODS §7 lists as a ~15x overstatement it corrected
- `README.md`:21
- verdict: CONFIRMED
- quoted: 2. **A molecule with no reported KD reads a 5.27× switch.** Of 11 molecules Adaptyv measured on this
   platform, the highest-scoring one on our own ranking metric is a design already measured
   **not to bind** — and it ranks 8th of 2,009 on the competition's primary objective.
- artifact: METHODS §7 line 701: "| a rank among 132 rankable molecules quoted against a denominator of 2,009 | overstated by ~15× (§4.4) |". METHODS §4.4 line 437-440: "(An earlier version of this document said \"rank 8 of 2,009\" ... That overstated it by about 15× and the derived \"top 0.4%\" was wrong; 2,009 is the scored universe, 132 is the rankable one.)" HANDOFF line 89-90 carries the same "rank 8 of 
- matters: This is one of the four headline claims the README offers as what the submission establishes. It states a rank the methods document formally retracts as a 15x overstatement, and describes the molecule with the status the reviewer specifically asked to be withdrawn. Both errors are in the same two sentences.
- gate: check_claims.py already has regex rules — add "of 2,009" / "of 2009" as a forbidden denominator for any rank claim, and extend the rule set to cover README.md and HANDOFF.md, not just submissions/.

### [contradiction] README headline claim 4 asserts the two-site route is closed; METHODS retracts that twice by name
- `README.md`:29
- verdict: CONFIRMED
- quoted: 4. **The single-site ceiling is 5.55× and the route past it is closed.** H433's free pKa is
   6.22, so no single-site design can beat 5.55× over a 0.9 pH-unit window. The only histidine
   pair close enough to bridge (H433+H370, 8.5 Å) is unreachable
- artifact: METHODS §8.1 line 747: "**2. The two-site H433+H370 route is demonstrated, not closed.**" METHODS §3.4 line 161-163: "**This section has been wrong twice.** It first read \"the two-site route is closed\" — it is not." METHODS §3.4 also reports rimA01_r15_L133E at 5.656x, above the 5.55x single-site ceiling, and §11.6 line 1360 reports 9 of 15 triad poses above it and 3 of 15 above the 7.94x one-pr
- matters: A grader reading the README's four-line summary takes away the opposite of the methods document's corrected position, and the opposite of the published control the external reviewer supplied. It is also the one claim the project most explicitly recanted.
- gate: Keep a retracted-claims phrase list ("route past it is closed", "two-site route is closed", "is unreachable") in bin/check_claims.py and run it over README.md, HANDOFF.md and outbox/ as well as submissions/.

### [contradiction] The CSV's rank-14 assessment cell carries three numbers METHODS §11.6 withdraws by name, plus a spread that contradicts its own column
- `submissions/01-egfr.csv`:15
- verdict: CONFIRMED
- quoted: Affinity HELD: human 0.594 -> 0.616, mouse 0.567 -> 0.435. WHY IT IS RANKED LAST: the pose spread is 4.24-11.64, i.e. 1.31x the median, against 0.04-0.90 for every other row. ... L133D at the same position reads 0.723x with human affinity 0.000 -- Asp spans ~2.5A from CB and LEU133 sits 6.28A from H370
- artifact: METHODS §11.6 line 1365-1377 withdraws all three: (1) "*\"affinity held: human went **up**, 0.594 → 0.616.\"* At n = 20 the human leg is **0.598** against the parent's 0.594 — **flat, not up**. The apparent gain was five-pose noise." — and the row's own ipsae_min_human column reads 0.5978; (2) "*`L133D` his-only \"0.723×\".* At n = 15 it is **2.221×**"; (3) "*The Glu/Asp reach argument was geometr
- matters: This is the graded upload, and the assessment column is the only prose Adaptyv reads. Five figures in one cell are the five-pose values the methods document explicitly names as not surviving deeper sampling, and two of them contradict numeric columns in the same row of the same file.
- gate: Have bin/emit_submission_csv.py assert that every numeric literal in an assessment string matches the row's own columns or a value in master_rank.json / ph_sensitivity.json, and fail on any number not found there.

### [contradiction] The CSV's rank-4 assessment calls 3.738x "the all-titratable-site basis" — the exact mislabel §11.7 exists to correct, and the row's all-site column reads 88.593
- `submissions/01-egfr.csv`:5
- verdict: CONFIRMED
- quoted: On the all-titratable-site basis (both partners, see METHODS) it reads 3.738x against a wild-type 0.574x -- swing 6.51x, i.e. the result survives the honest multi-site accounting, which most of this submission does not.
- artifact: 3.738 is the value of this row's `ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE` column; its `ph_ratio_allsite_SENSITIVITY` column reads 88.593. METHODS §11.7 line 1402-1411: "`bin/ph_gate_multisite.py` was described here and in the preregistration as composing \"every titratable site on both partners\". It does not. ... The shipped `ph_ratio_*` column is therefore a **two-partner histidine-only app
- matters: The cell claims the design's result "survives the honest multi-site accounting" using a number the methods document proves is blind to multi-site accounting — blind, specifically, to the ASP22 this design was built around. It is the one sentence in the graded file that asserts the correction was passed, and it asserts it with the pre-correction quantity under the post-correction label.
- gate: Add a check_claims.py rule that any assessment string containing "all-site", "all-titratable" or "multi-site" must quote the row's ph_ratio_allsite_SENSITIVITY value, not its his-only value.

### [arm-accountability] The Sec 11.7 pKa-perturbation table is a stale 17-design run; all six cells and the top-three conclusion fail to reproduce on the 18-design submission
- `submissions/01-egfr-METHODS.md`:1492
- verdict: CONFIRMED
- quoted: **No design holds a top-three position in more than 50% of draws.** The conclusion does not depend on σ: it already holds at the optimistic 0.4.
- artifact: `analysis/01-egfr/ph_pka_perturbation.json` contains 18 designs, and `rimA01_r15_L133E` has `"top3_fraction": 0.6`. Re-running `bin/ph_pka_perturbation.py --sigma 0.8 --draws 400` reproduces the committed artifact byte-identically and prints "18 submitted designs ... top3% 60" for that design. The σ table at lines 1485-1487 is also stale in every cell: document reads 0.4 -> "6 of 17 | 13 of 17", 0
- matters: The table is one of the document's two headline self-critical results and it is reported against a 17-design pool that no longer exists -- the 18th submitted design was added and the section was never regenerated. The error is not uniform in direction: "designs keeping their baseline rank" is reported as 4 of 17 when it is 2 of 18, and at σ=1.2 as 4 of 17 when it is 1 of 18, so the table overstate
- gate: Make the entire Sec 11.7 pKa table and both top-three sentences a `GENERATED:PKA-PERTURBATION` block in `bin/gen_methods_submission.py --write`, sourced from `ph_pka_perturbation.json`, with the σ=0.4/0.8/1.2 runs written to three distinct output paths instead of one overwritten file. Add a `check_c

### [overclaim] "No design holds a top-three position in more than 50% of draws" — one does, at 60%, and it is the design the submission leans on
- `submissions/01-egfr-METHODS.md`:1492
- verdict: CONFIRMED
- quoted: **No design holds a top-three position in more than 50% of draws.** The conclusion does not depend on σ: it already holds at the optimistic 0.4.
- artifact: analysis/01-egfr/ph_pka_perturbation.json (sigma 0.8, 400 draws, 18 designs): rimA01_r15_L133E has top3_fraction = 0.600. Next are c5_cf_short 0.487, rimA02_d3_rimA_14_vhh 0.427, d2c_mpnn13_S88D_serasp 0.275, rimA01_r15 0.263.
- matters: A blanket universal that a single row disproves, and it is restated as Limitation 17 (line 1817, "no design holds a top-three slot in more than half of draws") — the limitations list is the part a grader trusts most. The same paragraph (lines 1499-1501) then names the top set as "five designs (`rimA01_r15_L133E`, `c5_cf_short…_48`, `sd_d2c…_T65D`, `rimA02_d3_rimA_14_vhh`, `rimA01_r15`) hold a top-
- gate: check_claims RULE: assert max(d['top3_fraction'] for d in pert['designs']) against any "top-three ... more than N%" / "more than half" phrasing, and assert the named top-set membership is a subset of the CSV names.

### [overclaim] Three submitted designs contact the Asn420 glycosylation sequon, not one; §10b and Limitation 14 name only rank 1
- `submissions/01-egfr-METHODS.md`:1011
- verdict: CONFIRMED
- quoted: | **glycan** | **1 of 17** touches an N-glycosylation sequon — and it is the top-ranked design |
- artifact: analysis/01-egfr/finalist_footprints.json (regenerated 07:28, now 18 designs): glycan_sequon_hits == [420] for THREE designs — c5_cf_short__boltzgen_egfr_cropfree_short_48 (rank 1), bcr_d3acid3_l60_s647537_mpnn3 (rank 8) and bcr_d3acid3_l60_s647537_mpnn11 (rank 9).
- matters: A declared liability is understated three-fold, and Limitation 14 (line 1803, "**Rank 1 contacts a glycosylation sequon.**") names one design where three qualify. PK's 2026-10-03 reply made the glycan reassessment an explicit ask ("A sequon indicates potential occupancy... Include it when reassessing"). The two unnamed designs are the only members of the `d3acid3_l60_s647537` family, which §10/§11
- gate: gen_methods_submission.py: generate the §10b check table and Limitation 14 from finalist_footprints.json (counts AND names), the way LIMIT-FAMILY and LIMIT-AFFINITY are already generated.

### [overclaim] The §12 human-review and provenance attestations cover 17 named sequences, one of which is not in the submission, and two submitted rows are attested nowhere
- `submissions/01-egfr-METHODS.md`:1702
- verdict: CONFIRMED
- quoted: **Five sequences were added on 2026-10-05 and their review is recorded separately, because it is a human attestation and must not be inflated by restating a count.** They are `c5_cf_short__boltzgen_egfr_cropfree_short_48`, `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05`, `sd_d2c_101_l147_s144898_m_T65D`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA` and `cons_gap_h370_only__boltzgen_egfr_h370_018`.
- artifact: submissions/01-egfr.csv holds 18 rows. SIX were added on 2026-10-05 (their own assessment strings say "ADDED 2026-10-05" / "RESTORED 2026-10-05"): the four correctly named plus `bcr_d3acid3_l60_s647537_mpnn3` and `bcr_d3acid3_l60_s647537_mpnn11`. `sd_d2c_101_l147_s144898_m_T65D` is NOT in the CSV — §10 line 973-981 of this same document records that it was added and then removed, "replaced by `bcr
- matters: §12 is the one section that is explicitly a human attestation rather than a computation, and it attests review of a design that was withdrawn while leaving two shipped designs covered by no attestation at all. The same error runs through the Provenance paragraph (lines 1713-1715): "All seventeen sequences are de novo designs... This was checked by code over all seventeen, two ways" — the de-novo/n
- gate: check_claims RULE: the set of design names appearing in §12 must equal the CSV name set exactly (not a count — a set difference, reported by name in both directions).

### [overclaim] Slot accounting in §11 says five additions at ranks 1, 2, 3, 6, 8 and three slots unused; it is six additions at ranks 1, 2, 5, 7, 8, 9 and two slots unused
- `submissions/01-egfr-METHODS.md`:1135
- verdict: CONFIRMED
- quoted: They occupy ranks 1, 2, 3, 6 and 8. Nothing was displaced — the submission was at 12 of 20 and the five use free slots. Three slots remain unused.
- artifact: submissions/01-egfr.csv: 18 rows; the six 10-05 additions sit at ranks 1, 2, 5, 7, 8, 9 of the GENERATED rank table (§11.2). Rank 3 is `rimA01_r15_boltzgen_egfr_d3_rimA_20`, one of the original twelve. 20 − 18 = 2 slots unused, and §11.3 line 1286 of the same document states "The submission therefore stands at **18 of the 20 permitted**, with two slots deliberately unused rather than filled."
- matters: The opening paragraph of §11 — the first thing a grader reads about the submission's composition — gets the addition count, every added rank, and the remaining allocation wrong, and contradicts §11.3 155 lines later. PK's closing instruction on 2026-10-03 was "I would not fill the allocation simply to reach twenty"; the paragraph answering that ask is the one that cannot count the allocation. 12 +
- gate: gen_methods_submission.py: emit the whole "N designs / M added at ranks R / K slots unused" sentence as a generated block from the CSV plus the per-row ADDED/RESTORED markers.

### [overclaim] README declares three checks "Still open" that §10b and §4.4b in fact completed — and reports their results 14 lines earlier
- `README.md`:189
- verdict: CONFIRMED
- quoted: **Still open.** The partner-deletion free leg is a fixed-conformation diagnostic, not a measurement of the apo state. Full-ECD, glycan and receptor-state checks have not been applied to the finalist footprints. The rAC1 comparison needs structural contact recovery against 4UIP rather than predicted confidence attached to crystallographic coordinates. These are recorded in METHODS §13 rather than r
- artifact: All three were done. analysis/01-egfr/finalist_footprints.json carries n_outside_d3_crop, glycan_sequon_hits and hu_mo_identity_at_epitope for all 18 designs (METHODS §10b, lines 1007-1012, including the receptor-state verification at lines 1031-1041). analysis/01-egfr/rac1_contact_recovery.json + rac1_contact_recovery_chai.json are purely geometric contact recovery against 4UIP (METHODS §4.4b, li
- matters: This is the underclaim case the lens names, in the entry document a grader reads first — METHODS line 19 links the repository, and the README's "Still open" paragraph gives away three completed pieces of work, including the rAC1 contact recovery that is the strongest negative result in the submission and the receptor-state check that PK's 10-03 reply demanded. The apo free leg in the same sentence
- gate: Add a gate asserting that no item in README's "Still open" paragraph names an analysis artifact that exists under analysis/01-egfr/ — or generate that paragraph from a declared open-items list that §13 also reads.

### [overclaim] §1 presents a 5-row table as "a full PROPKA census of all 17 histidines" and concludes H433+H370 is the only bridgeable pair; the ectodomain coordinates contain four other pairs within 10.4 Å
- `submissions/01-egfr-METHODS.md`:91
- verdict: CONFIRMED
- quoted: A full PROPKA census of all 17 histidines in the EGFR ectodomain (6ARU, apo):
- artifact: targets/egfr/egfr_ecd_6aru.pdb chain A does hold exactly 17 histidines (mature 23, 121, 159, 209, 280, 334, 346, 359, 394, 409, 483, 535, 560, 566, 591, 594, 597 = canonical 47, 145, 183, 233, 304, 358, 370, 383, 418, 433, 507, 559, 584, 590, 615, 618, 621), but the table lists only the five inside the domain-III crop (targets/egfr/egfr_d3_6aru.pdb has exactly those 5). Minimum ring-nitrogen dista
- matters: Line 101 then states "H433 + H370 is the only pair close enough for one binder to bridge: a two-site ceiling of **43.1×**. Everything in §3 is an attempt to reach it" — an exclusivity claim over all 17 sites, asserted from a table showing 5, and false on the coordinates in this repo (one pair is closer than H370–H433). The whole two-site narrative in §3.4, §8.1 and §8.2 ("the cap is a property of 
- gate: Emit the census table from a script over targets/egfr/egfr_ecd_6aru.pdb (all 17 rows, or a row count that matches the header), and compute the "only pair within X Å" claim rather than asserting it.

### [overclaim] CONTROL-TABLE still calls the 10 no-KD molecules "measured negatives" and gitter-yolo10 a "measured non-binder" — the exact relabel PK required and the same file corrects elsewhere
- `outbox/CONTROL-TABLE.md`:29
- verdict: CONFIRMED
- quoted: 48–200 aa, from two independent groups. These are *measured* negatives, not presumed ones.
- artifact: The same file's §6 table heading (line 254) reads "**NO KD REPORTED (right-censored; affinity weaker than the quantifiable limit, or an expression/QC failure)**", and METHODS §4.4 line 364-373 records the correction: "Earlier versions of this section called them 'measured non-binders'. That is not what the data says... A molecule with no reported KD is one whose affinity is **right-censored**... W
- matters: CONTROL-TABLE is one of the two deliverables PK named as prerequisites ("Please send the scorer/fixtures and the new control table first"), and it carries the retracted label in two load-bearing places: line 29 introduces the whole panel as measured negatives, and line 363 — inside the bullet whose own parenthetical two lines up says "Corrected 2026-10-05: a missing KD is right-censoring, not a me
- gate: check_claims RULE over all DOCS: fail on /measured (non-?binder|negative)s?/ outside an explicitly quoted retraction block — the regex rules already cover phrasings of this kind for counts, and this is the same class.

### [adversarial-grader] Nine CSV rows quote the SUPERSEDED pH basis as the design's result, unlabelled, while their own graded column reads below the no-switch floor
- `submissions/01-egfr.csv`:18
- verdict: CONFIRMED
- quoted: pH 4.01x as the median over 6 refold poses; target's native H433.
- artifact: The same row's graded headline column ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE = 0.737, and ph_ratio_target_only_SUPERSEDED = 4.010 — i.e. the 4.01x quoted in the prose IS the superseded column. METHODS:1197 says rows below 1.20x 'carry no pH claim'. Same pattern, unlabelled, on CSV lines 4 (4.58x vs headline 4.256), 7 (4.57x vs 3.526), 11 (5.40x vs 1.835), 12 (5.43x vs 1.774), 13 (5.19x vs 4.8
- matters: The CSV is the graded upload and the assessment field is the only prose a grader who reads only the CSV sees. Four rows whose headline column is 0.737–1.062 (below the 1.20x PROPKA noise floor, i.e. no pH claim) describe themselves in the same cell as 4.01x–5.46x switches, one of them as '98% of the thermodynamic ceiling'. METHODS:1191 asserts the opposite — 'The superseded number ships as its own
- gate: Extend bin/check_claims.py to parse every float+'x' token in the CSV assessment field and require each to equal the same row's headline column within rounding, unless it is immediately preceded by a basis word ('superseded', 'target-only', 'all-site', 'partnered'). Better: generate the assessment's 

### [adversarial-grader] The L133E row ships three numbers that METHODS §11.6 formally withdraws, plus four self-contradictions against its own columns
- `submissions/01-egfr.csv`:15
- verdict: CONFIRMED
- quoted: Affinity HELD: human 0.594 -> 0.616, mouse 0.567 -> 0.435. WHY IT IS RANKED LAST: the pose spread is 4.24-11.64, i.e. 1.31x the median, against 0.04-0.90 for every other row. ... L133D at the same position reads 0.723x with human affinity 0.000 ... it needs 10-15 seeds and an isosteric L133Q control to be a result.
- artifact: METHODS:1368 withdraws 0.616 ('At n = 20 the human leg is 0.598 ... flat, not up') and the row's own ipsae_min_human column reads 0.5978. METHODS:1372 withdraws L133D 0.723x ('At n = 15 it is 2.221x'). METHODS:1379-1381 withdraws the spread ('reported as 1.31x its median on five poses. On twenty it is 4.38x, range 3.838 to 28.629') and the row's own ph_pose_spread_over_median column reads 4.38. 'a
- matters: §11.6 was explicitly rewritten on 2026-10-05 to retract these exact figures, and the retraction did not propagate to the graded artifact. A grader comparing the CSV against the methods finds the submission still shipping numbers its own methods document names as withdrawn — on the one row the document flags as its least reproducible. That is the difference between a declared limitation and an unco
- gate: Add a withdrawn-value blocklist to bin/check_claims.py: every number appearing in a METHODS 'what did not survive' / retraction block becomes a forbidden literal in submissions/01-egfr.csv and in README.md.

### [adversarial-grader] Both the CSV and METHODS assert the two swapped-in designs are "at most 0.467 identical to anything else submitted"; they are 0.867 identical to each other
- `submissions/01-egfr-METHODS.md`:1275
- verdict: CONFIRMED
- quoted: Both were removed and replaced by `bcr_d3acid3_l60_s647537_mpnn3` and `_mpnn11`, which are at most 0.467 identical to anything else submitted and open a backbone family that had no representation.
- artifact: Both sequences are exactly 60 aa and differ at 8 positions (1 V/I, 2 E/K, 4 E/K, 7 K/E, 15 K/E, 23 R/M, 40 R/E, 51 R/N) = 52/60 = 0.8667 positional identity; any gapped alignment can only raise this. The same false figure appears in the graded CSV on lines 9 and 10 ('Highest sequence identity to any other submitted design: <0.47.'), again at METHODS:979, and as the justification comment in bin/gen
- matters: This is the swap made specifically to satisfy PK's instruction, quoted in the document itself at METHODS:1268 — 'avoid filling available slots with nearly identical variants' — and the replacement pair is an 87%-identical sibling pair presented as the diverse alternative. Unlike the three >90% pairs, it is NOT declared a parent/mutant comparison, so by the document's own test at METHODS:1270 ('A n
- gate: Generate the identity claim rather than assert it: have bin/gen_methods_submission.py emit each design's max pairwise identity into a GENERATED block and into the CSV, and lower the check_claims.py near-duplicate disclosure threshold from 0.90 to ~0.80 so the sibling tier is reported, not just the p

### [adversarial-grader] The §12 human-review attestation names a sequence that is not in the submission and omits two that are; the provenance check covers 17 of 18 rows
- `submissions/01-egfr-METHODS.md`:1705
- verdict: CONFIRMED
- quoted: `sd_d2c_101_l147_s144898_m_T65D`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA` and
- artifact: submissions/01-egfr.csv holds 18 data rows and `sd_d2c_101_l147_s144898_m_T65D` is not one of them — METHODS:1275 records it as removed. Six rows carry a 2026-10-05 ADDED/RESTORED marker in the CSV, not five; the two the attestation omits are `bcr_d3acid3_l60_s647537_mpnn3` (CSV line 9) and `_mpnn11` (line 10). METHODS:1713-1715 then says 'All seventeen sequences are de novo designs ... checked by
- matters: Declarations are what the organisers ask for and the one section that is a human attestation rather than a computation. As written, the submitting researcher attests to having reviewed a sequence that was not uploaded, and attests to nothing for two rows that were. The substance of the provenance claim does survive independent recomputation (max identity to any known binder = 0.3967, rimA02 vs G53
- gate: Make the attestation list a GENERATED block keyed off the CSV: gen_methods_submission.py --check should fail if the set of names in §12 is not exactly the set of CSV rows, and if the stated count does not equal len(rows).

### [adversarial-grader] §10b's four reviewer-requested footprint checks were computed on a 17-design set that excludes two submitted designs and includes one that was withdrawn
- `submissions/01-egfr-METHODS.md`:1007
- verdict: CONFIRMED
- quoted: | check | result across all 17 |
- artifact: analysis/01-egfr/finalist_footprints.json contains exactly 17 entries. Set difference against the CSV: present in the artifact but not submitted = ['sd_d2c_101_l147_s144898_m_T65D']; submitted but absent from the artifact = ['bcr_d3acid3_l60_s647537_mpnn3', 'bcr_d3acid3_l60_s647537_mpnn11']. So '0 of 17 have any contact outside the 170 aa domain-III crop', '1 of 17 touches an N-glycosylation sequo
- matters: These are the four checks PK asked for by name on 2026-10-03 ('Compare the complete binder footprint in tethered and ligand-bound extended assemblies, including the second receptor, glycans and membrane-facing orientation'), and §10b opens by quoting that request. Presenting them as 'across all 17' invites a grader to read them as the submission's glycan, crop-adequacy and cross-species-epitope pr
- gate: Add a gate_sweep.py coverage check: every analysis JSON that METHODS cites as 'across all N' must have its key set equal the CSV name set, and the N in the prose must come from a GENERATED block rather than be typed.

### [adversarial-grader] README.md and HANDOFF.md — both named deliverables — describe a 10-design submission ranked on the retired metric
- `README.md`:13
- verdict: CONFIRMED
- quoted: | **[submissions/01-egfr.csv](submissions/01-egfr.csv)** | The submission: 10 designs, ranked. |
- artifact: submissions/01-egfr.csv holds 18 designs. README.md:34 repeats it ('We submitted **10 designs of the 20 allowed**'); README.md:17-20's headline claim reads 'took four independent ProteinMPNN sequences from 3.15-3.85x to 5.40-5.46x at no cost in predicted affinity', which are the ph_ratio_target_only_SUPERSEDED values (5.397, 5.428, 5.435, 5.461) for four rows whose graded headline column reads 1.0
- matters: README.md's 'Start here' table is the first thing a grader opens, and its four-line claim summary states the project's central causal result in a basis the methods document retired, with no label — the same defect as finding 1, at the front door. A reviewer who reads README then the CSV finds the design count, the headline ratios and the ranking all different, with no note saying the top of the fi
- gate: Have gen_methods_submission.py --write own the design count and the top-N table in README.md and HANDOFF.md as GENERATED blocks, and add a check_claims.py rule that no superseded-column value may appear in README.md or HANDOFF.md without the word 'superseded' in the same sentence.


## MED

### [pk-coverage] PREREGISTRATION records the instrument falsifier as surviving on the dual-species rescue PK ordered withdrawn
- `/Users/harish/code/adaptyv-2026/outbox/PREREGISTRATION.md`:402
- verdict: CONFIRMED
- quoted: * **The instrument fails** if ipSAE_min does not separate the 10 measured non-binders from the
  **one** measured binder ... **Status: fired on the human leg, survives on the
  dual-species leg — see §2.5, answered 2026-10-04 before outcomes.**
[§2.5:] Both molecules that beat it on human / are at **exactly 0.0000 on mouse, 5 of 5 dead seeds** — no interface, not a near miss. So the / separation r
- artifact: PK 2026-10-04, quoted verbatim in METHODS:402-405: 'a zero predicted interface is not an experimentally demonstrated specificity mechanism. Without matched mouse outcomes, this does not validate mouse binding or rescue the failed human control criterion.' METHODS:429-433 and CONTROL-TABLE:369-387 both withdraw the argument: 'That argument is withdrawn ... The human criterion failed and stays faile
- matters: This is not a stale number — it is the committed falsification verdict that the December analysis will be read against, and it records 'survives' on the one argument PK explicitly refused. The inbox status table for this email already says 'control recovery answered and its falsifier recorded as fired'. Left as is, the frozen plan tells a December reader that the instrument passed its own falsifie
- gate: Add an amendment note under §2.7 (the file's established pattern) stating the dual-species rescue is withdrawn and the falsifier is fired, and add a check_claims.py rule that any occurrence of 'survives on the dual-species' or 'mechanism rather than ... margin' must sit inside a withdrawal block.

### [pk-coverage] §10b concludes 'the crop is adequate' from footprints computed on crop-docked poses, contradicting limitation 15
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1010
- verdict: CONFIRMED
- quoted: | **full-ECD** | **0 of 17** have any contact outside the 170 aa domain-III crop (mature 311–480), so the crop is adequate and no footprint required the full ECD to assess |
- artifact: bin/finalist_footprints.py's own docstring describes this check as 'whether contacts fall OUTSIDE the 170 aa domain-III crop most of these binders were designed against -- i.e. whether the footprint is even assessable on the crop', and computes footprints 'from its own human-leg poses'. METHODS:1810 states '16 of 18 submitted designs are scored on the crop', so for 16 of 18 a contact outside matur
- matters: PK asked for the finalist footprints to be compared across assemblies, not for a self-consistency check; the result as written tells a grader the construct question is settled in the submission's favour while limitation 15 says the opposite on the only molecule with a solved complex. The '0 of 17' row can only ever read zero.
- gate: Have finalist_footprints.py refuse to emit the CROP verdict for designs whose poses were folded against the crop (report them as 'not assessable on this arm', the pattern §11.8 already uses for the relaxed-leg coverage table), and cross-reference limitation 15 in the generated block.

### [pk-coverage] PK's tethered-vs-extended footprint comparison was never done, and §10b presents 'the four checks' as the complete answer
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1014
- verdict: CONFIRMED
- quoted: **Domain II is not in play.** The earlier assessment concerned domain II; these binders do not
touch it. That resolves the question in the designs' favour but by irrelevance, not by passing.
[and 1031:] **Receptor state: the construct matches the assay, and I had this backwards for an hour.**
- artifact: PK 2026-10-03: 'Compare the complete binder footprint in tethered and ligand-bound extended assemblies, including the second receptor, glycans and membrane-facing orientation. Distance from one tether contact cannot settle this.' and 'A tethered-biased unliganded preparation is a reasonable working hypothesis; the proportions in Adaptyv's reagent are unknown.' finalist_footprints.json contains one
- matters: Three of the five things PK listed — extended/ligand-bound assembly, the second receptor, membrane-facing orientation — are absent from both the analysis and the limitations list, while the section reads as having discharged the request. Limitation 12 records only the Fab-templating caveat. A grader sees a four-row table marked 'result across all 17' with no indication that the comparison PK actua
- gate: List the unattempted checks explicitly in §10b and as a numbered limitation, and add a 'PK ask -> artifact' coverage manifest to gate_sweep.py that fails when an ask has no artifact file backing it.

### [pk-coverage] The glycan scan still misses the N-X-C sequon PK supplied, and §10b presents its eleven-sequon list as the ectodomain's full set
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1018
- verdict: PLAUSIBLE
- quoted: **Asn420**, one of eleven N-X-S/T sequons in the human ectodomain and one of four in domain III
(N328, N337, N389, N420).
- artifact: PK 2026-10-03: 'EGFR also has a documented atypical glycosylation site at mature N32, canonical N56, in an N-X-C motif, which your N-X-S/T scan misses. Include it when reassessing the domain-I fallback.' bin/finalist_footprints.py:56 implements only `re.finditer(r'N[^P][ST]', seq)` and prints 'N-glycosylation sequons in the human ECD: 11 (Asn at 104, 151, 172, 328, 337, 389, 420, 504, 544, 579, 59
- matters: Low direct consequence — all 18 designs sit in domain III and mature N32 is in domain I — but the document asserts a complete sequon count for the whole ectodomain that is known to be short by one, using the scan PK told us is incomplete, with no note of the gap. It also exposes the numbering problem: these positions are mature (mature 420 = canonical 444) while H433/H370 elsewhere in the document
- gate: Add the N-X-C motif to sequons() with a separate 'atypical' label, and have the printed line state the numbering convention; assert in --selftest that mature 32 is in the atypical set.

### [pk-coverage] The mature/canonical/PDB/mouse numbering table PK asked for twice does not exist, and the document mixes all three conventions
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:120
- verdict: CONFIRMED
- quoted: **ASP344 in PDB numbering, D368 canonical**, the convention H370 itself is quoted in — sits a
few Å from its ring. (The two numbering systems were mixed inside one sentence here ...)
- artifact: PK 2026-09-29, correction 1: 'Put mature, canonical, PDB-chain and aligned mouse positions in one table before generating constraints,' repeated in his closing ask: 'please send the numbering/alignment table'. He supplied the mappings (H409 = canonical H433; I467/S468 = canonical I491/S492; the seven acidic positions as D347, D368, D379, D388, D416, E391, E455; mouse Q01279 alignment). Grep across
- matters: This was the first of five pre-ranking corrections and remains unbuilt after three emails. Its absence is live in the deliverables: §1 and §3 use canonical (H433, H370), §2 mixes PDB and canonical in one sentence and says so, §4.4b and §10b use mature (epitope 411-489, footprint 316-474, sequons N328-N420), so the shared-epitope list at §10b:1055 contains '409' — which is mature numbering for the 
- gate: Generate one numbering table (mature / canonical P00533 / 6ARU chain A / aligned mouse Q01279) from targets/egfr/*.faa by alignment, publish it as a §1 GENERATED block, and have check_claims.py require a convention tag on any residue reference matching /[HDEN]\d{2,3}/.

### [pk-coverage] CONTROL-TABLE describes the affinity_above_null column as shipped, the exclusion count as 38, and the family split as unfixed — all three superseded
- `/Users/harish/code/adaptyv-2026/outbox/CONTROL-TABLE.md`:80
- verdict: CONFIRMED
- quoted: We have left the column in the submission rather than dropping it mid-flight, and said the same
thing in METHODS §4.5 and §11 so all three documents now agree.
[line 222:] §10 of the methods document lists 38 molecules rejected on it.
[line 394:] **One thing I have not fixed.** You asked twice for the historical EGFR data to be split by design family.
- artifact: submissions/01-egfr.csv has 15 columns and no affinity_above_null; METHODS:1229-1233 states 'The `affinity_above_null` column has therefore been REMOVED from the emitted CSV (2026-10-05)'. METHODS §10 now reports 15 distinct List-A sequences and a 75-molecule union, and §10 itself says the '38 distinct molecules' figure was an artefact of name-keyed double counting. The family split PK asked for I
- matters: This file is one of the two artifacts PK named as prerequisites and it is still unsent. In its current state it tells him a column is in the graded CSV that is not, quotes a molecule count the methods document retracted, and says a request is unfixed on the same page where it is answered — three ways to make the reply look careless about exactly the bookkeeping he has been pressing on.
- gate: Bring outbox/*.md inside check_claims.py's inventory check (CSV column names must match the emitted header; retired literals 38, 0.2218-as-shipped, 0.778 banned), and delete or date-stamp §5/§6's superseded paragraphs.

### [pk-coverage] §12 reports published-structure coverage as 10 of 12 when the submission is 18 and eight designs have no structure
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1731
- verdict: CONFIRMED
- quoted: **Coverage is 10 of 12**: the two designs added latest, `bc_s360518_mpnn9_A22D` and
`rimA01_r15_L133E`, have no published structure yet.
- artifact: submissions/structures/ contains 10 .cif files plus a README. The CSV has 18 rows, so eight designs have no published structure: bc_s360518_mpnn9_A22D, rimA01_r15_L133E, c5_cf_short__boltzgen_egfr_cropfree_short_48, c5_cr_crop_patch__boltzgen_egfr_crop_patch_05, ss_bc_s831683_mpnn6_S15D_S62H_routeA, cons_gap_h370_only__boltzgen_egfr_h370_018, bcr_d3acid3_l60_s647537_mpnn3 and bcr_d3acid3_l60_s6475
- matters: The sentence exists precisely because an earlier version overstated coverage ('an earlier version of this line claimed full coverage of "all ten designs" when the submission held twelve'), and it now understates the denominator by six. PK's 10-03 ask was a finalist table with model provenance; a grader who wants to inspect the six highest-ranked designs' poses finds that five of the top six have n
- gate: Make the coverage line a GENERATED block computed as len(glob('submissions/structures/*.cif')) over the CSV row count, naming the missing designs; fail gen_methods_submission.py --check on drift.

### [pk-coverage] HANDOFF.md's verified-state section describes a 10-design submission on the superseded ranking basis
- `/Users/harish/code/adaptyv-2026/HANDOFF.md`:48
- verdict: CONFIRMED
- quoted: **Submission: `submissions/01-egfr.csv`, 10 designs** (Track 3 allows 20; we gave half back
on purpose — see METHODS §11). `bin/check_discards.py` PASSES, exit 0.

     1  5.461  bc_s831683_mpnn8_S15D         hu 0.7760  mo 0.7637   65aa
[and line 25:] The cap is 20 designs per account per collection. We ship 10.
- artifact: The CSV holds 18 designs. bc_s831683_mpnn8_S15D is rank 18, not rank 1 — METHODS:1199-1200: 'bc_s831683_mpnn8_S15D led this submission at 5.461× before the correction and is now rank 11 at 1.023×' (rank 18 after the additions). Every number in the HANDOFF ladder is the target-only basis that §11.1 supersedes, under a heading reading 'STATE — everything below is verified, not remembered'.
- matters: HANDOFF is in the graded artifact set and is the document a reader picks up to learn the submission's current state. It asserts verification over a ten-row ladder whose top design is now last and whose ratios come from the basis the project renamed SUPERSEDED. It also still records 'All 10 of our designs switch on the target's native H433', where the generated block in METHODS says 17 of 18 switch
- gate: Generate HANDOFF's state ladder from the emitted CSV via gen_methods_submission.py, or mark the section with an as-of timestamp and point to §11.2 as the live table; add HANDOFF.md to the --check inventory.

### [pk-coverage] Ranks 1 and 2 are the two weakest predicted interfaces in the submission, ordered there by the pH ratio alone
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1233
- verdict: PLAUSIBLE
- quoted: **Assessable designs rank ahead of unassessable ones within tier 1.** ... We neither demote it on the pH axis nor score it at 0.0000
(the §4.2 error); we place it after the designs where both axes mean something.
... It is a judgement that credible-interface-first is the more
defensible frame, following the reviewer instruction to apply eligibility and interface checks
before the challenge priorit
- artifact: The CSV's rank-1 and rank-2 rows read ipsae_min_human 0.2415 / 0.1812 and 0.1324 / 0.2041 — the two weakest assessable interfaces of the sixteen. bin/emit_submission_csv.py rank_key() orders tier 1 by the pH ratio with only a format penalty and a SPREAD_BAR penalty; no interface quality enters the ordering, and MIN_AFFINITY is None. METHODS §10:966 flagged the eventual rank 1 on arrival: '`c5_cf_s
- matters: 'Credible-interface-first' is the stated ordering principle and the stated reason the two VHH rows sit at 12-13, but it is applied only to format-unassessable rows. The designs with the weakest readable interfaces lead the file on apparent selectivity, which is the specific failure mode PK named. Limitation 3 admits the objective does not discriminate binders in general; it does not say the orderi
- gate: Either state in §11.2 that tier-1 order is the pH ratio alone and that ranks 1-2 carry the set's weakest interfaces, or add an explicit interface floor to rank_key() and record it in the selection history; add a --selftest assertion that the top-ranked row is not the minimum of ipsae_min_human acros

### [pk-coverage] The seed-instability pilot PK specified (20-30 candidates, ~10 seeds each) was never run and is not recorded as outstanding
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:339
- verdict: CONFIRMED
- quoted: We therefore require **n ≥ 5 poses** before a ratio may rank a design, and we
report the per-pose spread
- artifact: PK 2026-10-03: 'To study seed instability, a practical pilot is 20–30 diverse candidates, enriched near decision boundaries, with approximately ten seeds each; assess rank changes, pose consistency and threshold crossings. This is a diagnostic starting point, not a powered validation sample.' The deepest sampling in the repository is the three-variant L133 triad at 15 seeds (runs/esmfold2/w3_triad
- matters: §11.6 is the evidence that this pilot was the right diagnostic: quadrupling the sampling on one design took the measured spread from 1.31x to 4.38x and moved it out of tier 1, and §6 reports that 5 of 12 single-pose switches fell below threshold on re-measurement. Fourteen of eighteen shipped rows still sit at 5-11 poses. The document reports the n>=5 floor as the answer to instability without not
- gate: Add it as a numbered limitation naming the design and seed counts actually achieved, and add a 'PK ask -> artifact or declared-open' manifest to gate_sweep.py so an unmet methodological ask cannot be silently absent from §13.

### [pk-coverage] §4.4b attributes the rAC1 control poses to ESMFold2-Fast; the folding path defaults to the Full model
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:523
- verdict: CONFIRMED
- quoted: **ESMFold2 is not the weak link here.** An architecturally independent model ... reproduces its
crop failure to within rounding on every column. Whatever is wrong is not specific to
ESMFold2-Fast, which removes the most convenient explanation for the zero-scoring positives.
- artifact: biomodals/modal_esmfold2.py:41 sets ESMFOLD2_HF_REPO default to 'biohub/ESMFold2' (its own comment: 'Full (default)'; 'Fast: ESMFOLD2_HF_REPO=biohub/ESMFold2-Fast'). bin/score-esmfold2.sh — the script that produced the design and control folds, including runs/esmfold2/w1_rac1 — sets no repo override, so it runs Full. The only Fast arms on disk are runs/esmfold2-fast/ipsae-probe-fast and runs/gate/
- matters: The sentence's whole load is which model arm produced the failure, in the project's only crystallographic control. Naming the Fast arm suggests the primary scoring ran on Fast — the configuration PK demoted to a robustness check — when the pipeline in fact followed his instruction. The Fast and Boltz robustness arms exist on disk (runs/gate/fast, runs/gate/boltz) and are reported nowhere, so a rea
- gate: Record the model repo and revision per run directory and have gen_methods_submission.py emit the predictor identity as a generated string, so no prose can name an arm the run did not use.

### [pk-coverage] Mechanism A is asserted twice to have been ruled out, with no mechanism-A result reported anywhere
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:1188
- verdict: CONFIRMED
- quoted: the same mechanism as
the 0.702× steric floor of §6, and the same physics that defeated mechanism A (§3.5).
[and 1448:] a desolvation shift with no electrostatic partner, the same artefact class
this project used to rule out Mechanism A
- artifact: PK 2026-09-29, answer 3: 'Give A most of the initial design effort, conditional on finding a suitable local acidic surface; retain B as a smaller exploratory branch.' METHODS §1 instead reads 'Mechanism B ... We used this', and §3.5 — the section cited as the evidence that mechanism A was defeated — contains no mechanism-A analysis at all; it is the retraction of the binding/switching trade-off. T
- matters: PK's primary mechanism recommendation was inverted, and the two sentences that justify the inversion point at a section that does not contain the evidence — while 540 scored mechanism-A complexes sit unreported in the analysis directory. Either the arm's result belongs in the document (it would be the direct answer to his answer 3) or the 'defeated / ruled out' claims need the data behind them.
- gate: Report the mechanism-A arm (n, ratio distribution, best designs) as a short §3 subsection sourced from ph_gate_mechA_all.json, and extend bin/check_references.py from 'section exists' to 'cited section contains the named quantity'.

### [pk-coverage] The burial atom count is still the sole support for the H370 conclusion; no solvent-accessible surface area was computed
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:206
- verdict: CONFIRMED
- quoted: H370
carries 153 heavy atoms within 10 Å of its ring against H433's 57, so reaching it costs the
interface area binding needs.
- artifact: PK 2026-10-03: 'The burial atom count is a useful proxy, not a substitute for solvent-accessible surface area.' The strings 'SASA' and 'solvent-accessible' appear in no deliverable; biomodals/modal_sasa.py exists and was never run against the EGFR targets (no SASA artifact in analysis/01-egfr/).
- matters: The 153-vs-57 count is load-bearing in three places — §3.4's conclusion that the two-site route cost interface area, §3.5's surviving site-specific trade-off, and §11.1's desolvation argument — and it is the only quantity offered for 'costs the interface area binding needs'. PK flagged the proxy explicitly and the substitute he named is one unrun script away. Neither the limitation nor the proxy's
- gate: Compute per-histidine SASA with modal_sasa.py (or gemmi/freesasa locally) on 6ARU and report it beside the atom count; failing that, label the count a proxy in §3.4/§3.5 and add it to §13.

### [arithmetic] All three quoted pKa-perturbation rank spans disagree with ph_pka_perturbation.json, and "1–17" is called "the entire submission" of 18
- `submissions/01-egfr-METHODS.md`:1489
- verdict: CONFIRMED
- quoted: `d2c_mpnn13_S88D_serasp`
spans ranks 1–17 — the entire submission. `bc_s360518_mpnn9_A22D` spans 2–15 and
`bc_s831683_mpnn6_S15D` spans 4–17.
- artifact: ph_pka_perturbation.json p5_rank/p95_rank: d2c_mpnn13_S88D_serasp 1–18; bc_s360518_mpnn9_A22D 1–16; bc_s831683_mpnn6_S15D 5–18. Same section's "bc_s831683_mpnn9_WT never leaves ranks 13–17" (line 1495) is 14–18 in the artifact.
- matters: Four separate central-90% intervals are quoted and none of them matches the file they are drawn from; all four are stated as a 17-design table while the submission is 18. Each is individually plausible, so regex matching against a value list cannot catch it — but a reviewer recomputing the interval gets a different number every time, which puts the whole §11.7 table in question rather than just on
- gate: Generate the whole §11.7 perturbation paragraph (table rows, the named spans, the top-set membership) from ph_pka_perturbation.json inside a GENERATED block, the way §11.1–11.3 already are.

### [arithmetic] The §11 opening paragraph's counts describe a 17-design submission under an "18 designs" heading, and the five additions' ranks are wrong
- `submissions/01-egfr-METHODS.md`:1135
- verdict: CONFIRMED
- quoted: They occupy ranks 1, 2, 3, 6 and 8. Nothing was displaced — the
submission was at 12 of 20 and the five use free slots. Three slots remain unused.
- artifact: The CSV has 18 rows; 12 + 5 = 17, so the restored `ss_bc_s831683_mpnn6_S15D_S62H_routeA` is unaccounted for and only two slots remain unused — as §11.3 line 1283 states: "The submission therefore stands at **18 of the 20 permitted**, with two slots deliberately unused". The five designs whose assessment says "ADDED 2026-10-05 from the reopened exclusion pool" occupy ranks 1, 2, 7, 8 and 9 in this 
- matters: This is the first paragraph of the section that describes the graded upload, and it both miscounts the submission (17 vs 18, three unused slots vs two) and misplaces the five additions against the generated table 70 lines below it. Track 3 allows 20 and the number of used slots is a graded fact; two statements of it in one section disagree.
- gate: Emit the whole opening paragraph's arithmetic (total designs, designs added, ranks occupied, slots unused) as a GENERATED block keyed on the CSV row count and the "ADDED 2026-10-05" marker in the assessment field.

### [arithmetic] A graded CSV assessment says six of twelve rows share the backbone; the methods' generated family list says seven of eighteen
- `submissions/01-egfr.csv`:18
- verdict: CONFIRMED
- quoted: Six of the twelve submitted rows share this backbone: that is a deliberate trade of panel diversity for wet-lab replication of the only causal result we have.
- artifact: Seven of the 18 shipped sequences are on `d3acid_l65_s831683`: ss_..._S62H_routeA, mpnn6_S15D, mpnn19_S15D, mpnn9_S15D, mpnn9_WT, bc_d3acid_l65_s831683_mpnn11, mpnn8_S15D. METHODS' GENERATED FAMILY-LIST (line 1290) reads "`d3acid_l65_s831683` **x7** (ranks 5, 10, 11, 15, 16, 17, 18)" and generated limitation 19 reads "7 of the 18 submitted designs sit on one backbone".
- matters: The clustering of the submission on one backbone is the effective-n claim and the single biggest caveat on the only causal result. The graded file understates it on both numerator and denominator while the methods document states it correctly twice from generated blocks — so the CSV reader and the METHODS reader get different effective-n.
- gate: Generate the family-share sentence of each assessment string from the same family partition that feeds the FAMILY-LIST block, rather than writing it by hand.

### [arithmetic] §11.7 states three mutually incompatible human-leg pose totals (75, 158, 165) for the same analysis
- `submissions/01-egfr-METHODS.md`:1425
- verdict: CONFIRMED
- quoted: recomputes
every submitted design on the same 75 human-leg poses through one code path ... all **17** are on
one footing; across the 17 the analysis covers **158 human-leg poses**
- artifact: ph_sensitivity.json n_poses sums to 165 over the 18 submitted designs (and 160 over the 17 excluding routeA); no subset gives 158, and no subset of the submission gives 75. §11.7's own later sentence (line 1517) uses the right figure: "Across all **165 poses of the eighteen submitted designs**". §11.1 repeats the stale pair twice in six lines: "Measured over 75 poses, n = 5–11 per design" (line 11
- matters: The pose count is the sample size of the whole pH analysis, and the document gives it as 75, 76, 158 and 165 in four places. The "n = 5–11 per design" range is contradicted by the generated table immediately below it (20, 20, 26). None of this is a rounding question; one of the four is right and three are leftovers from earlier submission states.
- gate: Emit every pose total and per-design n range as generated blocks summed from ph_sensitivity.json; add a check_claims rule that any "N human-leg poses" literal equals the artifact sum.

### [arithmetic] §10b's shared-epitope footprint numbers do not match finalist_footprints.json: union 58 vs 62, and the "20 residues ≥80%" list contains a residue at 77.8%
- `submissions/01-egfr-METHODS.md`:1054
- verdict: CONFIRMED
- quoted: Across all 17 designs the union of contacted residues is only **58 distinct positions
(mature 316–474)**, and **20 residues are contacted by at least 80% of the designs**: 323, 325,
348, 349, 350, 353, 355, 357, 382, 384, 408, 409, 411, 412, 417, 418, 438, 440, 465, 467.
- artifact: Recomputed from finalist_footprints.json over all 18 designs: the union is 62 distinct positions (range 316–474 is right), and 19 residues are contacted by ≥80% — the printed list minus 323, which is contacted by 14 of 18 designs (77.8%). The same figures also appear as "All 17 designs ... 20 residues are shared by ≥80% of them" in limitation 11 (line 1791) and in README line 170.
- matters: This is the correlated-failure finding the section was added to surface, and both of its quantities are off: the epitope union is 7% larger than stated and the shared-core list includes a member below its own stated threshold, so the list as printed is 20 items against a 19-item fact. The direction overstates the concentration slightly, but the mismatch means the number cannot be reproduced from t
- gate: Make the union size, the mature range, and the ≥80% residue list a GENERATED block computed from finalist_footprints.json, with the threshold applied in code rather than in prose.

### [arithmetic] §11.7's apo arm claims it folded "every submitted binder" and that d2c rises "from fifth"; one design has no apo structure and d2c is seventh
- `submissions/01-egfr-METHODS.md`:1568
- verdict: CONFIRMED
- quoted: The largest single mover is `d2c_mpnn13_S88D_serasp`, **3.526 → 5.041 (+43%)**, which rises
from fifth to third.
- artifact: ph_apo_freeleg.json: on the deletion basis d2c_mpnn13_S88D_serasp is 7th of the 17 designs with apo structures (5.6559, 5.5459, 4.838, 4.8125, 4.2561, 3.7377, then 3.5262) and 8th of 18 in the CSV; it rises to 3rd on apo, so the move is seventh-to-third. Its +43% is also not the largest fold move — bc_s831683_mpnn8_S15D moves 2.018×, mpnn9_WT 1.853×, mpnn9_S15D 1.756×. And line 1526's "every submi
- matters: The apo arm is the half of PK's free-leg request that §11.7 presents as completed, and its coverage statement ("every submitted binder") is false for one shipped design while its internal check (8 + 9 = 17) silently drops that design. The rank-movement claim is two places off. Together they make the arm look more complete and its largest effect better placed than the artifact supports.
- gate: Generate the apo coverage line (n folded, n with/without structure, named gaps) and the rank-movement sentence from ph_apo_freeleg.json, mirroring the §11.8 coverage table that already does this correctly.

### [arithmetic] Limitation 24 says seven designs then eight designs for the same set in one sentence; nine actually carry no binder histidine
- `submissions/01-egfr-METHODS.md`:1843
- verdict: CONFIRMED
- quoted: Seven designs carry
    no binder histidine — their switch is target-borne — and a binder-only relaxation cannot
    move them, so the eight designs whose signal sits on EGFR's own H370/H433 are exactly the
    ones this sensitivity check cannot cover.
- artifact: ph_relaxed_freeleg.json summary: n_compared 10, n_no_movable_site 7, n_no_relaxed_structure 1. Counting histidines in the shipped sequences, 9 of 18 designs carry none — the 7 "no movable site" plus rimA01_r15_L133E and bcr_d3acid3_l60_s647537_mpnn3, both of which are covered via a relaxed target chain. The 8th uncovered design is ss_..._S62H_routeA, which carries 4 binder histidines, so it is not
- matters: A limitation whose job is to state the reach of a sensitivity arm gives two different counts for the same set one clause apart, and the larger count mislabels a design with four binder histidines as target-borne. §11.8's own coverage table (line 1596) has the partition right, so the limitation contradicts the section it summarises.
- gate: Generate limitation 24's counts from the ph_relaxed_freeleg.json summary block (n_compared / n_no_movable_site / n_no_relaxed_structure) the way limitations 19 and 22 are already generated.

### [arithmetic] §4.4b's refold accounting does not sum: "one design already uses it, the other seventeen" against "16 of the 18 scored on the crop"
- `submissions/01-egfr-METHODS.md`:541
- verdict: PLAUSIBLE
- quoted: We are not able to rescore the submission on full ECD
before the deadline — one design already uses it, the other seventeen would need refolding and
re-gating
- artifact: The same section's heading sentence (line 537) and limitation 15 (line 1810) both say "16 of 18 submitted designs are scored on the crop" / "16 of the 18 submitted designs were scored against that crop", which leaves 2 on the full ECD, not 1 — so the remainder is sixteen, not seventeen. ph_sensitivity.json corroborates two full-ECD designs (d2c_mpnn13_S88D_serasp and sd_d2c..._T65D carry 123 titra
- matters: The construct choice is the strongest negative finding in the document (neither predictor docks the one solved complex correctly on the crop), so how much of the submission is exposed to it is the number that matters. Two statements of it four lines apart disagree, and the one in the mitigation sentence is the one that understates the work outstanding.
- gate: Derive the crop-vs-ECD split once from the pose index / n_sites in ph_sensitivity.json and emit both the "16 of 18" and the "would need refolding" counts from it.

### [arithmetic] §4.4 calls a no-KD molecule "already measured not to bind" in the same section that corrects exactly that label
- `submissions/01-egfr-METHODS.md`:441
- verdict: CONFIRMED
- quoted: A human-leg-only pipeline would have submitted a
molecule already measured not to bind.
- artifact: The same section, lines 364–371: "Earlier versions of this section called them 'measured non-binders'. That is not what the data says ... A molecule with no reported KD is one whose affinity is **right-censored** ... It is not a measurement of zero affinity ... We do not know the true affinity of any of the ten." PK's 2026-10-03 reply: "Distinguish assay failure and missing KD from no binding dete
- matters: gitter-yolo10 is a right-censored observation; the section says so three times and then states the retracted reading as fact in its punchline, which is the sentence most likely to be quoted. It is the exact labelling error PK corrected by name, surviving in the paragraph that was rewritten to fix it, and it converts a bound-from-above into a measured negative.
- gate: Add a check_claims regex forbidding "measured not to bind", "measured non-binder" and "non-binder" within §4.4/§4.5 and the CSV, with the approved phrasing "no KD reported" / "right-censored" as the only substitutes.

### [arithmetic] §9 attributes ranks 1–4 to the BindCraft pool; ranks 1–3 are BoltzGen designs
- `submissions/01-egfr-METHODS.md`:801
- verdict: CONFIRMED
- quoted: **114 of 238 clear Level ≥ 3. The old single bar passed 0 of 238** — including the entire
BindCraft pool, which supplies ranks 1–4, 8 and 9 of this submission.
- artifact: Ranks 1, 2 and 3 in this document's GENERATED rank table are c5_cf_short__boltzgen_egfr_cropfree_short_48, c5_cr_crop_patch__boltzgen_egfr_crop_patch_05 and rimA01_r15_boltzgen_egfr_d3_rimA_20 — all BoltzGen. BindCraft (bc_/bcr_) supplies ranks 4, 5, 8, 9, 10, 11, 15, 16, 17 and 18.
- matters: The claim is doing real work: it says the novelty-implementation bug would have voided the top of the submission. The top three rows are from the other generator, so the stated consequence attaches to the wrong designs, and the ten rows BindCraft actually supplies are undercounted. The same "ranks 1–4, 8 and 9" pattern looks like the pre-restoration addition list rather than a generator partition.
- gate: Generate generator-provenance statements from the CSV name prefixes (boltzgen / bc_ / bcr_ / vhh) rather than writing rank lists by hand.

### [arithmetic] README carries four counts the methods document supersedes, including one METHODS lists as a known miscount
- `README.md`:182
- verdict: CONFIRMED
- quoted: the count of
  designs carrying binder histidines is ten of seventeen, not the eight that was carried forward
  without recounting.
- artifact: Counting histidines in the 18 shipped sequences gives 9 of 18, which is what METHODS' GENERATED block states (line 1144); METHODS line 1147 explicitly lists "ten of seventeen" as one of the five hand counts that were wrong. README also still says "the minimum of the three for all twelve designs" (line 99, now 18 of 18 — bin/emit_submission_csv.py prints n_env of len(scored)), "12 → 17 of the 20 pe
- matters: README is a graded artifact and is the first file a reviewer opens. It states as a correction the very figure METHODS names as an error, and its submission size, envelope count and near-duplicate count all describe the pre-restoration 17-design state. A reviewer comparing README to METHODS finds them disagreeing about how many designs were submitted.
- gate: Make the README "what changed" bullets reference generated METHODS blocks, or extend gen_methods_submission.py --check to cover README's count literals (design total, envelope n of n, binder-histidine count, near-duplicate count) against the emitted CSV.

### [unsupported] The graded CSV's two-coupled-sites claim rests on a premise the artifact refutes: H370 is not inert in "every other design" — rank 1 reads it at 5.819x
- `submissions/01-egfr.csv`:15
- verdict: CONFIRMED
- quoted: It is the only design here with evidence of TWO coupled sites: all-site median 5.659x, above the 5.55x thermodynamic ceiling for H433 alone, with H433 steady at 4.458-4.675x across all five poses while H370 -- which reads 0.92-0.94x in every other design in this submission -- pulls 0.921 to 2.584x.
- artifact: analysis/01-egfr/ph_sensitivity.json, per-design median of target:H370 ratio over that design's own human-leg poses: c5_cf_short__boltzgen_egfr_cropfree_short_48 (rank 1) = 5.819 (max 5.945, n=6). The other 16 submitted designs span 0.846-0.954, not 0.92-0.94 (bc_s360518_mpnn9_A22D 0.846, ss_bc_s831683_mpnn6_S15D_S62H_routeA 0.880, bc_d3acid_l65_s831683_mpnn11 0.868). METHODS.md line 67, a GENERAT
- matters: The uniqueness claim ("the only design here with evidence of TWO coupled sites") is argued by contrast: H370 is supposedly inert everywhere else, so L133E moving it must be a designed second site. The submission's own rank-1 design moves H370 by 5.819x — four times harder than L133E's 1.308x median — and §1's generated block names it. The contrast that licenses the mechanism claim does not exist, 
- gate: In gen_methods_submission.py / emit_submission_csv.py: for every submitted row, recompute each target histidine's median ratio from ph_sensitivity.json and fail if an assessment string asserts a range for a site that any other submitted row falls outside. Cheaper partial gate: forbid the literal str

### [unsupported] The §12 human-review and provenance attestations name a design that was not submitted and omit two that were
- `submissions/01-egfr-METHODS.md`:1703
- verdict: CONFIRMED
- quoted: They are
`c5_cf_short__boltzgen_egfr_cropfree_short_48`, `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05`,
`sd_d2c_101_l147_s144898_m_T65D`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA` and
`cons_gap_h370_only__boltzgen_egfr_h370_018`. [...] **Provenance.** All seventeen sequences are de novo designs from this project's own generation
runs [...] This was checked by code over all seventeen
- artifact: submissions/01-egfr.csv has 18 rows. sd_d2c_101_l147_s144898_m_T65D is not one of them; bcr_d3acid3_l60_s647537_mpnn3 and bcr_d3acid3_l60_s647537_mpnn11 are (rows 9 and 10). §10 line 973-981 and §11.3 line 1270-1276 both record that swap. So 12 attested + 5 named = 17 attestations against 18 graded rows, with one attestation pointing at a design that does not ship.
- matters: This is the one claim in the document that is a human attestation rather than a computation, and §12 itself says it "must not be inflated by restating a count". As written, two graded sequences carry no human review and no provenance/de-novo check, and the attestation covers a sequence a grader cannot find in the file. Related stale counts downstream: §11.3 line 1292 "all 17 are first-time submiss
- gate: Have gen_methods_submission.py emit the §12 added-design list and the provenance count as GENERATED blocks keyed on the CSV, so a swap cannot leave the attestation behind; and add the name-set equality check above so any artifact whose design set differs from the CSV fails the sweep.

### [unsupported] L133E's "second site" is a +0.54 pKa-unit median shift — inside the ±0.8 PROPKA error §11.7 adopts as noise, sign-reversed in 2 of 20 poses — and §11.7 says so while §11.6 and the CSV say the opposite
- `submissions/01-egfr-METHODS.md`:1360
- verdict: CONFIRMED
- quoted: **It is still a two-site reading.** Over the 15 triad poses, **9 of 15 exceed the 5.55×
single-site ceiling** for H433 alone (§1) and **3 of 15 exceed the 7.94× one-proton bound**.
With two sites moving the ceiling is 7.94² = 63×, so these are not ceiling violations — they are
values one site cannot produce.
- artifact: analysis/01-egfr/ph_sensitivity.json, rimA01_r15_L133E, target:H370 across all 20 poses: pKa_bound − pKa_free = −0.40, −0.39, +0.22, +0.25, +0.26, +0.28, +0.40, +0.46, +0.51, +0.52, +0.56, +0.69, +0.70, +0.76, +0.78, +0.88, +1.22, +1.25, +1.26, +2.56. Median +0.54 pKa units. PROPKA's own RMSD, which §11.7 line 1481 adopts as σ, is ~0.8. Six of the 20 H370 ratios fall below the 1.20× bar §6 line 66
- matters: "Values one site cannot produce" is only true if the gate's per-site factors are calibrated linkage ratios; §11.7's own ASP33 example shows they are not, and §1 line 81-88 already argues that designs piling up on an analytic extreme signal "protonation-model saturation [...] rather than evidence that those designs genuinely occupy a physical extreme". The document applies that scepticism to the lo
- gate: In ph_gate_multisite.py / ph_sensitivity_multisite.py, emit per-site ΔpKa alongside the ratio and flag any site whose |ΔpKa| is below PROPKA's RMSD, or whose sign flips across poses, as NOT-RESOLVED; then a check_claims rule that fails any "two-site"/"second site"/"coupled sites" assertion naming a 

### [contradiction] HANDOFF's "verified, not remembered" state block is a 10-design submission with 6 stale affinity values
- `HANDOFF.md`:48
- verdict: CONFIRMED
- quoted: **Submission: `submissions/01-egfr.csv`, 10 designs** (Track 3 allows 20; we gave half back
on purpose — see METHODS §11). `bin/check_discards.py` PASSES, exit 0.
- artifact: The CSV holds 18 designs. The rank list at HANDOFF lines 51-60 also disagrees with the CSV on 6 of 10 affinity pairs: rimA02_d3_rimA_14_vhh "hu 0.2275 mo 0.4511" (CSV 0.2186/0.4468), h370_020_vhh "hu 0.4409" (CSV 0.4171), bc_s831683_mpnn9_S15D "hu 0.8025 mo 0.8026" (CSV 0.8040/0.8037), rimA01_r15 "hu 0.5765 mo 0.5603" (CSV 0.5938/0.5668), bc_d3acid_l65_s831683_mpnn11 "hu 0.7949" (CSV 0.7963), bc_s
- matters: The section is headed "## 1. STATE — everything below is verified, not remembered" and the file carries a 2026-10-05 edit (line 4), so it reads as current. It states the submission size, ranking and affinity numbers all differently from the file that is actually uploaded.
- gate: Replace the hand-typed HANDOFF state block with a generated block read from submissions/01-egfr.csv, or add a gate that fails if HANDOFF's design count or any numeric in its rank list does not match the CSV.

### [contradiction] The "81% mature human EGF" figure METHODS §4.1 corrects to 32% is still asserted in README and CONTROL-TABLE
- `README.md`:200
- verdict: CONFIRMED
- quoted: - Our negative control turned out to be **81% mature human EGF** — the agonist. Every bar in
  the project had been calibrated against it. The thresholds are retired, not replaced.
- artifact: METHODS §4.1 line 298-301: "Stated precisely, because an earlier version said \"81% mature human EGF\" and that is the percentage **of EGF**, not of the molecule: `NEG_nonbinder` is **134 residues**, so 43 are EGF-derived and **91 are not** — 32% of the molecule, not 81%. \"Two insertions in the loops\" also undersold 81 extra residues." outbox/CONTROL-TABLE.md line 38 carries both retracted forms
- matters: README presents this as one of the two errors "worth knowing before you read any number here", using the wrong denominator. CONTROL-TABLE, the artifact going to the external reviewer, repeats it plus the "two insertions" framing the same correction calls an undercount of 81 residues.
- gate: Forbid "81% EGF", "81% mature human EGF" and "two insertions" in any file under submissions/ or outbox/ plus README.md and HANDOFF.md, via a check_claims.py regex rule.

### [contradiction] Two generated blocks in the same subsection give d2c's binder-histidine count as 14 and as 2
- `submissions/01-egfr-METHODS.md`:1144
- verdict: CONFIRMED
- quoted: **nine of the eighteen submitted designs carry at least one histidine of their own**: `d2c_mpnn13_S88D_serasp` (14); `ss_bc_s831683_mpnn6_S15D_S62H_routeA` (4); `bc_d3acid_l65_s831683_mpnn11`, `bc_s831683_mpnn19_S15D`, `bc_s831683_mpnn6_S15D`, `bc_s831683_mpnn8_S15D`, `bc_s831683_mpnn9_S15D`, `bc_s831683_mpnn9_WT` (3 each); `bc_s360518_mpnn9_A22D` (1).
- artifact: The GENERATED:BASIS-TABLE 28 lines later (line 1172) gives the same design 2 in its "binder histidines" column: "| d2c_mpnn13_S88D_serasp | 4.572 | **3.526** | 2 | 0.979 |". §11.7 line 1569 agrees: "It carries two binder histidines and the full 621 aa ECD as its target." The CSV sequence contains 2 H. analysis/01-egfr/ph_sensitivity.json gives d2c n_his=19 with 17 target histidines, i.e. 2 on the 
- matters: bin/gen_methods_submission.py --check passes, so the generated-block gate cannot see this: both blocks regenerate from the same code and the code has two different definitions. The inflated 14 is produced by identifying the target construct by its histidine count — the same construct-identification bug class §7 line 689 lists ("target construct identified by chain length | silently dropped 32 desi
- gate: Make binder_his_sentence() call the same per-design `nb` that load() computes from partner=="binder", and add a selftest asserting the BINDER-HIS parentheticals equal the BASIS-TABLE "binder histidines" column design-for-design.

### [contradiction] §11.2's interpretive text places rimA02 at rank 6; the generated rank table in the same subsection places it at rank 12
- `submissions/01-egfr-METHODS.md`:1239
- verdict: CONFIRMED
- quoted: **The cost, stated:** rimA02 carries the second-highest honest pH ratio in the submission and
sits at rank 6, below a design at 1.774×. If the organisers rank strictly on the primary
objective, this ordering costs us.
- artifact: The GENERATED:RANK-TABLE at line 1224 reads "| 12 | `rimA02_d3_rimA_14_vhh` | nanobody | rimA02_d3_rimA_14 (VHH) | 129 | **4.838** | ... |". Its ratio of 4.838x is also the third-highest in the CSV, not the second: rimA01_r15_L133E reads 5.656 and c5_cf_short reads 5.546.
- matters: The sentence is the submission's stated justification for where its highest-ratio antibody-format design sits in the graded order. It names a rank six positions off the uploaded order and miscounts the design's own standing on the primary objective, in a paragraph whose whole purpose is to own the ranking cost.
- gate: The document already warns against rank references (line 1307: "Designs are referred to here by NAME rather than by rank"). Extend bin/gen_methods_submission.py --check to parse every "rank N" / "ranks N" in the prose and fail when N does not match the CSV position of the design named in the same se

### [contradiction] §4.4b gives the crop/ECD split two ways in four lines; the artifact says 17 of 18, and the wrong number propagates to limitation 15
- `submissions/01-egfr-METHODS.md`:537
- verdict: CONFIRMED
- quoted: **Why this matters to the submission: 16 of the 18 submitted designs were scored against that
crop.** The one molecule in this project with a solved complex is never docked correctly on the
crop by either predictor ... We are not able to rescore the submission on full ECD
before the deadline — one design already uses it, the other seventeen would need refolding and
re-gating
- artifact: 16 of 18 on the crop implies 2 on full ECD; "one design already uses it, the other seventeen" implies 1 and 17. analysis/01-egfr/ph_sensitivity.json resolves it: exactly one submitted design, `d2c_mpnn13_S88D_serasp`, has 17 target histidines (full 621 aa ECD); the other 17 all have 5 (the d3 crop). So 17 of 18 are on the crop. §13 limitation 15 line 1810 repeats the wrong figure: "16 of 18 submit
- matters: This is the headline exposure of the project's strongest negative result — neither ESMFold2 nor Chai-1 recovers a single crystallographic contact on the crop. The submission understates how much of itself is exposed to it, in the section that raises the concern and again in the limitation that records it.
- gate: Derive the crop/ECD split from the target histidine count or chain length in ph_sensitivity.json and emit it as a generated block; fail if any prose states a different split.

### [contradiction] §11's account of the five additions names the wrong ranks and leaves a different number of slots unused than §11.3
- `submissions/01-egfr-METHODS.md`:1132
- verdict: CONFIRMED
- quoted: Twelve were submitted on 2026-10-04; **five were added on 2026-10-05 from the reopened
exclusion pool of §10**, by a rule fixed before the result was examined ... They occupy ranks 1, 2, 3, 6 and 8. Nothing was displaced — the
submission was at 12 of 20 and the five use free slots. Three slots remain unused.
- artifact: The five CSV rows whose assessment strings say "ADDED 2026-10-05 from the reopened exclusion pool" are c5_cf_short (rank 1), c5_cr_crop_patch (rank 2), cons_gap_h370_only (rank 7), bcr_d3acid3_l60_s647537_mpnn3 (rank 8) and bcr_d3acid3_l60_s647537_mpnn11 (rank 9) per the GENERATED:RANK-TABLE — ranks 1, 2, 7, 8, 9, not 1, 2, 3, 6, 8. Ranks 3 and 6 are rimA01_r15_boltzgen_egfr_d3_rimA_20 and d2c_mpn
- matters: The paragraph is the audit trail for how the submission went from 12 to its current size without displacing anything. Its rank attribution is wrong for three of five additions and its arithmetic accounts for 17 designs, so a reviewer cannot reconstruct which rows were added when.
- gate: Generate the added-design rank list and the unused-slot count from the CSV; make gen_methods_submission.py --check fail when the prose slot count differs from 20 minus the CSV row count.

### [contradiction] §10 cites §9 for a 2% round-2 novelty baseline that §9 does not contain, and gives the same figure as 2.8% twelve lines later
- `submissions/01-egfr-METHODS.md`:936
- verdict: CONFIRMED
- quoted: Novelty is also a severe filter on this target —
§9 of this document measures Adaptyv's own round-2 set at 2% of *binders* clearing the strict
reading — so a high ratio is no guarantee any of these is submittable.
- artifact: §9 (lines 774-848) contains no measurement of Adaptyv's round-2 set and no 2% figure of any kind; its only counts are the 82/42/114/0 level table over this project's own 238 designs. The only round-2 baseline anywhere in the document is §10's own novelty table at line 948: "| Level 4 (fully de novo) | 0 — Adaptyv's round-2 baseline was 2.8% |".
- matters: The references gate confirms §9 exists but not that it supports the claim. Here the cited section has no such content, and the citing sentence and its own neighbouring table disagree on the number, in the paragraph that decides whether 25 reopened designs are submittable.
- gate: Extend the citation check so a "§N.M" pointer must find at least one distinctive numeric or phrase from the citing sentence inside the target section's text; otherwise fail.

### [contradiction] §11.7 reports the matched negative control at 5.344x on the partnered basis; the CSV and the artifact say 5.226
- `submissions/01-egfr-METHODS.md`:1442
- verdict: CONFIRMED
- quoted: The only matched negative control we have —
`bc_s831683_mpnn9_WT`, the parent of `mpnn9_S15D`, carrying no designed acid — reads
**5.344× on the partnered basis**, above two shipped designs.
- artifact: analysis/01-egfr/ph_sensitivity.json gives bc_s831683_mpnn9_WT partnered_median = 5.226. The CSV's ph_ratio_partnered_SENSITIVITY column for that row reads 5.226, and the GENERATED:RANK-TABLE at line 1228 prints 5.226.
- matters: This single number is the stated reason the partnered basis is rejected as a ranking basis and the his-only basis retained. It is quoted 0.118 above the value in the shipped CSV and in the artifact it is computed from, in the same section whose generated table prints the correct figure.
- gate: Have bin/gen_methods_submission.py --check scan §11.7 prose for the pattern `<design> ... <number>× on the partnered basis` and compare against ph_sensitivity.json.

### [contradiction] §11.6 and §11.7 give GLU133's contribution as 1.141x and 1.30x, and §10c and §11.6 give L133D as 2.116 and 2.221
- `submissions/01-egfr-METHODS.md`:1448
- verdict: CONFIRMED
- quoted: **9.84 Å** away — a desolvation shift with no electrostatic partner, the same artefact class
this project used to rule out Mechanism A — while the actual designed `GLU133` contributes
only **1.30×**.
- artifact: ph_sensitivity.json for rimA01_r15_L133E over its 20 poses gives binder:GLU133 a median ratio of 1.1409 (range 0.992-1.589). §11.6 states the pooled value twice: line 1343 ("| **Glu** (submitted) | acid, reaches | **1.141×** |") and line 1394 ("the designed GLU133 contributes 1.141×"). 1.2966 is the first single pose. Separately, §10c line 1088 gives L133D as "| B | L133D | 4.256 | 2.116 | 0.50 |"
- matters: §11.6 is built on the claim that the deeper 15-20 pose sampling superseded the five-pose numbers, and lists by name what did not survive. Two sections away the five-pose figure is still in use, and the paired L133D control — the mechanism's negative control — has two values in two sections.
- gate: Register the per-site contributions and the triad variant ratios as generated values pulled from ph_sensitivity.json, so §10c, §11.6 and §11.7 cannot drift apart.

### [contradiction] §11.5 says check_discards.py reads the candidate JSON and not the emitted CSV; §10 says it was fixed to read the CSV, and the code reads the CSV
- `submissions/01-egfr-METHODS.md`:1323
- verdict: CONFIRMED
- quoted: One consistency note recorded rather than papered over: `bin/check_discards.py` reads the
30-design candidate JSON, not the emitted CSV, and still compares on the target-only ratio. It
therefore flags a superset of what the shipped bar would flag — the safe direction — and we left
it alone rather than edit a gate at submission time.
- artifact: bin/check_discards.py line 50: `SUB_CSV = "submissions/01-egfr.csv"`; line 35 comments "Until 2026-10-05 this read analysis/01-egfr/submission_final.json"; submission_final.json is retained only as "source of per-design metadata" (line 51). METHODS §10 line 883 states the opposite: "Fixed; it now reads the graded CSV and matches by sequence, and reports 23 warns against a 2.289× bar rather than 45
- matters: The pre-submission gate is the project's keystone safeguard, and two sections of the methods document describe its input differently — one of them claiming a known defect was left unfixed at submission time when the code shows it was fixed.
- gate: Have gate_sweep.py print the file check_discards.py actually opened and have the generator inject that path into the prose, rather than describing the gate by hand.

### [contradiction] §4.4b cites §12 for a warning §12 does not contain and claims submitted designs score 0.0000 on one species; none do
- `submissions/01-egfr-METHODS.md`:496
- verdict: CONFIRMED
- quoted: 1. **A score of 0.0000 on this instrument does not mean "no interface".** It can mean "a
   correctly reproduced crystallographic interface that the predictor is not confident about".
   This matters directly: several submitted designs carry 0.0000 on one species, and §12 of
   this document already warns that an absent measurement must not read as a measured zero.
- artifact: §12 is "Declarations" (AI assistance, human review, provenance, tools, structures, funding) and contains no statement about absent measurements or zeros; the warning lives in §4.5 and in §13's limitations. And no submitted row carries 0.0000: the CSV's ipsae_min_human values run 0.1320-0.8077 and ipsae_min_mouse 0.1681-0.8037, with no zero in either column.
- matters: Both halves of the sentence fail. The pointer sends a reviewer to a section that says nothing of the kind (it is off by one — §13 is Limitations), and the factual premise that makes the §4.4b result "matter directly" to this submission is contradicted by the graded CSV.
- gate: Add two checks: a citation check that the target section contains supporting text, and an emitter assertion that any claim about zero-valued affinity cells matches the count of 0.0000 cells actually in the CSV.

### [contradiction] §4.2 says two of the three VHH candidates are novelty-ineligible; the novelty artifact clears two of three and the submission ships both
- `submissions/01-egfr-METHODS.md`:329
- verdict: CONFIRMED
- quoted: any VHH design is **inadequately assessed by this pipeline**. Our three VHH-framework
candidates are marked as such rather than rejected for low scores. **§4.5 closes this: against
framework-preserving CDR decoys the blind spot is total, not partial.** (Two are separately
ineligible on the organisers' novelty rule at 77.5% and 74.5% sequence identity to solved
structures — an upload-time hard gate
- artifact: analysis/01-egfr/novelty_ANTIBODY.tsv holds exactly three VHH records and clears two: boltzgen_egfr_h370_020 (antibody True, CDRH3 0.2727, level 4, clears_gate True) and rimA02__d3_rimA_14 (antibody True, CDRH3 0.1364, level 3, clears_gate True). Only boltzgen_egfr_2site_009 fails (antibody False, CDRH3 undelimited, global_id 0.7593). The 74.5% figure appears in outbox/PREREGISTRATION.md line 137 
- matters: If two of three VHH candidates were ineligible, only one could ship — two do. The section states a hard upload-time gate as having disqualified designs that the artifact clears, which mischaracterises the submission's actual eligibility exposure (which §9 locates in h370_020's classification, not in two identity failures).
- gate: Generate the VHH eligibility sentence from novelty_ANTIBODY.tsv (n records, n clearing, and the global_id of each failure), so the count and the identities cannot be typed by hand.

### [contradiction] §11.7 states the min-of-three-bases check over seventeen designs and reports "17 of 17"; the emitter's own data gives 18 of 18, and README says twelve
- `submissions/01-egfr-METHODS.md`:1434
- verdict: CONFIRMED
- quoted: **Why the shipped order is nevertheless retained.** The histidine-only value is the
**minimum of the three bases for all seventeen designs** — re-checked by the emitter on every
run rather than remembered, and it reports 17 of 17.
- artifact: Recomputed from analysis/01-egfr/ph_sensitivity.json over the 18 CSV rows, his-only is the minimum of {his-only, all-site, partnered} for 18 of 18. The same subsection contradicts its own denominator twice: line 1424-1425 "so all **17** are on one footing; across the 17 the analysis covers **158 human-leg poses**" against line 1470 "Across all **165 poses of the eighteen submitted designs** the co
- matters: This check is the stated justification for the shipped ranking basis — that it is the conservative envelope under one uniform rule. It is quoted with three different denominators across two documents, and the sentence claims a live emitter re-check whose reported count (17 of 17) does not match the data the emitter reads.
- gate: Have the emitter write its min-of-three pass count to a JSON artifact and have gen_methods_submission.py inject "N of N" and the pose total as generated values.

### [contradiction] README's "Still open" list says two checks are not done that its own bullets and METHODS report as completed
- `README.md`:188
- verdict: CONFIRMED
- quoted: **Still open.** The partner-deletion free leg is a fixed-conformation diagnostic, not a
measurement of the apo state. Full-ECD, glycan and receptor-state checks have not been applied
to the finalist footprints. The rAC1 comparison needs structural contact recovery against
4UIP rather than predicted confidence attached to crystallographic coordinates. These are
recorded in METHODS §13 rather than r
- artifact: METHODS §10b ("Added 2026-10-05 at the reviewer's request") applies exactly the full-ECD, glycan and receptor-state checks to the finalist footprints, with analysis/01-egfr/finalist_footprints.json as its artifact. METHODS §4.4b ("Added 2026-10-05") runs the rAC1 structural contact recovery against 4UIP via bin/rac1_contact_recovery.py, explicitly "purely geometric; no ipSAE or PAE value is attach
- matters: A method described as done in one place and as not done in another, 15 lines apart in the same file. These are three of the external reviewer's named requests, and the README tells him they are unresolved while the methods document reports results for all three.
- gate: Derive README's "Still open" list from METHODS §13 by cross-referencing which limitations carry a "not done / unresolved" marker; fail if a bullet names work for which an analysis/ artifact and a §N.M write-up both exist.

### [contradiction] README asserts binder-histidine and near-duplicate-pair counts that METHODS lists as superseded hand counts
- `README.md`:179
- verdict: CONFIRMED
- quoted: - **Two corrections to this day's own work**, both caught before release: 6ARU is the
  cetuximab-Fab complex in the *tethered* conformation ... and the count of
  designs carrying binder histidines is ten of seventeen, not the eight that was carried forward
  without recounting.
- **The five additions created two new near-duplicate pairs**, so 8 of 17 designs now sit in a
  pair differing by one 
- artifact: METHODS §11.1 line 1146-1149 lists "ten of seventeen" as one of five wrong hand counts and replaces it with a generated block reading "nine of the eighteen submitted designs carry at least one histidine of their own". METHODS §11.3 line 1253-1264 records three pairs, not four: "**Three pairs of submitted designs exceed 90% sequence identity, and all three are declared parent/mutant comparisons** .
- matters: README's correction log presents as the fixed values two counts the methods document moved past the same day, including the eligibility-relevant near-duplicate count that drives the "avoid filling slots with nearly identical variants" instruction. A reviewer auditing the near-duplicate exposure gets 8 of 17 from README and 6 of 18 from METHODS.
- gate: Emit the binder-histidine census and the >90%-identity pair census into README from the same generator functions that write the METHODS blocks, and add README to gen_methods_submission.py --check.

### [contradiction] CONTROL-TABLE cites METHODS §10 for 38 rejected molecules, a figure §10 was rebuilt to retract
- `outbox/CONTROL-TABLE.md`:220
- verdict: CONFIRMED
- quoted: **And the ordering consequence we have not honoured.** Your instruction was to assess pose and
protonation separately *before* using the pH gate to discard candidates. The gate has been used
as a discard filter throughout — §10 of the methods document lists 38 molecules rejected on it.
- artifact: METHODS §10 line 857-864: "It held two lists and never reconciled them. One was 45 run names that `bin/check_discards.py` warns on ... described here as \"38 distinct molecules\" because the gate is name-keyed on the discard side and alias pairs double-count." The rebuilt ledger (analysis/01-egfr/exclusion_ledger.json, verified: list_a_names 23, list_a_sequences 15, list_b_sequences 60, overlap 0,
- matters: CONTROL-TABLE goes to the external reviewer as the answer to his gate-discard concern, and it quantifies the concern with the alias-inflated name count the methods document rebuilt itself to eliminate — understating the gate-only exclusions by 25 molecules.
- gate: Have bin/exclusion_ledger.py emit a one-line summary ("N gate-only exclusions of an M-molecule union") and inject it as a generated block in both METHODS §10 and outbox/CONTROL-TABLE.md; forbid the literal "38" as an exclusion count.

### [contradiction] HANDOFF asserts all designs switch on H433; the generated census says one switches on H370 and it is rank 1
- `HANDOFF.md`:38
- verdict: CONFIRMED
- quoted: evidence will be weighted more heavily to catch it. **All 10 of our designs switch on the
  target's native H433** (48 of 50 pool-wide; zero tag). Say so.
- artifact: METHODS §1's GENERATED:SWITCH-SITE block at line 67: "17 of the 18 submitted designs switch on **H433**. The remaining 1: **H370** -- `c5_cf_short__boltzgen_egfr_cropfree_short_48` (rank 1)." analysis/01-egfr/master_rank.json records that sequence with "site": "H370".
- matters: HANDOFF line 38-39 ends "Say so" — it is an instruction to make this claim to the organisers, in the context of their His-tag warning. The claim is false for the top-ranked design, which is also the design the glycan check flags (limitation 14). Saying "all" where the generated census says 17 of 18 is the kind of blanket claim the reviewer asked to be qualified in his 2026-09-29 reply ("qualify th
- gate: Reuse the GENERATED:SWITCH-SITE block verbatim in HANDOFF instead of a hand-typed "All N", and add a check_claims.py rule forbidding "all N of our designs switch on" unless the switch-site census is unanimous.

### [contradiction] §10 says rimA01_r15 is shipped at rank 1; the rank table puts it at rank 3, and §9/README attribute ranks 1-4 to the wrong generator
- `submissions/01-egfr-METHODS.md`:877
- verdict: CONFIRMED
- quoted: Two of the 45 warn names resolved to `rimA01_r15`, which is **shipped at rank 1** — including
its own exact name. That is the same error §10 previously caught once by hand (citing 5.630×,
which is `bc_s360518_mpnn9_A22D` at rank 2) now found systematically.
- artifact: The GENERATED:RANK-TABLE places rimA01_r15_boltzgen_egfr_d3_rimA_20 at rank 3 (line 1215) and bc_s360518_mpnn9_A22D at rank 4 (line 1216). METHODS §1's generated block names c5_cf_short__boltzgen_egfr_cropfree_short_48 as rank 1, and §10 line 921-922 itself says c5_cf_short "would rank second of everything on this basis" though it is shipped first. §9 lines 801-802 and 811 carry the same stale att
- matters: Four separate rank attributions across §9, §10 and README describe an ordering the submission no longer has. §9's version matters most: it is the sentence that tells a grader which ranks would have been lost to the novelty-implementation bug, and it names ranks held by designs from a different generator.
- gate: Enforce the document's own stated convention (line 1307: refer to designs by name, not rank) with a check that fails on any "rank N" in prose outside the generated blocks, or that validates each against the CSV.

### [arm-accountability] The Chai-1 independent-validation arm ran 13 complexes; the document reports 1, and the omitted values are the unfavourable ones
- `submissions/01-egfr-METHODS.md`:1283
- verdict: CONFIRMED
- quoted: It also carries the strongest independent corroboration in the submission: Chai-1 ipTM **0.838** with a **39-residue** interface, the joint-largest in the validation set and above both working positives (cetuximab scFv 0.793, human EGF 0.500)
- artifact: `runs/chai1/w4_indep/` holds 13 complexes x 5 models = 65 .cif files (mtimes 2026-10-04 23:55 to 2026-10-05 00:38). Re-running `bin/chai_interface_summary.py` emits all 13 ipTM medians and interface sizes: fin_ss_bc_s831683...routeA 0.838/39res (the one quoted), fin_bc_s360518_mpnn9_A22D 0.788/39, fin_d2c_mpnn13_S88D_serasp 0.817/30, fin_rimA02_d3_rimA_14_vhh 0.440/22, fin_rimA01_r15 0.332/24, **f
- matters: `bin/chai_interface_summary.py`'s docstring states the arm's purpose in the reviewer's own words: "run a bounded independent structure-prediction check on the failed positives and a diverse finalist subset; ESMFold2-Fast alone is not an independent validation," and its title asks "Does Chai-1 form an interface where ESMFold2 did not?" The arm answers yes on both failed positives -- Chai builds 35-
- gate: Add a `GENERATED:CHAI` block to `bin/gen_methods_submission.py` that emits the full `chai_interface_summary.json` table (tag, n_models, ipTM median, interface residues, clash count) and add a `gate_sweep.py` gate asserting (a) the artifact is regenerated from the current `runs/chai1/w4_indep/*` dire

### [arm-accountability] A design that is not in the submission is named in the pKa top-set and in the human-review attestation
- `submissions/01-egfr-METHODS.md`:1499
- verdict: CONFIRMED
- quoted: five designs (`rimA01_r15_L133E`, `c5_cf_short…_48`, `sd_d2c…_T65D`, `rimA02_d3_rimA_14_vhh`, `rimA01_r15`) hold a top-three slot in 25–50% of draws and the rest hold it in 0–18%
- artifact: `sd_d2c_101_l147_s144898_m_T65D` appears 0 times in `submissions/01-egfr.csv` (grep -c = 0) and is absent from `ph_pka_perturbation.json`. The artifact's five designs with top3_fraction >= 0.25 are `rimA01_r15_L133E` 0.600, `c5_cf_short...48` 0.487, `rimA02_d3_rimA_14_vhh` 0.427, **`d2c_mpnn13_S88D_serasp` 0.275** and `rimA01_r15_boltzgen_egfr_d3_rimA_20` 0.263 -- so the actual range is 26-60%, no
- matters: Sec 12 is the human-attestation section, written because the organisers ask for it, and it explicitly warns that the count "must not be inflated by restating a count." As written it attests review of a sequence that was withdrawn and leaves two of the actually-submitted rows (`bcr_d3acid3_l60_s647537_mpnn3`, `_mpnn11`) outside both the twelve-design and five-design attestations. In Sec 11.7 the sa
- gate: Add a `check_claims.py` rule that extracts every backticked design name from METHODS and asserts each is present in `submissions/01-egfr.csv`, in `master_rank.json`, or on an explicit `HISTORICAL_NAMES` allowlist annotated with the section that legitimately discusses it -- and that no name on that a

### [arm-accountability] Expression QC covers 10 of the 18 submitted rows; Sec 12 declares the gap for five sequences, one of which is not submitted
- `submissions/01-egfr-METHODS.md`:1711
- verdict: CONFIRMED
- quoted: What has **not** been done for them: expression QC, which was only ever run on the original candidate set.
- artifact: `analysis/01-egfr/express_qc.tsv` has 20 rows. Joined against the 18 names in `submissions/01-egfr.csv`, **8 submitted designs have no express-QC row**: `c5_cf_short__boltzgen_egfr_cropfree_short_48`, `c5_cr_crop_patch__boltzgen_egfr_crop_patch_05`, `bc_s360518_mpnn9_A22D`, `ss_bc_s831683_mpnn6_S15D_S62H_routeA`, `cons_gap_h370_only__boltzgen_egfr_h370_018`, `bcr_d3acid3_l60_s647537_mpnn3`, `bcr_d
- matters: Expression is one of the three outcomes PREREGISTRATION Sec 2.1 prespecifies separately, on PK's instruction to "prespecify expression / binding / pH selectivity separately," and Sec 2 of the pipeline (line 113) lists express QC as a pipeline stage. A grader reading Sec 12 concludes the gap is confined to the five late additions; in fact it reaches four more rows, including `bc_s360518_mpnn9_A22D`
- gate: Add a `GENERATED:EXPRESS-QC-COVERAGE` block to `bin/gen_methods_submission.py` that joins `submissions/01-egfr.csv` against `analysis/01-egfr/express_qc.tsv` and emits the coverage count plus the explicit list of uncovered submitted names, so the sentence cannot drift from the join. Generalise to on

### [overclaim] "All 8 designs carrying binder histidines move" — nine carry them, and the ninth has no apo structure at all
- `submissions/01-egfr-METHODS.md`:1537
- verdict: CONFIRMED
- quoted: **Internal check first.** All **8** designs carrying binder histidines move; all **9** carrying none move by exactly nothing. The free-leg choice affects only the designs whose own histidines enter the product, which is what it should do and is evidence the comparison isolates what it claims to.
- artifact: The GENERATED:BINDER-HIS block (line 1144) says nine of the eighteen designs carry at least one binder histidine, and analysis/01-egfr/ph_sensitivity.json confirms 9 with n_his above the target's 5. analysis/01-egfr/ph_apo_freeleg.json has 18 rows: 8 move, 10 do not — and `ss_bc_s831683_mpnn6_S15D_S62H_routeA`, which carries 4 binder histidines, is not a non-mover but an error row: {"error": "no a
- matters: The sentence is offered as positive evidence that the apo comparison isolates the binder leg, and the one design that would test it hardest (4 binder histidines) is silently counted out: 8 + 9 = 17 against a shipped 18. §11.8 discloses the identical gap for the relaxed arm (line 1616, "entered the submission after the relax queue was built — a real gap, not a filtered one") but §11.7 does not disc
- gate: ph_apo_freeleg.py should emit explicit counts (n_moved / n_zero_binder_his / n_error) and gen_methods_submission.py should generate that sentence; a row with an "error" key must never fall into a "moved by exactly nothing" bucket.

### [overclaim] §11.7 understates basis instability: "up to 8 ranks" is 12, and the named example's own rank range contradicts the sentence
- `submissions/01-egfr-METHODS.md`:1429
- verdict: CONFIRMED
- quoted: is **+0.000** — the two orderings are uncorrelated. Designs move by up to **8 ranks** (`rimA02_d3_rimA_14_vhh`: 2nd on his-only, 10th on partnered).
- artifact: Ranking the 18 CSV rows on analysis/01-egfr/ph_sensitivity.json: the largest his-only→partnered shift is `bc_d3acid_l65_s831683_mpnn11`, 17th → 5th = 12 ranks. `rimA02_d3_rimA_14_vhh` is 3rd on his-only and 13th on partnered (10 ranks), which is exactly what its own CSV ph_rank_range_across_bases cell says: "3-13". Kendall tau over the 18 is +0.046, not +0.000 (the +0.000 is a hardcoded string in 
- matters: The error runs in the submission's favour on the one axis §11.7 exists to expose — the paragraph argues the ranking basis is not established, then quotes a ceiling 4 ranks too low. The cited example is self-refuting against its own CSV row, which is the first thing a reviewer checking §11.7 would compare. The same paragraph's "+0.000" is a derived number shipped inside the graded CSV.
- gate: gen_methods_submission.py: generate the max-shift value, the exemplar design, both its ranks, and the Kendall tau from ph_sensitivity.json; drop the hardcoded 0.000/8-ranks strings from emit_submission_csv.py.

### [overclaim] The pKa-perturbation table and its three span examples are computed over 17 designs and contradict Limitation 17's gate-checked 17-of-18
- `submissions/01-egfr-METHODS.md`:1486
- verdict: CONFIRMED
- quoted: | **0.8 (PROPKA's own RMSD)** | **4 of 17** | **16 of 17** |
- artifact: analysis/01-egfr/ph_pka_perturbation.json (sigma 0.8, 18 designs): 17 of 18 span ≥5 ranks between p5 and p95 — which is what Limitation 17 (line 1816) states and what check_claims' "designs span five or more ranks" rule verifies. The three examples at lines 1489-1491 are also stale: d2c_mpnn13_S88D_serasp is p5–p95 1–18, not "ranks 1–17"; bc_s360518_mpnn9_A22D is 1–16, not "2–15"; bc_s831683_mpnn6
- matters: Two numbers for the same quantity in one document, 330 lines apart, and the one in the body is the stale 17-design version while the gate-checked one sits in the limitations. "spans ranks 1–17 — the entire submission" is now literally wrong about what the entire submission is. The rule that should catch the table cell (bin/check_claims.py:265, pattern `\*\*(\d+) of 1[78]\*\*\s*\|?\s*$`) anchors on
- gate: Fix the anchor in that rule (drop `$`, or match against the unflattened text), and generate the σ-table rows and the three span examples from ph_pka_perturbation.json.

### [overclaim] README restates the binder-histidine count as "ten of seventeen" — the count that broke five times — against the generated nine of eighteen
- `README.md`:182
- verdict: CONFIRMED
- quoted: designs carrying binder histidines is ten of seventeen, not the eight that was carried forward without recounting.
- artifact: The GENERATED:BINDER-HIS block in METHODS (line 1144) reads "**nine of the eighteen submitted designs carry at least one histidine of their own**", computed from ph_sensitivity.json over the CSV; bin/check_claims.py recomputes it (9) and verifies both digit and word forms.
- matters: METHODS line 1146-1149 documents this exact number breaking five times by hand and generates it for that reason — and the README then publishes a sixth wrong value, in a bullet whose subject is "two corrections to this day's own work". The gate's word-form rule only fires on the phrase "N of the M submitted designs carry at least one histidine", so the README's paraphrase passes. README line 184 c
- gate: Broaden the binder-histidine and near-duplicate rules to any /(\w+|\d+) of (the )?(\d+|<word>)/ within a window of "binder histidine(s)" or "sit in a pair", and put the README's summary bullets inside generated markers.

### [overclaim] 17 of 18 designs are scored on the domain-III crop, not 16 — Limitation 15 and §4.4b understate the construct exposure and contradict §4.4b's own next sentence
- `submissions/01-egfr-METHODS.md`:1810
- verdict: CONFIRMED
- quoted: explanation. 16 of 18 submitted designs are scored on the crop (§4.4b).
- artifact: analysis/01-egfr/ph_sensitivity.json: exactly one submitted design has a 17-histidine target leg (`d2c_mpnn13_S88D_serasp`, n_his 19 with 2 binder histidines); the other 17 all have a 5-histidine target, i.e. the domain-III crop (targets/egfr/egfr_d3_6aru.pdb holds exactly those 5 histidines). METHODS says so twice itself: line 541 "one design already uses it, the other seventeen would need refold
- matters: §4.4b is the section showing that neither ESMFold2 nor Chai-1 docks the one molecule with a solved complex correctly on the crop, so the crop count is the size of the submission's largest structural caveat. Stating 16 rather than 17 shrinks it, and the same section states 1-plus-17 four lines below the "16 of the 18" sentence at line 537.
- gate: gen_methods_submission.py: generate the crop/ECD split from the per-design target construct, and use the same generated value in §4.4b and Limitation 15.

### [overclaim] Limitation 24 says seven designs carry no binder histidine, then says eight, where the artifact says nine
- `submissions/01-egfr-METHODS.md`:1843
- verdict: CONFIRMED
- quoted: Seven designs carry no binder histidine — their switch is target-borne — and a binder-only relaxation cannot move them, so the eight designs whose signal sits on EGFR's own H370/H433 are exactly the ones this sensitivity check cannot cover.
- artifact: The GENERATED:BINDER-HIS block (line 1144) lists nine designs with no binder histidine of their own. analysis/01-egfr/ph_relaxed_freeleg.json summary: n_compared 10, n_no_movable_site 7, n_no_relaxed_structure 1 — the 7 are the zero-binder-histidine designs MINUS the two (`rimA01_r15_L133E`, `bcr_d3acid3_l60_s647537_mpnn3`) whose target chain was also relaxed, and the 8th uncovered design (`ss_bc_
- matters: Two different numbers for the same set inside one sentence, and the stated reason is wrong for the 8th member. The quantitative claim the limitation exists to support — that the relaxed arm's 0.842× median cannot speak for the target-borne designs — is correct; only the counts and the attribution are wrong, which is the kind of thing a reviewer who checks one number stops trusting the rest for.
- gate: Generate Limitation 24's counts from ph_relaxed_freeleg.json's summary block (it already carries n_compared / n_no_movable_site / n_no_relaxed_structure) rather than restating them.

### [overclaim] CONTROL-TABLE asserts the instrument is "anti-correlated" on antibodies two paragraphs after withdrawing exactly that inference
- `outbox/CONTROL-TABLE.md`:208
- verdict: CONFIRMED
- quoted: measured example we have, it is **anti-correlated**.
- artifact: The same file, lines 169-173: "*Second, four rows cannot establish an inverse relationship.* The observation that the tightest binder scores lowest is a four-point ordering with no replication across molecules and no error model... on its own it does not show that score runs *opposite* to affinity, and it is not offered as such." The underlying data is 4 molecules from one published series (G532 0
- matters: The correction and the claim it withdraws sit 40 lines apart in the document PK named as a prerequisite for his next decision, and the stronger wording is the one in the numbered "Three consequences" list a reader will quote. METHODS §4.5 line 587 inherits it ("on the single measured example available it points the **wrong way**"), and CONTROL-TABLE line 216-218 repeats it as "not merely uninforma
- gate: No mechanical gate fits a self-contradicting inference; the checkable version is a RULE that fails when /anti-?correlated|inverted|points the wrong way/ appears in a file that also contains its own "cannot establish an inverse relationship" retraction.

### [overclaim] README publishes the superseded leave-one-family-out range, 0.778–1.000, which CONTROL-TABLE names as the figure it replaced
- `README.md`:176
- verdict: CONFIRMED
- quoted: 0.778–1.000. No interval is reported and no effective-n is substituted into Clopper–Pearson.
- artifact: bin/control_family_balance.py, run now: "SENSITIVITY: leave-one-family-out range [0.867, 1.000] around 0.889." CONTROL-TABLE line 309-310 says so explicitly: "the sensitivity is much tighter: leave-one-family-out now spans **0.867–1.000** rather than 0.778–1.000, because no single family carries nine molecules any more."
- matters: README quotes the two-family version that PK ruled out ("A shared submitting group is a clue, not a family definition"), i.e. the number produced by the construction he rejected — and it understates the panel's robustness, so it is an underclaim that also advertises the discarded method. The adjacent raw 8/10 = 0.800 and family-balanced 0.889 in the same bullet are both current, which makes the st
- gate: check_claims RULE: assert the three control-balance figures (raw, family-balanced, LOFO range) in every document against bin/control_family_balance.py's computed values, as the "control panel human-leg count" rule already does for 8/10.

### [overclaim] Stale rank citations throughout, in the one document that declares a rule against them
- `submissions/01-egfr-METHODS.md`:1240
- verdict: CONFIRMED
- quoted: sits at rank 6, below a design at 1.774×. If the organisers rank strictly on the primary
- artifact: The GENERATED:RANK-TABLE at lines 1213-1230 puts `rimA02_d3_rimA_14_vhh` at rank 12, not 6 (its "below a design at 1.774×" is consistent with 12, since 1.774× is rank 11). Same class elsewhere, all against that generated table: line 877 "`rimA01_r15`, which is **shipped at rank 1**" (it is rank 3); line 879 "`bc_s360518_mpnn9_A22D` at rank 2" (rank 4); line 954 "a variant of shipped rank 6" (rimA0
- matters: §11.4 line 1307-1309 states the rule these violate — "Designs are referred to here by NAME rather than by rank. Rank references inside a rank-ordered file drift every time the file changes, and an audit found six of the twelve assessment strings citing the wrong design by rank for exactly that reason" — and seven prose rank references then drift anyway. Two are consequential rather than cosmetic: 
- gate: check_claims RULE: flag any /rank[- ]?\d+|ranks \d+ and \d+|ranks \d+(st|nd|rd|th)/ in prose and resolve it against the CSV order where a design name sits in the same window — the citation pass already does exactly this for pH ratios and affinities.

### [overclaim] §11.7's negative-control figure for the partnered basis is 5.344× against an artifact value of 5.226×, and "above two shipped designs" is six
- `submissions/01-egfr-METHODS.md`:1443
- verdict: CONFIRMED
- quoted: **5.344× on the partnered basis**, above two shipped designs. A basis on which the
- artifact: analysis/01-egfr/ph_sensitivity.json: bc_s831683_mpnn9_WT partnered_median = 5.2255, which is what the GENERATED:RANK-TABLE prints (5.226). Ranking all 18 rows on the partnered column, six shipped designs fall below it: rimA01_r15 4.843, cons_gap 3.859, bcr_mpnn3 3.106, bcr_mpnn11 2.960, rimA02 5.183, h370_020 2.267. The figure traces to the hardcoded docstring at bin/emit_submission_csv.py:234, w
- matters: This is the single number that disqualifies the partnered basis from ranking the submission, so its provenance matters; and "above two shipped designs" understates the evidence three-fold, which is an underclaim against the submission's own argument. The citation pass does not reach it because the number is not one of the five fields in FIELD_PATTERNS.
- gate: Add partnered_median / allsite_median to check_claims' FIELD_PATTERNS citation window, and generate this sentence's value and the below-count from ph_sensitivity.json.

### [overclaim] "Nine backbone families" against the generated ten, and §10b's footprint totals are the pre-swap 17-design values
- `submissions/01-egfr-METHODS.md`:1053
- verdict: CONFIRMED
- quoted: **The finding the checks surfaced: this submission has nine backbone families and one epitope.**
- artifact: The GENERATED:FAMILY-LIST block (line 1248) and the §11.3 heading both give TEN families over 18 designs, and bin/check_claims.py verifies the heading against gen_methods_submission.FAMILY. Nine is the count before `d3acid3_l60_s647537` entered on 10-05. The footprint totals in the same section are likewise pre-swap: line 1054 "only **58 distinct positions (mature 316–474)**" is 62 positions on th
- matters: "Nine backbone families, one epitope" is the headline of Limitation 11 (line 1791) and of README line 170, and it is the submission's own statement of its correlated-failure risk — so it should not disagree with the gate-checked family count in the next section. The footprint drift is the smaller half: the conclusion (one epitope) survives and is if anything stronger at 19 residues contacted by ≥8
- gate: Generate the §10b paragraph's union size, ≥80% count and residue list from finalist_footprints.json, and extend check_claims' family-count rule beyond the §11.3 heading to any /(\w+) backbone families/ in any document.

### [adversarial-grader] Seven prose rank references in METHODS contradict the generated rank table, in the document that says it switched to names to prevent exactly that
- `submissions/01-egfr-METHODS.md`:1240
- verdict: CONFIRMED
- quoted: sits at rank 6, below a design at 1.774x. If the organisers rank strictly on the primary objective, this ordering costs us.
- artifact: The GENERATED:RANK-TABLE eleven lines above (METHODS:1222) puts `rimA02_d3_rimA_14_vhh` at rank 12, not 6; it is also the third-highest his-only ratio (5.656, 5.546, 4.838), not the 'second-highest' the same sentence claims. Other stale references against the same table: METHODS:1199 '`bc_s831683_mpnn8_S15D` ... is now rank 11' (actual 18); METHODS:164 and 1379 call `rimA01_r15_L133E` the 'rank-1 
- matters: §11.4 states the document deliberately refers to designs by name because 'an audit found six of the twelve assessment strings citing the wrong design by rank', and the prose then does the same thing seven times — twice within eleven lines of the generated table it contradicts. A grader checking the 'cost, stated' passage finds the stated cost understated by six rank positions, which makes the docu
- gate: Add a check_claims.py rule: any `rank <n>` in METHODS prose within N characters of a backticked design name must match the RANK-TABLE, and flag bare rank references that name no design.

### [adversarial-grader] Two CSV rows state a pose count in prose that contradicts their own ph_poses_n column
- `submissions/01-egfr.csv`:17
- verdict: CONFIRMED
- quoted: It reads 3.52x as the median over 11 refold poses against that design's 5.43x
- artifact: The same row's ph_poses_n column reads 26. CSV line 16 (`bc_s831683_mpnn9_S15D`) likewise says 'Median over 5 refold poses' against its own ph_poses_n of 20. Every other row agrees with its column. README.md:138-139 records the recount that created the gap ('three submitted designs had more poses on disk than the submission was counting (n=5 -> 20, 5 -> 20, 11 -> 26)') but the CSV prose was not up
- matters: PK asked twice, in writing, for pose consistency and seed counts to be reported per finalist ('Keep five seeds and their median as the provisional primary summary, with individual values and pose consistency retained'). A grader who spot-checks the n a row claims against the n the row reports finds two of eighteen disagreeing by factors of 2.4 and 4, which undercuts the pose-count audit the projec
- gate: Generate the pose-count clause of the assessment string from ph_poses_n in bin/emit_submission_csv.py, and add a check_claims.py rule that any 'over N ... poses' in an assessment must equal that row's ph_poses_n.


## LOW

### [pk-coverage] §3.7 describes d2c_mpnn13 as binding, which PK asked to be phrased as a predicted candidate
- `/Users/harish/code/adaptyv-2026/submissions/01-egfr-METHODS.md`:239
- verdict: CONFIRMED
- quoted: `d2c_mpnn13` binds and does not switch (0.70×). A single Ser→Asp at position 88 gives
**4.57×** over 5 poses.
- artifact: PK 2026-10-03: 'Describe d2d_mpnn9 and d2c_mpnn13 as predicted binding candidates with no supported pH switch, not established cross-species binders.' The project complied everywhere else: all 18 CSV assessment strings open 'COMPUTATIONAL CANDIDATE -- not shown to bind', and the d2c row reads 'parent wild-type 0.70x (no switch)'. The only measurement behind 'binds' is an ipSAE_min of 0.603 human /
- matters: It is the opening sentence of §3.7, the section the document calls 'the only causal result in the project', and it states as fact the one thing PK named this design in order to have withdrawn. Small in isolation; it is the kind of unqualified verb the rest of the submission was systematically scrubbed of, surviving in the highest-traffic result section.
- gate: Add 'binds' / 'is a binder' as banned verbs in check_claims.py when the subject is a design name, requiring 'predicted'/'computational candidate' phrasing — the rule already exists in spirit for the CSV assessment strings.

### [arithmetic] §11.7 says the emitter "reports 17 of 17" for the conservative-envelope check; it reports 18 of 18
- `submissions/01-egfr-METHODS.md`:1438
- verdict: CONFIRMED
- quoted: The histidine-only value is the
**minimum of the three bases for all seventeen designs** — re-checked by the emitter on every
run rather than remembered, and it reports 17 of 17.
- artifact: bin/emit_substitution_csv.py lines 519–541 compute n_env over every scored row and print "min(his_only, allsite, partnered) for {n_env} of {len(scored)}"; len(scored) is 18. Recomputed from ph_sensitivity.json, the his-only median is the minimum of the three for all 18 shipped designs, so the emitter now prints 18 of 18.
- matters: The sentence's whole point is that the figure is machine-re-checked rather than remembered, and the quoted output is a remembered one from before the 18th design was restored. The underlying property still holds, which is why no gate fires — but a reviewer who runs the emitter to verify the one claim the document says is auto-verified gets a different number than the document reports.
- gate: Capture the emitter's envelope line into a GENERATED block in §11.7 instead of transcribing it; the emitter already prints exactly the sentence needed.

### [contradiction] §4.3 blames §3.1 for the 609-residue ECD figure; §3.1 contains no residue count
- `submissions/01-egfr-METHODS.md`:342
- verdict: CONFIRMED
- quoted: *(This project writes the full ectodomain as 609 residues in some places and 621 in others —
including §3.1 above, against this section's own argument that a bar must come from the matching
construct. **621 is correct** for what was actually folded
- artifact: §3.1 (lines 134-139) reads in full: "Against the full ectodomain, 4/30 designs touched the declared patch. Against a 170-residue domain III crop, **27/30 and 29/30** did." It gives no residue count for the ectodomain. A grep for "609" across METHODS, README, HANDOFF and CONTROL-TABLE returns only lines 342 and 346 — this parenthetical itself — so the "some places" the note describes do not exist i
- matters: The note is a self-correction that names a specific location for an inconsistency that is not there, which is how a reader loses confidence in the correction history. It also asserts a pattern ("in some places") that no other line in the deliverable set instantiates.
- gate: Make the references gate assert that when a pointer accuses a target section of containing a specific value, that value appears in the target section; otherwise fail.

### [contradiction] §3.5 cites §3.7 for the S60D result; §3.7 does not mention S60D
- `submissions/01-egfr-METHODS.md`:204
- verdict: CONFIRMED
- quoted: the one design built deliberately to put a carboxylate on H370, `S60D`, returned
1.11× — inside PROPKA's noise — while dropping human ipSAE from 0.654 to **0.012** (§3.7). H370
- artifact: S60D appears in METHODS at exactly three lines: 205 (this one), 749 (§8.1) and 1409 (§11.7). §3.7 (lines 234-287) is "One mutation converts a binder into a switch, and it replicates four times" and covers only d2c_mpnn13 S88D and the four d3acid_l65_s831683 S15D sequences. It contains neither S60D nor the 1.11x / 0.654 / 0.012 figures.
- matters: §3.5 is the retraction of the project's former central result, and the narrower claim it says survives ("The trade-off is real at **H370** and only there") rests on this S60D datum. The pointer offered as its support leads to a section with no S60D content, so the one surviving claim from a retracted result is uncited.
- gate: Same citation-content check: require the target section to contain the design name or a figure from the citing clause.

### [contradiction] Limitation 24 gives the no-binder-histidine count as seven and then eight in one sentence; the generated census says nine
- `submissions/01-egfr-METHODS.md`:1843
- verdict: CONFIRMED
- quoted: 24. **The relaxed arm reaches 10 of 18 designs, and the gap is not random.** Seven designs carry
    no binder histidine — their switch is target-borne — and a binder-only relaxation cannot
    move them, so the eight designs whose signal sits on EGFR's own H370/H433 are exactly the
    ones this sensitivity check cannot cover.
- artifact: The GENERATED:BINDER-HIS block at line 1144 ends "The other nine carry none." Recomputed from ph_sensitivity.json, 9 of the 18 submitted designs have zero binder histidines. §11.8's own coverage table (line 1612-1616) puts 7 in the "relaxed, no movable site" bucket — correct, because 2 of the 9 received a relaxed target instead — and 1 in "no relaxed structure". §11.7 line 1537 gives yet another p
- matters: One sentence states the same quantity as seven and as eight, and both differ from the generated census in the same document. The limitation's point — that the designs the sensitivity check cannot reach are exactly the target-borne ones — depends on getting this partition right.
- gate: Make the §11.8 coverage buckets and limitation 24's counts generated from the same per-design binder-histidine census that feeds the BASIS-TABLE.

### [contradiction] Limitation 17 cites §11.7 for "17 of 18 designs span five or more ranks"; §11.7's table says 16 of 17
- `submissions/01-egfr-METHODS.md`:1815
- verdict: CONFIRMED
- quoted: 17. **The pH ratio cannot order this submission.** Under PROPKA's own reported accuracy
    (±0.8 pKa units) 17 of 18 designs span five or more ranks and no design holds a
    top-three slot in more than half of draws (§11.7).
- artifact: §11.7's perturbation table at line 1486 reads "| **0.8 (PROPKA's own RMSD)** | **4 of 17** | **16 of 17** |", and line 1480 says the analysis "re-ranks all 17". 17 of 17 is the σ=1.2 row, not the σ=0.8 row. §13 limitation 9 (line 1773-1774) has a matching denominator failure in the other direction: "two of ten rows depend on Adaptyv's ANARCI calling them antibodies (§9)" against an 18-row CSV — an
- matters: The limitations list is where a grader looks for the project's own honest bounds, and these two entries restate their cited sections with changed numerators and denominators. Limitation 17 is the submission's central caveat about its own ranking.
- gate: Make §13's numeric limitations generated from the artifacts they cite (ph_pka_perturbation output, the CSV row count), the way limitations 19 and 22 already are.

