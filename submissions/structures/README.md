# Predicted complexes — 9 of the 16 submitted designs, plus one withdrawn

**Chain A is the binder, chain B is the target.** These ten files were cut when the
submission stood at ten designs and have not been regenerated since, so coverage is an
artifact of that history rather than a judgement about any design.

⚠️ **`rimA01_r15_boltzgen_egfr_d3_rimA_20.cif` is no longer a submitted design.** It and
its mutant `rimA01_r15_L133E` were removed on 2026-10-06 when the organisers’ novelty
check rejected them — having accepted the same two sequences the day before (METHODS
limitation 37). The file is kept because the pose and its scores are real and were cited
in the methods document; it is simply not part of the submission. So nine of these ten
cover a submitted design, and seven submitted designs have no structure here.

| file | target construct | ipSAE_min of this pose | poses available |
|---|---|---|---|
| `bc_s831683_mpnn8_S15D.cif` | d3 crop, 170 aa | 0.7760 | 5 |
| `bc_s831683_mpnn19_S15D.cif` | d3 crop | 0.8077 | 5 |
| `bc_s831683_mpnn9_S15D.cif` | d3 crop | 0.8025 | 5 |
| `bc_s831683_mpnn6_S15D.cif` | d3 crop | 0.7803 | 5 |
| `rimA02_d3_rimA_14_vhh.cif` | d3 crop | 0.2275 — **not interpretable, VHH format** | 6 |
| `rimA01_r15_boltzgen_egfr_d3_rimA_20.cif` | d3 crop | 0.6110 | 6 |
| `d2c_mpnn13_S88D_serasp.cif` | **full ECD, 621 aa** | 0.6031 | 5 |
| `bc_d3acid_l65_s831683_mpnn11.cif` | d3 crop | 0.7978 | 6 |
| `bc_s831683_mpnn9_WT.cif` | d3 crop | 0.7833 | 11 |
| `h370_020_vhh.cif` | d3 crop | 0.4171 — **not interpretable, VHH format** | 11 |

## Three things to know before using these

**Each file is the MEDIAN pose by ipSAE_min, not the best.** Shipping the best pose of five would
be selection on the outcome, which is the error this project spent a week correcting. Where the
count is even, the upper median is taken. The per-design values in `../01-egfr.csv` are medians
over all poses and will therefore differ slightly from the single pose here.

**These are ESMFold2 predictions without an MSA.** They are a hypothesis about geometry, not a
measurement. The methods document is explicit that a generator's or predictor's structure is not
evidence about where a binder sits (§2), and §4.4 shows this instrument ranks a measured
non-binder above our measured positive on the human leg.

**The two VHH-format rows should not be read from these files at all.** On the one antibody in
this project with a measured KD — G532, 294 nM, 13.26× pH switch — this pipeline scores the real
binder at 0.0135 and its non-switching comparator 18× higher, while folding the Fv itself at 0.85
(METHODS §4.5). ESMFold2 builds antibodies and fails to dock them. The 0.2275 and 0.4171 above
are reported for completeness and mean nothing.

**Numbering.** On the d3 crop, target H433 is residue **99** and H370 is residue **36**
(crop + 334 = canonical). On the full ECD, H433 is residue **409** (ECD + 24 = canonical). The
switching residue for all ten designs here is the target's native H433; across the
16-design submission 15 switch on H433 and one on H370.

The pose cache these were drawn from (6,011 scored outputs) is not published; see the
reproducibility note in the methods document. Ask if you want it.
