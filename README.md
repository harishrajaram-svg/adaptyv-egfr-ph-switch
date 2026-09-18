# Adaptyv 2026 — working directory

Code and compute for the Anthropic × Adaptyv protein design competition, Track 3.

**Notes, strategy, and the per-problem workflow live in the vault**, not here:
`~/code/context-directory/projects/anthropic-adaptyv-2026/`

This directory is for things that are large, generated, or not worth committing to a personal journal repo.

## Layout

```
biomodals/     Modal GPU wrappers for every design tool (cloned, MIT)
reference/     Anthropic's 16k-word binder-design protocol prompt
targets/       Target structures (1ALU = human IL-6, the test target)
runs/          Output from design runs
```

## Status as of 2026-09-18

| Piece | State |
|---|---|
| `adaptyv@protein-design-skills` plugin | ✅ installed, 24 skills, user scope |
| `biomodals` GPU wrappers | ✅ cloned |
| Modal CLI | ✅ v1.5.5 via `uv tool install modal` |
| Anthropic protocol prompt | ✅ `reference/`, 111 KB |
| Test target 1ALU | ✅ `targets/` |
| Modal authentication | ✅ workspace `harishrajaram-svg`, verified, $0.00 spent |
| **Proteinbase account** | 🚫 **BLOCKED ON ADAPTYV, 2026-09-18.** Both Google and GitHub return "Signups not allowed for this instance" — signups are disabled at their project level. Support email drafted to proteinbase@adaptyvbio.com. Not urgent: the account is only needed to submit, and problem 1 closes Oct 4. |
| **Compute cap** | ❌ **dashboard only — modal.com → Settings → Usage & Billing** |

## Still needed from you

**1. Set the spend cap.** There is no CLI for this — modal.com, Settings, Usage and Billing. `modal billing summary` reads current spend, `modal billing rates` reads pricing, but neither sets a limit.

**2. Create a Proteinbase account** at https://proteinbase.com/login. Google or GitHub sign-in. This is the only thing standing between you and being able to submit.

## Set the spend cap first

Track 3 is self-funded. Set a limit in the Modal dashboard under workspace settings before running anything real.

Rates confirmed live on 2026-09-18: L40S $1.95/hr, A10G $1.10, A100-80 $2.50, H100 $3.95, L4 $0.80, T4 $0.59, CPU $0.047/core/hr.

Defaults worth knowing: `modal_boltzgen.py` runs on an L40S at $1.95/hr with a 120-minute timeout, so one unattended run that hangs costs about $3.90. Both are overridable:

```
GPU=A10 TIMEOUT=30 modal run modal_boltzgen.py ...
```

Real budget is roughly $500–1,500 per target. Trimming the target before sampling is the difference between about $100 and about $30,000 on the same problem.

## First run, once Modal is authenticated

```
cd biomodals
modal run modal_boltzgen.py --help
curl -o ../targets/1ALU.pdb https://files.rcsb.org/download/1ALU.pdb   # already done
```

Then design against 1ALU end to end. It is a known target with published results, so it tells you whether the toolchain works before a real problem exists.

## Available wrappers

`boltzgen` `bindcraft` `chai1` `boltz` `esmfold2` `esmfold2_binder_design` `protenix` `alphafold` `ligandmpnn` `germinal` `afdesign` `af2rank` `rso` `anarci` `sasa` `usalign` `pdb2png` `esm2_predict_masked` `diffdock` `faspr` `iggm` `mber` `tmol` `minimap2`

**No RFdiffusion or plain ProteinMPNN wrapper exists here.** `modal_ligandmpnn.py` runs ProteinMPNN mode via `--model_type protein_mpnn`. An RFdiffusion wrapper is the one thing you would write yourself, and it is optional — PXDesign, BoltzGen and BindCraft cover generation.

## Ranking instrument

The thing that decides which designs get ordered. Score `ipSAE_min`, max over ≥5 seeds, across three arms, z-scored within target:

- `modal_esmfold2.py` — ESMFold2 Full and Fast
- `modal_protenix.py` — Protenix v2

Set `TORCH_EXTENSIONS_DIR` to a persistent Modal Volume before the first Protenix call. It JIT-compiles a CUDA kernel for 4–6 minutes and will redo it in every container otherwise.

Full rationale and thresholds: `methods-stack.md` in the vault project folder.
