# biomodals wrapper gotchas

Verified by reading the scripts on 2026-09-18. These are the things that cost money or silently produce the wrong result.

## The env-var names are inconsistent

| Script | GPU override | Timeout override | Default GPU | Default timeout |
|---|---|---|---|---|
| `modal_boltzgen.py` | `GPU` | `TIMEOUT` | L40S | 120 min |
| `modal_bindcraft.py` | `GPU` | `TIMEOUT` | L40S | 300 min |
| `modal_protenix.py` | `GPU` | `TIMEOUT` | L40S | 60 min |
| **`modal_esmfold2.py`** | **`MODAL_GPU`** | **`MODAL_TIMEOUT`** | **A100-40GB** | 30 min |
| `modal_esmfold2_binder_design.py` | `GPU` | `TIMEOUT` | **H100** | 60 min |

**`GPU=L40S modal run modal_esmfold2.py ...` is silently ignored.** It reads `MODAL_GPU`, so the override does nothing and the job runs on an A100-40GB at $2.10/hr instead of an L40S at $1.95. Not ruinous per run, but it is the kind of thing that quietly compounds across a campaign, and the same mistake in the other direction would be worse.

`modal_esmfold2_binder_design.py` defaults to an **H100 at $3.95/hr with a 60-minute ceiling**, so one hung job is about $4.

## ESMFold2 defaults do not match Anthropic's protocol

| Parameter | Wrapper default | Anthropic's protocol |
|---|---|---|
| `num_loops` | 3 | **10** |
| `num_sampling_steps` | 50 | **68** |
| `num_diffusion_samples` | 1 | 1 ✓ |

Running the wrapper as shipped gives a **weaker instrument than the published one**, and the ranking instrument is where the competition is decided. Override both explicitly on every scoring run.

## Protenix

- `seeds` takes a comma-separated list, which is how to satisfy the five-seed requirement in one call.
- `use_msa` defaults to **true**. Partially resolved 2026-09-18: **the wrapper itself contains no hardcoded MSA endpoint** and simply passes `--use_msa` through to the Protenix CLI, so any server call happens inside the upstream package, not here. Still worth confirming against Protenix's own docs before a scaled run, since a public MMseqs2 backend would risk an IP ban under heavy use. Safest default either way: stage target MSAs once per target and reuse them. The binder chain is single-sequence regardless.
- First call JIT-compiles a CUDA kernel, 4–6 minutes on an H100. Point `TORCH_EXTENSIONS_DIR` at a persistent Modal Volume or every container pays it again.

## Logging

**Do not pipe `modal run` through `tail`.** It buffers until the process exits, so the log stays empty for the whole run and any monitor watching it is blind. Stream `modal app logs <app-id>` instead, and note that stream can drop on its own without the job failing.

## Modal's default Python is now 3.14, and it breaks builds

**This killed the first 1ALU run.** `Image.debian_slim()` with no `python_version` gets whatever Modal currently defaults to, which is **3.14.2**. BoltzGen depends on `numba`, which supports `>=3.10,<3.14`, so the image build failed with:

```
RuntimeError: Cannot install on Python version 3.14.2; only versions >=3.10,<3.14 are supported.
ERROR: Failed to build 'numba' when getting requirements to build wheel
```

Audit of all 24 wrappers as of 2026-09-18:

| Status | Scripts |
|---|---|
| **Was unpinned, in our critical path** | `modal_boltzgen.py`, `modal_chai1.py` — **both patched locally to 3.12** |
| Unpinned, not needed | `faspr`, `minimap2`, `nextflow_example`, `rso`, `sasa` |
| Already pinned | alphafold 3.11, bindcraft 3.11, boltz 3.11, diffdock 3.10, esm2 3.10, ligandmpnn 3.11, germinal 3.10, mber 3.11, pdb2png 3.11, iggm 3.10, tmol 3.12, usalign 3.11, esmfold2_binder_design 3.12 |
| Immune, fixed registry base | `modal_esmfold2.py`, `modal_protenix.py`, `modal_esmfold2_binder_design.py` |

Originals kept as `*.py.orig`. These patches are local to our clone, so **a `git pull` in `biomodals/` will revert them.** Re-apply after any update, and check for newly unpinned scripts.

## `modal run` exit codes lie when piped

`modal run ... | tail` reported **exit code 0 on a failed build**, because the shell reports the pipe's status. Never pipe it, or set `pipefail`. This is how a hard failure gets mistaken for success.
