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
- `use_msa` defaults to **true**. **Verify what MSA backend it uses before any scaled run.** If it reaches a public MMseqs2 server, heavy use risks an IP ban. Stage target MSAs once per target; the binder is single-sequence regardless.
- First call JIT-compiles a CUDA kernel, 4–6 minutes on an H100. Point `TORCH_EXTENSIONS_DIR` at a persistent Modal Volume or every container pays it again.

## Logging

**Do not pipe `modal run` through `tail`.** It buffers until the process exits, so the log stays empty for the whole run and any monitor watching it is blind. Stream `modal app logs <app-id>` instead, and note that stream can drop on its own without the job failing.
