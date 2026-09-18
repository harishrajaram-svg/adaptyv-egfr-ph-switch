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

## ESMFold2: two defects, both patched

### 1. The pinned `esm` commit cannot build (fixed)

`ESMFOLD2_GIT_REF = "c94ed8d"` fails at image build:

```
failed to fetch commit 3a8956fb4d4ea16b0ec8e71deef2c2909b6a5cbf
Terminating task due to error: failed to run builder command
  uv pip install 'esm @ git+https://github.com/Biohub/esm.git@c94ed8d' ...
```

Traced it: `esm` at that commit declares
`transformers @ git+https://github.com/Biohub/transformers.git@3a8956fb...`,
and **`github.com/Biohub/transformers` returns 404** — private or deleted. No retry fixes this.

Upstream already solved it. The current `esm` main branch uses plain `transformers>=4.57.6,<5.0.0`. **Patched the ref to `43b4548b86762edfa747b07d5f440aad3c33acee`** (esm 3.4.1.post1, 2026-09-16).

The HF weights repo `biohub/ESMFold2` is fine and always was — only the code package was broken.

### 2. It emits no PAE, so ipSAE cannot score it (patched)

The wrapper hand-builds its scores JSON with only `plddt` (a scalar mean), `ptm`, `iptm`, `chain_pair_iptm`. **ipSAE needs the full PAE matrix**, and it reads `plddt` as a per-token array.

So out of the box, the wrapper cannot feed the best-known ranking method — and it fails *silently*, leaving you to fall back on ipTM, which the post-mortems put barely above chance.

Patched to also emit `<name>_sample_<i>_ipsae.json` containing `pae`, per-token `plddt`, `ptm`, `iptm`. It probes four attribute names for the matrix and, if none exist, prints the sample object's full attribute list so the gap is obvious rather than silent.

**If ESMFold2 genuinely exposes no PAE, swap that arm** for Chai-1 or Boltz-2, both of which emit one.

## Only the Full model is reachable

`ESMFOLD2_HF_REPO` is hardcoded to `biohub/ESMFold2`. The protocol's three-arm ensemble also wants **ESMFold2-Fast** (`biohub/ESMFold2-Fast`), which needs either a second patched copy of the script or an env var. Not done yet.

## Scoring is free

ipSAE is pure NumPy and runs locally on CPU. Cloned to `ipsae/`, driven by `bin/ipsae_min.py`.

**`bin/ipsae_min.py` takes the MINIMUM over both alignment directions**, then the **max over seeds**. ipsae.py also prints a row labelled `max` — that is *not* the metric. On the bundled example the directions are 0.449 and 0.866; ipSAE_min is **0.449**. Using the `max` row would systematically overrate every design.

## Fix #3: the model import also pointed at the dead fork

Bumping the `esm` ref was **not sufficient**. The wrapper also did:

```python
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
```

That module exists only in the deleted Biohub transformers fork and in upstream transformers 5.16.0.dev0+. With mainline 4.x it raises `ModuleNotFoundError`. The `esm` package says so in its own `hf_adapter.py`.

Correct path, per the package's own docstring — note the capitalisation:

```python
from esm.models.esmfold2 import EsmFold2Model   # NOT ESMFold2Model
model = EsmFold2Model.from_pretrained("biohub/ESMFold2")
```

**ESMFold2 needed three separate patches to run at all**: the dependency ref, the model import, and the PAE sidecar.

## Verified working, 2026-09-18

ESMFold2 → PAE → ipSAE → ipSAE_min runs end to end.

- PAE matrix emitted: **247 × 247**, per-token pLDDT 247 values. The patch works.
- Model load **116s** per cold container, fold **15.9s** for a 247-residue complex. **Batch many designs per run**; the load dominates.
- ipSAE consumed the sidecar directly in AF3 mode (`.cif` + `.json`). The missing AF3 summary file is only a warning, not a failure.

**Gotcha in `bin/ipsae_min.py`, fixed:** ipsae.py is invoked with `cwd=struct.parent`, so paths must be resolved to absolute first or it exits silently with no output and no error.

## Protenix: the wrapper defaults to v1, not the v2 the protocol specifies

`DEFAULT_MODEL = "protenix_base_20250630_v1.0.0"`.

Anthropic's three-arm ensemble — the one that reached **0.66 macro-AP** against AlphaFold3's 0.55 — used ESMFold2-Full, ESMFold2-Fast and **Protenix v2**. Running the shipped default silently substitutes a different, older model for the third arm, and nothing warns you.

`model_name` is exposed as a CLI flag, so override it: `--model-name <v2 name>`. Confirm the exact identifier against Protenix's model registry before trusting it — a wrong name may fall back rather than error.

### Protenix cannot feed ipSAE without v2 — third arm switched to Boltz-2

The wrapper's Protenix run succeeds but writes **only** `*_summary_confidence_*.json` (pTM, ipTM, per-chain scalars). **No PAE matrix**, so ipSAE cannot score it.

The protocol's fix is `--need_atom_confidence true`. **That option does not exist in v1**:

```
Error: No such option '--need_atom_confidence'
```

The config key `configs.need_atom_confidence` is real (`runner/dumper.py`), but v1's `pred` subcommand does not expose it. It is a v2 flag, consistent with the protocol having been written for v2. Patch reverted; the wrapper is back to stock.

**Decision 2026-09-18: use Boltz-2 as the third arm instead.**

- `modal_boltz.py` passes arbitrary params straight through, so `--write_full_pae` needs **no patch**.
- **ipSAE supports Boltz natively** — `ipsae.py pae_*.npz model_0.cif 10 10` is one of its three documented input formats.
- MIT, open weights, already pinned to Python 3.11, so none of today's other defects apply.

The protocol's requirement is a multi-arm ensemble of *independent co-folders*, not those three specific models. Boltz preserves that. Protenix stays available as a fourth arm if the v2 identifier is ever resolved — its image builds and runs fine.

## Boltz-2 does not build — root cause found, fix NOT applied (timeboxed)

`modal_boltz.py` fails during image build, inside the **ColabFold** install:

```
ModuleNotFoundError: No module named 'pkg_resources'
ERROR: Failed to build 'pandas' when getting requirements to build wheel
```

**The real cause is Modal's Python 3.14 default again, one layer deeper than it looks.** The traceback path is:

```
/tmp/pip-build-env-ekk729pc/overlay/lib/python3.14/site-packages/setuptools/build_meta.py
```

Note **python3.14**, even though the image pins `python_version="3.11"`. pip's **build-isolation overlay** is built independently of the image environment, and its setuptools is new enough to have dropped `pkg_resources`, which pandas' legacy `setup.py` imports.

So adding `setuptools<81` to the image does **not** fix it — the overlay never sees the image's packages. That patch is kept (harmless, likely still necessary) but is **not sufficient on its own**.

**Candidate fix, untested:** pass `--no-build-isolation` to that specific `pip_install` so it uses the image's pinned setuptools, or pin a pandas version with a PEP 517 build backend. Either needs verification.

**Decision:** stopped here on a stated timebox. Two arms (ESMFold2 Full + Fast) are verified and sufficient to run the validation gate, which is the higher-value work. Boltz stays a known-diagnosed, unfixed third arm.

Note this is the **third distinct failure** traceable to Modal's 3.14 default (BoltzGen, Chai pre-emptively, Boltz). Any new wrapper should be assumed guilty until checked.
