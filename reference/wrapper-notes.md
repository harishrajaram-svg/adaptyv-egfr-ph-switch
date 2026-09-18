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

## Boltz-2 — FIXED 2026-09-18. Root cause was one missing argument.

**`.micromamba()` called with no args builds its conda env with Modal's CURRENT DEFAULT Python — now 3.14 — silently overriding the `python_version="3.11"` on the line directly above it.**

```python
Image.debian_slim(python_version="3.11")   # decorative
.micromamba()                              # <-- actually installs Python 3.14
```

Fix: `.micromamba(python_version="3.11")`.

Every symptom traced to this: pip ran from `/opt/conda/lib/python3.14/`, whose setuptools >=81 has dropped `pkg_resources`, which pandas' legacy `setup.py` imports. Once the interpreter was actually 3.11, the conda env shipped setuptools **80.9.0** and the original build path worked unmodified.

**Two of my earlier patches were treating symptoms.** Adding `setuptools<81` to the image did nothing (the overlay never saw it) and `--no-build-isolation` did not fix it either — it only produced a traceback that exposed the real interpreter path. An error message disappearing is not the same as a cause being found.

**Runtime config, separate issue:** Boltz refuses to run without alignments. The wrapper's FASTA→YAML converter sets no MSA field, so pass a YAML with `msa: empty` per chain for single-sequence mode. That is the correct setting anyway — binders are scored single-sequence, and it avoids `--use_msa_server` and its IP-ban risk entirely.

**Verified mechanically:** Boltz-2 emits `pae_*.npz`, which **ipSAE reads natively** (no patch needed).

### ⛔ But it FAILS the validation gate — do not ensemble it

| Complex | ESM2-Full | ESM2-Fast | **Boltz-2** |
|---|---|---|---|
| barnase + barstar | 0.8887 | 0.8888 | **0.0000** |
| four shuffles | 0.0000 | 0.0000 | 0.0000 |

It scores one of the tightest complexes known **identically to random shuffles**.

**Cause: `msa: empty`.** Boltz-2 is a co-folding model that depends on alignments. ESMFold2 is a language-model folder built for single sequences (the Fast checkpoint has no MSA encoder at all). The single-sequence setting that is correct for ESMFold2 is **crippling** for Boltz.

**Cost of including it:** the 3-arm mean drags the true binder to **0.5925 — below the 0.61 "worth ordering" threshold**. Worse than either ESMFold2 arm alone.

**Every mechanical signal was green** — build, run, PAE file, ipSAE parse, number returned. Only a known answer exposed it.

**Rule: validate every arm independently before ensembling. An arm that produces numbers is not an arm that produces signal.**

To use Boltz at all, it needs real MSAs (staged per target, not `--use_msa_server` at scale).

⚠️ **Latent conflict:** a later layer upgrades numpy to 2.4.6, breaking ColabFold's `numpy<2.0` pin. Harmless for us — we never invoke ColabFold — but it would matter if target MSAs are ever staged through it.

### (superseded) original diagnosis

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


## Multi-seed: patched, and measured

The wrapper took **one seed per call**, so the protocol's "max over >=5 distinct seeds" meant five containers and five model loads -- a real 5x.

**Patched to accept `--seed "1,2,3,4,5"` and loop inside one container**, reusing the loaded model. Outputs are tagged `_seed<N>_` so they don't collide; `bin/ipsae_min.py` groups them and reports max, min and spread.

Measured on barnase/barstar, 2026-09-18:

| | |
|---|---|
| Model load | 137s (paid once) |
| Fold, seed 1 | 5.5s |
| Folds, seeds 2-5 | 4.0s each |

**5 seeds costs ~1.16x one seed**, not 5x. There is little reason to screen at one seed except on very large pools.

### Seed noise, measured

| Seed | ipSAE_min |
|---|---|
| 1 | 0.8860 |
| 2 | 0.8848 |
| 3 | 0.8859 |
| 4 | 0.8800 |
| 5 | 0.8899 |

**Max 0.8899, min 0.8800, spread 0.0099.**

Tighter than the ~0.07 discrepancy Adaptyv reported between entrant and organiser runs -- so most of that gap was likely hardware/software, not seed choice.

**Interpretation.** 0.01 is negligible against the binder/non-binder gap (0.89 vs 0.00), so one-seed screening is defensible for a bulk pool. It is **not** negligible between candidates near a threshold: two designs at 0.62 and 0.61 are indistinguishable at one seed. Five seeds is what makes a shortlist ranking mean anything.

Measured on one strong complex. **Re-measure on a marginal design in week 1** -- spread may be wider exactly where it matters most.


## Batching: one complex per container was the real bottleneck

The wrapper folded **one complex per container**. Scoring 100 designs meant 100 model loads at ~137s each -- about **4 GPU-hours of pure loading** before any useful work. That dwarfs the seed question.

**Patched: `--input-faa` now accepts a DIRECTORY** and folds every `.faa` in it inside one container, reusing the loaded model. Outputs go to per-complex subdirectories.

Measured, 5 complexes, 1 seed, one container: **load 115.8s once**, then folds at 5.5 / 4.0 / 4.0 / 4.0 / 4.0s.

| 100 designs, 1 seed | GPU time | Cost |
|---|---|---|
| One container per design | ~3.9 h | ~$7.60 |
| Batched | ~9 min | ~$0.29 |

**~20x.** This also reframes the seed question: once batched, seeds are the marginal cost, so 5 seeds really is ~4-5x the *folding* -- but folding is now cheap in absolute terms.

## Two-tier scoring — `bin/two_tier.py`

Screen the pool at 1 seed, confirm the shortlist at 5. Verified end to end 2026-09-18:

| | Screen (1 seed) | Confirm (5 seeds) | delta | moved |
|---|---|---|---|---|
| barnase/barstar | 0.8883 | 0.8903 | +0.0020 | 0 |
| shuffles x2 | 0.0000 | 0.0000 | 0.0000 | 0 |

**The delta is expected POSITIVE.** Max-over-seeds is a *biased* estimator -- more draws means a higher max -- so a 5-seed score is systematically above a 1-seed score for the same design. Consequences:

- **Never compare raw scores across seed counts.** Compare ranks. The script says so in its output.
- **Freeze the seed count** alongside the rest of the instrument, for the same reason the instrument itself is frozen.

The `moved` column is the diagnostic: large rank movement between tiers means the 1-seed screen was noisy *at that score range*, which tells you empirically whether screening is safe on a given problem instead of assuming.

**Standing rule, in the script's docstring:** scoring never eats the optimization rounds. If a run is behind, cut seeds or sampling breadth before predict-then-redesign cycles -- those are the cheapest hit-rate gain in the literature.
