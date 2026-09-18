# Patches for `biomodals`

`biomodals/` is a git clone and is gitignored, so edits there are **not version-controlled**. These diffs are the record. A `git pull` in that directory silently reverts everything and the failures come back.

## Re-apply after any clone or pull

```bash
cd ~/code/adaptyv-2026/biomodals
for p in ../patches/modal_*.patch; do patch -p0 < "$p"; done
../bin/check-pins.sh
```

## What each one fixes

**`modal_boltzgen.patch`** — pins Python 3.12. Modal's default is now 3.14, and BoltzGen's `numba` dependency supports `>=3.10,<3.14`, so the image build fails outright.

**`modal_chai1.patch`** — same Python pin, applied pre-emptively. Chai sits in the documented BoltzGen → Chai → QC pipeline and would have failed identically at the next step.

**`modal_esmfold2.patch`** — three separate fixes, all required to run at all:
1. **Dependency ref bumped** from `c94ed8d` to `43b4548b…` (esm 3.4.1.post1). The old commit depends on `github.com/Biohub/transformers`, which returns **404** — private or deleted. Upstream replaced it with PyPI `transformers>=4.57.6`.
2. **Model import** changed from `transformers.models.esmfold2.modeling_esmfold2.ESMFold2Model` to `esm.models.esmfold2.EsmFold2Model`. The transformers path exists only in that dead fork and in upstream 5.16.0.dev0+. Note the capitalisation differs.
3. **ipSAE sidecar emitted.** The wrapper writes only a scalar mean pLDDT, pTM and ipTM, and discards the PAE matrix that ipSAE requires. Adds `<name>_sample_<i>_ipsae.json` with the full PAE and per-token pLDDT. Falls back to printing the sample object's attributes if no PAE is found.
4. Also makes `ESMFOLD2_HF_REPO` / `ESMFOLD2_HF_REVISION` env-selectable so the protocol's **Fast** arm can run, not just Full.

## Upstreaming

All four are genuine bugs in `biomodals` or stale pins, not local preferences. Worth filing upstream — the Python 3.14 default will break these for everyone.
