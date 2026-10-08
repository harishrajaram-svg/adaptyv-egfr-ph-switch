# Item 2d — is the 2026-10-07 ESMFold2 run comparable to the calibration it was chosen from?

**Verdict: NO, not directly. The numbers are valid; the AUC that justified the ranker is not
transferable to them without one more step.** Checked 2026-10-08 00:20 ET.

## What we ran

| | |
|---|---|
| wrapper | `biomodals/modal_esmfold2.py`, via `bin/score-esmfold2.sh` |
| checkpoint loaded | `biohub/ESMFold2` — the **Full** repo (confirmed from the run log's HF snapshot path) |
| `ESMFOLD2_HF_REPO` override | none; `score-esmfold2.sh` never sets it, so the default Full repo loaded |
| MSAs | **none, for any chain.** `grep -ci msa` is **0** in both `modal_esmfold2.py` and `patches/modal_esmfold2.patch` — this wrapper has no MSA path at all |
| params | `--num-loops 10 --num-sampling-steps 68 --num-diffusion-samples 1`, seed 42, L40S |
| construct | 1 binder : 3 protomers, 555 residues, 35 complexes, 35/35 folded |

## What the calibration's `ef2full` row was

`reference/anthropic-binder-design-protocol.md`: *"ESMFold2-Full uses target-chain MSAs with
msa_max_seq=2048; binder single-sequence."*

So the released campaign's `ef2full` column — the one carrying **pae_interface_min AUC 0.901**,
which is why §4.5 selected that metric — was measured **with target-chain MSAs**. Our run has none.

## The consequence, stated precisely

1. **Our configuration is neither of the two calibrated arms.** It is the **Full checkpoint in
   single-sequence mode**. `ef2full` in the calibration is Full + target MSAs; `ef2fast` is a
   *different checkpoint* (`biohub/ESMFold2-Fast`), selected by env var, not by MSA presence.
   An earlier reading of mine that "no MSA ⇒ this is the Fast arm" was **wrong**: the two arms are
   distinguished by checkpoint, and the notes say so at `reference/wrapper-notes.md:149`.
2. **So AUC 0.901 does not attach to these numbers.** It was measured on a configuration this
   repository cannot currently produce, because no MSA staging exists for ESMFold2 here. (MSA
   staging *does* exist for Boltz — `reference/wrapper-notes.md:265` — and is described there as
   reusable for any model that needs alignments.)
3. **The numbers are still internally valid and still comparable to one thing: this project's own
   12 ipSAE fixtures**, which were produced by the same MSA-less wrapper. On those,
   `pae_interface_min` reads **0.332 for barnase/barstar** and **19.137 for the shuffled negative**,
   with renumber- and chain-swap-invariance at 0.332. That is a same-config positive and negative,
   and it is the only calibration these 35 values currently sit inside.

## What would close this properly, cheapest first

| | option | cost | what it buys |
|---|---|---|---|
| **A** | Score the project's own on-construct controls through the identical path — `targets/validation-tnf/pos_tnf_tnfr2.faa` plus the four `neg_shuffled_tnfr2_*.faa` | **~5 folds, ~6 min, ~$0.25** | a same-config, same-target, same-stoichiometry positive and negative band. Makes the 35 values interpretable **without** any cross-run transfer |
| B | Add MSA staging to the ESMFold2 wrapper and re-run | hours + ~$1.30 | reproduces the protocol's `ef2full` exactly, so AUC 0.901 transfers |
| C | Report the values with the non-comparability stated and no AUC claim | $0 | honest, but leaves the ranker switch justified only by external data on a config we did not run |

**Recommendation: A.** It is the project's own blocked dependency (`arms-backlog.md` 3a: a positive
control per target) applied to this configuration, it costs about a quarter of a dollar, and the
fixtures already exist and are construct-audited. B is the only option that makes the 0.901 claim
literally true for our numbers, and it is not a four-days-out job.

## What must change in the documents either way

- **METHODS §4.5** currently implies the selected metric's AUC was measured on our arm. It was
  measured on the protocol's MSA-fed Full arm. The sentence needs the qualifier.
- **Limitation 51** says the ranker is "declared and plumbed, not yet applied." It is now *applied*,
  but on a configuration whose discriminative power is unmeasured. That is a different, and more
  specific, limitation than the one currently written.
- **The methodology box, FIELD 2 point 4** cites AUC 0.901 "on our own arm". That phrase is not
  supportable as written and should say which configuration the figure comes from.

## Provenance

Run `ap-6IaWZdHjrGore0gL6EhZAO`, 35/35 folds, mean 54.5 s/fold, ~$1.30. Outputs in
`runs/p2-esmfold2-out/p2rank/`, 105 files also persisted to the `esmfold2-runs` Volume — the first
GPU exercise of that per-fold Volume write, which `challenges/02-tnf-alpha.md` flagged as untested.
Scored by `bin/pae_interface.py --dir runs/p2-esmfold2-out/p2rank`.
