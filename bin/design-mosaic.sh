#!/usr/bin/env bash
# Gradient design with Mosaic (step 4 of biomodals/modal_mosaic.py).
#
# Steps 1-3 are environment probes and already PASSED (2026-10-02, ~$0.45). This runs the
# design. The non-GPU logic has a free self-test -- it runs first, below, every time,
# because a wrong residue index here costs a GPU run, and that has happened.
#
# usage: bin/design-mosaic.sh <run_name> [n_seeds] [lengths]
#
#   SMOKE=1     3+2 steps, 1 seed. Proves the whole path (features -> loss -> APGM ->
#               predict -> re-predict -> geometry -> files) for a few minutes of GPU.
#               RUN THIS FIRST. Mosaic's own README says it needs hand-holding.
#   DETACH=1    pass --detach, so a dead local client does not kill the app. Four runs
#               died to a DNS blip on 2026-10-03. Results are also written to the
#               mosaic-weights Volume under runs/<run_name>/, so a detached run whose
#               client dies is still recoverable -- poll the Volume, not the log.
#   CROP=none   full 609-residue ECD instead of domain III. Expect an OOM on 48 GB.
#   GPU=        default L40S (48 GB). TIMEOUT= minutes, default 90.
#   EXTRA_ARGS= passed through to the entrypoint, e.g. "--w-acid 4.0 --steps-soft 120"
#
# After it lands, in this order and no other (the geometry check caught a 3.59x "switch"
# whose binder was 27 A away):
#   1. geometry   already in designs.tsv: acid_O_to_his_N, geometry_pass. >4.0 A is dead.
#   2. fold+bind  bin/score-esmfold2.sh then bin/ipsae_min.py. ipSAE_min 0.0000 is dead.
#   3. pKa        bin/ph_gate_target_his.py (paired), then bin/knockout_control.py.
#      NOTE: that gate's CANON is BoltzGen-output numbering (H433 -> 409). These files are
#      written in ECD-positional, so it needs CANON = {406: "H433"}. See numbering.json.
set -euo pipefail
: "${1:?usage: design-mosaic.sh <run_name> [n_seeds] [lengths]}"
RUN="$1"; SEEDS="${2:-2}"; LENGTHS="${3:-76}"

cd "$(dirname "$0")/.."
python3 bin/mosaic_selftest.py

cd biomodals
GPU="${GPU:-L40S}" TIMEOUT="${TIMEOUT:-90}" \
modal run ${DETACH:+--detach} modal_mosaic.py \
  --step 4 \
  --run-name "$RUN" \
  --n-seeds "$SEEDS" \
  --lengths "$LENGTHS" \
  --crop "${CROP:-284:453}" \
  --out-dir "${OUT_DIR:-../runs/mosaic}" \
  ${SMOKE:+--smoke} \
  ${EXTRA_ARGS:-}
