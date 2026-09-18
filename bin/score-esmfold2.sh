#!/usr/bin/env bash
# Score with ESMFold2 using Anthropic's PUBLISHED protocol parameters.
# The wrapper defaults (loops=3, steps=50) are weaker than the protocol (loops=10, steps=68).
# Defaults to 5 seeds (the protocol's FINAL tier). The wrapper is patched to loop
# seeds inside ONE container, so 5 seeds costs ~1.16x a single seed, not 5x --
# model load is ~137s and a fold ~4s. Override with SEEDS=42 for a screen tier.
# Note: this script reads MODAL_GPU / MODAL_TIMEOUT, NOT GPU / TIMEOUT. See reference/wrapper-notes.md.
set -euo pipefail
: "${1:?usage: score-esmfold2.sh <input.faa> [run_name]}"
IN="$1"; RUN="${2:-$(date +%y%m%d%H%M)}"; SEEDS="${SEEDS:-1,2,3,4,5}"
cd "$(dirname "$0")/../biomodals"
MODAL_GPU="${MODAL_GPU:-L40S}" MODAL_TIMEOUT="${MODAL_TIMEOUT:-30}" \
modal run modal_esmfold2.py \
  --input-faa "$IN" \
  --seed "$SEEDS" \
  --num-loops 10 \
  --num-sampling-steps 68 \
  --num-diffusion-samples 1 \
  --out-dir ../runs/esmfold2 \
  --run-name "$RUN"
