#!/usr/bin/env bash
# Score with ESMFold2 using Anthropic's PUBLISHED protocol parameters.
# The wrapper defaults (loops=3, steps=50) are weaker than the protocol (loops=10, steps=68).
# Note: this script reads MODAL_GPU / MODAL_TIMEOUT, NOT GPU / TIMEOUT. See reference/wrapper-notes.md.
set -euo pipefail
: "${1:?usage: score-esmfold2.sh <input.faa> [run_name]}"
IN="$1"; RUN="${2:-$(date +%y%m%d%H%M)}"
cd "$(dirname "$0")/../biomodals"
MODAL_GPU="${MODAL_GPU:-L40S}" MODAL_TIMEOUT="${MODAL_TIMEOUT:-30}" \
modal run modal_esmfold2.py \
  --input-faa "$IN" \
  --num-loops 10 \
  --num-sampling-steps 68 \
  --num-diffusion-samples 1 \
  --out-dir ../runs/esmfold2 \
  --run-name "$RUN"
