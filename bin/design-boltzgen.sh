#!/usr/bin/env bash
# Generate binders with BoltzGen. Bounded by default so a hung job can't run for two hours.
set -euo pipefail
: "${1:?usage: design-boltzgen.sh <spec.yaml> [num_designs] [run_name]}"
SPEC="$1"; N="${2:-10}"; RUN="${3:-$(date +%y%m%d%H%M)}"
# OUT_DIR   : where the arm lands (default ../runs/boltzgen)
# DETACH    : set to 1 to pass --detach, so a dead local client (e.g. a network blip)
#             does not kill the Modal app. Two runs died that way on 2026-10-03.
# EXTRA_ARGS: passed through to the boltzgen CLI, e.g. "--diffusion_batch_size 1"
#             (without that, every design in a batch shares ONE sampled length, so a
#              pinned residue lands at the SAME position in all of them)
# This script cd's into biomodals/, so a path relative to the repo root would resolve
# wrong there -- the same footgun that killed a launch on 2026-10-02 (see
# bin/score-esmfold2.sh, bug 1). Resolve the spec to an absolute path FIRST.
if [[ ! -f "$SPEC" ]]; then
  echo "error: spec not found: $SPEC" >&2
  echo "  (resolved from your current directory, not from biomodals/)" >&2
  exit 1
fi
SPEC="$(cd "$(dirname "$SPEC")" && pwd)/$(basename "$SPEC")"

cd "$(dirname "$0")/../biomodals"
GPU="${GPU:-L40S}" TIMEOUT="${TIMEOUT:-30}" \
modal run ${DETACH:+--detach} modal_boltzgen.py \
  --input-yaml "$SPEC" \
  --protocol "${PROTOCOL:-protein-anything}" \
  --num-designs "$N" \
  --out-dir "${OUT_DIR:-../runs/boltzgen}" \
  --run-name "$RUN" \
  ${EXTRA_ARGS:+--extra-args "$EXTRA_ARGS"}
