#!/usr/bin/env bash
# Generate binders with BoltzGen. Bounded by default so a hung job can't run for two hours.
set -euo pipefail
: "${1:?usage: design-boltzgen.sh <spec.yaml> [num_designs] [run_name]}"
SPEC="$1"; N="${2:-10}"; RUN="${3:-$(date +%y%m%d%H%M)}"
cd "$(dirname "$0")/../biomodals"
GPU="${GPU:-L40S}" TIMEOUT="${TIMEOUT:-30}" \
modal run modal_boltzgen.py \
  --input-yaml "$SPEC" \
  --protocol "${PROTOCOL:-protein-anything}" \
  --num-designs "$N" \
  --out-dir ../runs/boltzgen \
  --run-name "$RUN"
