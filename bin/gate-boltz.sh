#!/usr/bin/env bash
# Validation gate, Boltz-2 arm. One YAML per call (this wrapper does not batch),
# so each complex gets its own container. Image is cached, so it's mostly predict time.
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
FAILED=0
for f in targets/validation_boltz/*.yaml; do
  name=$(basename "$f" .yaml)
  out="runs/gate/boltz/$name"
  if [ -d "$out" ]; then echo "skip (exists): $name"; continue; fi
  echo ">>> boltz / $name"
  ( cd biomodals && GPU=L40S TIMEOUT=25 modal run modal_boltz.py \
      --input-yaml "../$f" \
      --params-str "--write_full_pae --seed 42 --recycling_steps 10 --diffusion_samples 1 --output_format mmcif" \
      --out-dir "../runs/gate/boltz" --run-name "$name" ) \
    || { echo "FAILED: $name"; FAILED=1; }
done
echo
echo "=========== BOLTZ ARM — GATE ==========="
for d in runs/gate/boltz/*/; do
  name=$(basename "$d")
  pae=$(find "$d" -name "pae_*_model_0.npz" | head -1)
  cif=$(find "$d" -name "*_model_0.cif" | head -1)
  if [ -n "$pae" ] && [ -n "$cif" ]; then
    printf "%-26s " "$name"
    "$ROOT/.venv/bin/python" bin/ipsae_min.py "$pae" "$cif" 2>/dev/null | grep -oE "ipSAE_min=[0-9.]+" | head -1
  else
    echo "$name  (no pae/cif)"
  fi
done
exit $FAILED
