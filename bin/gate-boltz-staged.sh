#!/usr/bin/env bash
# Boltz gate using a STAGED target MSA -- zero calls to the public MSA server.
# The target's alignment was generated once (--use_msa_server) and saved; every
# design against that target reuses the file.
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
FAILED=0
for f in targets/validation_boltz_staged/neg_*.yaml; do
  name=$(basename "$f" .yaml)
  [ -d "runs/gate/boltz_msa/$name" ] && { echo "skip: $name"; continue; }
  echo ">>> boltz-staged / $name"
  ( cd biomodals && GPU=L40S TIMEOUT=25 modal run modal_boltz.py \
      --input-yaml "../$f" \
      --msa-file "../targets/msa/barnase.csv" \
      --params-str "--write_full_pae --seed 42 --recycling_steps 10 --diffusion_samples 1 --output_format mmcif" \
      --out-dir "../runs/gate/boltz_msa" --run-name "$name" ) \
    || { echo "FAILED: $name"; FAILED=1; }
done
echo; echo "====== BOLTZ (target MSA) — GATE ======"
for d in runs/gate/boltz_msa/*/; do
  name=$(basename "$d")
  pae=$(find "$d" -name "pae_*_model_0.npz" | head -1)
  cif=$(find "$d" -name "*_model_0.cif" | head -1)
  [ -n "$pae" ] && printf "%-26s %s\n" "$name" \
    "$("$ROOT/.venv/bin/python" bin/ipsae_min.py "$pae" "$cif" 2>/dev/null | grep -oE 'ipSAE_min=[0-9.]+' | head -1)"
done
exit $FAILED
