#!/usr/bin/env bash
# Validation gate: does the instrument separate a real binder from shuffled negatives?
#
# Positive: barnase + barstar (native, KD ~1e-14 M, non-antibody -- co-folders
#           systematically underperform on antibody-antigen, so an antibody
#           control would understate the instrument).
# Negatives: 4 composition-matched shuffles of barstar. Same amino acids,
#           destroyed order -- isolates interface prediction from composition.
#
# PASS = the positive scores clearly above all four negatives.
# If it does not, the instrument has no discriminative power and nothing
# downstream of it means anything.
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
ARMS="${ARMS:-full fast}"
FAILED=0

for arm in $ARMS; do
  case "$arm" in
    full) REPO=biohub/ESMFold2;      REV=1afea82e432079d9af2ebd71d1e4c339ecca2ff0 ;;
    fast) REPO=biohub/ESMFold2-Fast; REV=main ;;
    *) echo "unknown arm: $arm"; exit 2 ;;
  esac
  for f in targets/validation/*.faa; do
    name=$(basename "$f" .faa)
    out="runs/gate/$arm/$name"
    if [ -d "$out" ]; then echo "skip (exists): $arm/$name"; continue; fi
    echo ">>> $arm / $name"
    ( cd biomodals && \
      ESMFOLD2_HF_REPO=$REPO ESMFOLD2_HF_REVISION=$REV \
      MODAL_GPU=L40S MODAL_TIMEOUT=25 \
      modal run modal_esmfold2.py \
        --input-faa "../$f" \
        --num-loops 10 --num-sampling-steps 68 --num-diffusion-samples 1 \
        --out-dir "../runs/gate/$arm" --run-name "$name" ) \
      || { echo "FAILED: $arm/$name"; FAILED=1; }
  done
done

echo
echo "================ GATE RESULTS ================"
for arm in $ARMS; do
  echo "--- arm: $arm ---"
  "$ROOT/.venv/bin/python" bin/ipsae_min.py --dir "runs/gate/$arm" 2>&1 | grep -E "ipSAE_min|FAILED" || echo "  (no scores)"
done
echo
echo "PASS if pos_barnase_barstar > every neg_shuffled_*."
exit $FAILED
