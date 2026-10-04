#!/usr/bin/env bash
# Score with ESMFold2 using Anthropic's PUBLISHED protocol parameters.
# The wrapper defaults (loops=3, steps=50) are weaker than the protocol (loops=10, steps=68).
# Defaults to 5 seeds (the protocol's FINAL tier). The wrapper is patched to loop
# seeds inside ONE container, so 5 seeds costs ~1.16x a single seed, not 5x --
# model load is ~137s and a fold ~4s. Override with SEEDS=42 for a screen tier.
# Note: this script reads MODAL_GPU / MODAL_TIMEOUT, NOT GPU / TIMEOUT. See reference/wrapper-notes.md.
#
# Three bugs cost real money on 2026-10-02; all three are handled below:
#   1. This script cd's into biomodals/, so a path relative to the repo root
#      resolved wrong and the run died instantly. Input is now made absolute first.
#   2. MODAL_TIMEOUT is in MINUTES and defaulted to 30. Every prior wave ran
#      10 folds (~28 min) and cleared it by ~2 minutes. A 20-fold run hit the
#      1800s wall at fold 14 and -- because ESMFold2 returns results only at the
#      END -- lost all 13 completed folds. The timeout is now sized from the
#      actual workload, and the estimate is printed before launch.
#   3. --out-dir was hardcoded, so every run landed in the same tree regardless
#      of intent. It is now overridable via OUT_DIR.
set -euo pipefail
: "${1:?usage: score-esmfold2.sh <input.faa|dir-of-faa> [run_name]}"

IN_RAW="$1"; RUN="${2:-$(date +%y%m%d%H%M)}"; SEEDS="${SEEDS:-1,2,3,4,5}"

# --- footgun 1: resolve the input BEFORE cd'ing, so repo-relative paths work ---
if [[ ! -e "$IN_RAW" ]]; then
  echo "error: input not found: $IN_RAW" >&2
  echo "  (paths are resolved from your current directory, not from biomodals/)" >&2
  exit 1
fi
IN="$(cd "$(dirname "$IN_RAW")" && pwd)/$(basename "$IN_RAW")"

# --- footgun 2: size the timeout from the workload instead of hoping ---
if [[ -d "$IN" ]]; then
  N_COMPLEX=$(find "$IN" -maxdepth 1 -name '*.faa' | wc -l | tr -d ' ')
  [[ "$N_COMPLEX" -gt 0 ]] || { echo "error: no .faa files in $IN" >&2; exit 1; }
else
  N_COMPLEX=1
fi
N_SEEDS=$(awk -F, '{print NF}' <<<"$SEEDS")
N_FOLDS=$(( N_COMPLEX * N_SEEDS ))
# Fold time scales with complex SIZE, and the flat 160s/fold constant this used to
# carry is calibrated on the 620+150 residue full-ECD complexes. On the 170-residue
# domain III construct, folds measured 6-14s -- so the estimate overshot by ~13x and
# REFUSED TO LAUNCH a run that needed 55 minutes (2026-10-04, 8 apps x 254 folds).
# That is the third time this constant has misled a launch decision, so it is now
# derived from the actual residue count in the input.
#   measured anchors: 770 residues ~ 140s/fold ; 240 residues ~ 10s/fold
# => roughly quadratic in length; a linear fit through those two points is enough
#    for a timeout, and it is rounded UP.
N_RES=$(awk '/^[A-Z]/ {n+=length($0)} END {print n}' "$(find "$IN" -maxdepth 1 -name '*.faa' | head -1)" 2>/dev/null)
N_RES=${N_RES:-770}
SEC_PER_FOLD=$(( 10 + (N_RES - 240) * 130 / 530 ))
(( SEC_PER_FOLD < 8 )) && SEC_PER_FOLD=8
EST_MIN=$(( (N_FOLDS * SEC_PER_FOLD + 180) / 60 + 1 ))
echo "[score-esmfold2] ${N_RES} residues/complex -> ~${SEC_PER_FOLD}s per fold"
NEED_MIN=$(( EST_MIN * 2 ))                     # 2x headroom for a slow GPU draw
MODAL_TIMEOUT="${MODAL_TIMEOUT:-$NEED_MIN}"

echo "[score-esmfold2] $N_COMPLEX complex(es) x $N_SEEDS seed(s) = $N_FOLDS folds"
echo "[score-esmfold2] estimate ~${EST_MIN} min; MODAL_TIMEOUT=${MODAL_TIMEOUT} min"
if (( MODAL_TIMEOUT < EST_MIN )); then
  echo "[score-esmfold2] REFUSING TO LAUNCH: timeout ${MODAL_TIMEOUT}min is below the" >&2
  echo "  ~${EST_MIN}min estimate. ESMFold2 writes output only at the END, so a timeout" >&2
  echo "  loses every completed fold and you pay for all of it. Raise MODAL_TIMEOUT" >&2
  echo "  (minutes) or split the input into smaller batches." >&2
  exit 1
fi

cd "$(dirname "$0")/../biomodals"
# --- footgun 3: out-dir is now overridable ---
MODAL_GPU="${MODAL_GPU:-L40S}" MODAL_TIMEOUT="$MODAL_TIMEOUT" \
modal run modal_esmfold2.py \
  --input-faa "$IN" \
  --seed "$SEEDS" \
  --num-loops 10 \
  --num-sampling-steps 68 \
  --num-diffusion-samples 1 \
  --out-dir "${OUT_DIR:-../runs/esmfold2}" \
  --run-name "$RUN"
