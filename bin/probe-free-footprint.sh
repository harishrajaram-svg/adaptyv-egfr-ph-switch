#!/usr/bin/env bash
# Problem 2, s41 family B: the FREE-FOOTPRINT probe.
#
# Tests the one remaining testable cause of the flat interface trajectory (METHODS s6.3):
# that the nine pinned epitope positions admit no gradient path to a bound pose. Everything
# else is held identical to the s6.2 probe -- same weights, same lengths, same seeds, same
# target and anchor -- so the footprint restriction is the only variable.
#
# 100 steps, not the 50 the pinned probe used, for one reason: 50 steps cannot distinguish
# "the free footprint does not help" from "it never had time to find a site", and a flat
# result is the likely outcome. The soft phase does not anneal against its own length
# (stepsize = 0.1*sqrt(L), momentum 0.9, scale 1.0 are all constants), so a 100-step run
# CONTAINS the matched 50-step comparison as a prefix. Measured detectable rise over 4 runs:
# 0.017 at 50 steps, 0.012 at 100. Neither is the binding constraint -- a rise worth acting
# on is ~+0.12 -- so the extra steps buy interpretability, not sensitivity.
#
# BOTH legs must be free. --epitope none alone leaves the species leg pinned to 9 mouse
# positions, and the symmetry check then refuses the run. The ANCHOR stays pinned in both:
# the pH mechanism needs its one cation contact, which is what makes this "free footprint,
# pinned anchor" rather than Boyd's unconstrained problem.
#
# usage: bin/probe-free-footprint.sh <tag> <length> <seed>
#        SMOKE=1 bin/probe-free-footprint.sh smoke 76 0     # 3+2 steps, prove the path first
set -euo pipefail
: "${1:?usage: probe-free-footprint.sh <tag> <length> <seed>}"
TAG="$1"; LEN="${2:-76}"; SEED="${3:-0}"

cd "$(dirname "$0")/.."
python3 bin/mosaic_selftest.py >/dev/null
echo "[probe] selftest PASS"

# An array, not a string. An unquoted $COMMON collapsed into one giant option twice in this
# project; zsh does not word-split, which is the opposite of the bug and just as fatal.
ARGS=(
  --step 4
  --run-name "p2free-$TAG"
  --out-dir ../runs/mosaic-p2
  --target ../targets/tnf/tnf_trimer_renum.pdb
  --target-chain A,B,C
  --anchor-chain B
  --crop none
  --epitope none                 # FREE on human
  --anchor 27
  --mechanism his_near_cation
  --his-d0 6.5
  --target2 ../targets/tnf/tnf_mouse_trimer_renum.pdb
  --target2-chain A,B
  --target2-anchor 24
  --target2-epitope none         # FREE on mouse -- required, see above
  --w-species 1.0
  --w-acid 0.5                   # the s6.2 ratio, the best of four conditions
  --lengths "$LEN"
  --n-seeds 1
  --seed0 "$SEED"
  --steps-soft 75
  --steps-sharp 25
)

cd biomodals
GPU="${GPU:-L40S}" TIMEOUT="${TIMEOUT:-120}" \
modal run ${DETACH:+--detach} modal_mosaic.py "${ARGS[@]}" ${SMOKE:+--smoke}
