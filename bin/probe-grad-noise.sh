#!/usr/bin/env bash
# Problem 2, s45: the GRADIENT-NOISE probe.
#
# Tests the one axis all sixteen trajectories held fixed. Weights, step budget and footprint were
# each varied and each came back flat; the ESTIMATOR never was. Every run before 2026-10-07 drew
# ONE sample per gradient evaluation and integrated it under momentum 0.9, on a measured per-step
# interface SD of 0.0345 against a signal of ~0.13 -- a mechanism for a flat curve on its own.
#
# This is NOT ArcRefine's structural carryover. Carryover needs a parent that already binds and 0
# of our 25 clear the gate. Both axes varied here are first-class MIT API in the commit the image
# already pins (b94b9d4): build_multisample_loss re-runs structure+confidence from ONE trunk
# output and averages, so 4 samples cost well under 4x. No parent, no warm start, no fork, no
# patent licence.
#
# MATCHED TO CONDITION 2 -- the best of the five -- so the estimator is the only variable:
# pinned 9-position epitope, anchor 27 on chain B, w_acid 0.5, lengths 76/84, seeds 0/1,
# steps_soft 38 / steps_sharp 12, species leg on. Identical to the p2probe2-* runs in every
# respect except --grad-samples and --momentum-soft.
#
# BAR, PRE-REGISTERED in s45 BEFORE this was ever run: four-run mean rise in iptm_repred across
# the soft phase, read by analysis/02-tnf/loss_traj.py --block 13. The free-footprint probe gave
# +0.004 +/- 0.013. A rise worth acting on is >= +0.12. Under +0.02 ELIMINATES cause 2 and ends
# the search -- it is not grounds for a seventh condition.
#
# MEMORY is the engineering risk, not cost. The 4x vmap is memory-hungry by its own upstream
# docstring and an L40S already cannot hold three protomers -- which is why the MONOMER leg stays
# single-sample (it carries no interface or pH term) and why the fallback is --grad-samples 2.
# SMOKE FIRST. The free-footprint path had never been executed and was broken in two places.
#
# usage: bin/probe-grad-noise.sh <tag> <length> <seed>
#        SMOKE=1 bin/probe-grad-noise.sh smoke 76 0     # 3+2 steps, prove the path first
set -euo pipefail
: "${1:?usage: probe-grad-noise.sh <tag> <length> <seed>}"
TAG="$1"; LEN="${2:-76}"; SEED="${3:-0}"

cd "$(dirname "$0")/.."
python3 bin/mosaic_selftest.py >/dev/null
echo "[probe] selftest PASS"

# An array, not a string. An unquoted $COMMON collapsed into one giant option twice in this
# project; zsh does not word-split, which is the opposite of the bug and just as fatal.
ARGS=(
  --step 4
  --run-name "p2grad-$TAG"
  --out-dir ../runs/mosaic-p2
  --target ../targets/tnf/tnf_trimer_renum.pdb
  --target-chain A,B,C
  --anchor-chain B
  --crop none
  --anchor 27
  # PASSED EXPLICITLY, both legs. The smoke run died here: DEFAULT_EPITOPE in the wrapper is
  # problem 1's EGFR site ("403,...,409" = TKQHGQF, anchor H433), so omitting --epitope on a TNF
  # run does not mean "the pinned 9 positions", it means "design against the wrong protein".
  # It failed closed -- STEP 4 refused because 403-409 are not in the TNF chain -- but nothing
  # local would have caught it, and a launcher that relies on this default is a trap for the
  # next one too. These are condition 2's values, read back from runs/mosaic-p2/p2probe-a.
  --epitope 16,27,28,70,72,81,82,85,86
  --target2-epitope 13,24,25,66,68,77,78,81,82
  --mechanism his_near_cation
  --his-d0 6.5
  --target2 ../targets/tnf/tnf_mouse_trimer_renum.pdb
  --target2-chain A,B
  --target2-anchor 24
  --w-species 1.0
  --w-acid 0.5                   # the s6.2 ratio, the best of four conditions
  --lengths "$LEN"
  --n-seeds 1
  --seed0 "$SEED"
  --steps-soft 38                # condition 2's split exactly, not the free probe's 75/25
  --steps-sharp 12
  --grad-samples 4               # THE VARIABLE. 1 everywhere before 2026-10-07
  --momentum-soft 0.0            # THE OTHER VARIABLE. 0.9 everywhere before 2026-10-07
)
# NOTE: this is the matched arm to CONDITION 2 (runs p2probe-a..d: w_acid 0.5 with the binding
# weights at 1.0), NOT to condition 3 (p2probe2-e..h, which raised iptm and contact to 3.0) and
# NOT to the free-footprint probe. Every other value here is design()'s default and already
# agrees with p2probe-a/config.json: iptm 1.0, contact 1.0, pae 0.05, within 1.0, plddt 1.0,
# glob 0.3, mpnn 3.0, comp 10.0, monomer 0.5, caps V.12/G.08/H.08, sharpen 1.3,
# sampling_steps 10, recycling_steps 1, stepsize 0.0 (-> 0.1*sqrt(L)).

cd biomodals
GPU="${GPU:-L40S}" TIMEOUT="${TIMEOUT:-120}" \
modal run ${DETACH:+--detach} modal_mosaic.py "${ARGS[@]}" ${SMOKE:+--smoke}
