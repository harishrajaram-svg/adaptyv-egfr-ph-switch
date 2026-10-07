#!/usr/bin/env bash
# Problem 2, s47: the EPITOPE probe -- Region I instead of the R108 receptor site.
#
# Sixteen trajectories varied weights, steps, footprint and the gradient estimator on the R108
# site and every one was flat. The epitope was never a variable. Three independent sources say it
# is the suspect: AlphaProteo's stated reason for 0 of 54 ("a flat, highly polar binding site at
# an interface between 2 subunits in a homotrimer" -- our site), Boyd's Nipah post-mortem ("it
# probably comes down to the epitope"), and Glogl et al. measuring that this surface class has
# hydrophobics "separated by distances of up to 28 A -- too small for 65 residue proteins to
# simultaneously engage". And TNF-alpha IS bindable: Chen et al. got 0.55 nM by avoiding the
# receptor site (Commun Biol 2025, 10.1038/s42003-025-09030-7).
#
# REGION I, verified here by analysis/02-tnf/region1_verify.py -- located by MOTIF, not by the
# residue labels the source supplied, because five numbering schemes are live in this project:
#   hotspots  positional 69,70,92 (one protomer) + 109,110 (adjacent) = mature 74/75/97/114/115
#   mouse     65,66,88,105,106 -- 4 of 5 conserved (92 I->V), stable at gap penalties -4..-12
#   seam      2.97 A between 92.B and 110.C, so it genuinely spans two protomers
#   anchor    positional 93 (mature K98), 1.31 A from the patch, K in BOTH species, mouse 89
#   avoid     68 (H->Y) and 133 (R->L) are NOT conserved -- do not anchor on them
#
# MATCHED TO CONDITION 2 with the EPITOPE AS THE ONLY VARIABLE. That includes keeping the OLD
# estimator (grad_samples 1, momentum_soft 0.9): arm A returned a null, so condition 2 is the
# baseline, and changing two things at once would waste the run.
#
# BAR, pre-registered in s47 BEFORE this ran: four-run mean rise in iptm_repred across the soft
# phase, loss_traj.py --block 13. >= +0.12 acts. Under +0.02 is a null and the epitope hypothesis
# dies. The R108 null distribution is seventeen trajectories: +0.004 +/- 0.013 (free footprint),
# +0.008 +/- 0.012 (arm A), and twelve flat tuning trajectories. HARD STOP: this is the last arm.
#
# usage: bin/probe-region1.sh <tag> <length> <seed>
#        SMOKE=1 bin/probe-region1.sh smoke 76 0
set -euo pipefail
: "${1:?usage: probe-region1.sh <tag> <length> <seed>}"
TAG="$1"; LEN="${2:-76}"; SEED="${3:-0}"

# Fail closed on the caller's shell, per the zsh collapse that cost four runs on 2026-10-07.
case "$TAG" in *[[:space:]]*) echo "[probe] REFUSING: tag '$TAG' has whitespace -- the caller's"\
     " shell did not split its arguments, so LEN and SEED silently took defaults." >&2; exit 1;; esac
case "$LEN"  in ''|*[!0-9]*) echo "[probe] REFUSING: length '$LEN' is not a number" >&2; exit 1;; esac
case "$SEED" in ''|*[!0-9]*) echo "[probe] REFUSING: seed '$SEED' is not a number"  >&2; exit 1;; esac
echo "[probe] tag=$TAG L=$LEN seed=$SEED  epitope=RegionI anchor=93"

cd "$(dirname "$0")/.."
python3 bin/mosaic_selftest.py >/dev/null
python3 bin/check_species_map.py >/dev/null
echo "[probe] selftest + species map PASS"

ARGS=(
  --step 4
  --run-name "p2reg1-$TAG"
  --out-dir ../runs/mosaic-p2
  --target ../targets/tnf/tnf_trimer_renum.pdb
  --target-chain A,B,C
  --anchor-chain B
  --crop none
  --mechanism his_near_cation
  --his-d0 6.5
  --target2 ../targets/tnf/tnf_mouse_trimer_renum.pdb
  --target2-chain A,B
  --w-species 1.0
  --w-acid 0.5                       # condition 2's ratio, the best of five
  --lengths "$LEN"
  --n-seeds 1
  --seed0 "$SEED"
  --steps-soft 38                    # condition 2's split exactly
  --steps-sharp 12
  # THE VARIABLE: Region I, both legs, with its own conserved anchor.
  --epitope 69,70,92,109,110
  --anchor 93
  --target2-epitope 65,66,88,105,106
  --target2-anchor 89
)

cd biomodals
GPU="${GPU:-L40S}" TIMEOUT="${TIMEOUT:-120}" \
modal run ${DETACH:+--detach} modal_mosaic.py "${ARGS[@]}" ${SMOKE:+--smoke}
