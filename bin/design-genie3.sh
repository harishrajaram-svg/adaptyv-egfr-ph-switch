#!/usr/bin/env bash
# Genie 3 binder-design GENERATION for problem 2. Phase 0 measures; it does not claim.
#
# WHY GENIE 3. On Anthropic's 150 measured TNF-alpha designs, Genie3 scored 8/23 = 34.8% against
# 4/127 = 3.1% for every other generator combined (Fisher exact, two-sided, p = 0.00003), and it
# is the ONLY generator in that set that produced mouse cross-reactivity (3 of its 8 binders).
# Caveat kept in view: those 8 are one backbone family -- all 83 aa, pairwise identity 0.58-0.81
# inside a set spanning 0.02-0.90 -- so the honest unit is one backbone found in 23 samples.
#
# WHAT PHASE 0 IS FOR. Nothing in the repo or the paper states a generation cost, and
# max_n_chain is 4 while our complex is exactly 4 chains. This measures seconds per sample and
# peak CUDA memory at a ~550-token, 4-chain target. That is the whole deliverable.
#
# TIMEOUT. Unlike bin/score-esmfold2.sh this CANNOT size the timeout from a measured rate,
# because producing that rate is the point. The default is deliberately generous: MODAL_TIMEOUT
# is a ceiling and Modal bills actual use, so an over-long ceiling costs nothing while an
# under-long one loses the run and bills for it anyway.
#
# usage: bin/design-genie3.sh <tag> [n_sample] [cond_strategy]
#        cond_strategy: hotspot | extended | common | region1
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
PY="$ROOT/.venv/bin/python"

: "${1:?usage: design-genie3.sh <tag> [n_sample] [cond_strategy]}"
TAG="$1"; N_SAMPLE="${2:-2}"; STRATEGY="${3:-hotspot}"; BINDER_LEN="${BINDER_LEN:-83}"

# Fail closed on the caller's shell, per the zsh word-splitting collapse that silently turned
# four distinct runs into four identical ones on 2026-10-07.
case "$TAG"       in *[[:space:]]*) echo "[genie3] REFUSING: tag '$TAG' contains whitespace --"\
  " the caller's shell did not split its arguments." >&2; exit 1;; esac
case "$N_SAMPLE"  in ''|*[!0-9]*) echo "[genie3] REFUSING: n_sample '$N_SAMPLE' is not a number" >&2; exit 1;; esac
case "$BINDER_LEN" in ''|*[!0-9]*) echo "[genie3] REFUSING: BINDER_LEN '$BINDER_LEN' is not a number" >&2; exit 1;; esac
case "$STRATEGY" in hotspot|extended|common|region1) ;; *)
  echo "[genie3] REFUSING: cond_strategy '$STRATEGY' is not one of hotspot|extended|common|region1" >&2
  exit 1;; esac

# Pre-flight, the way score-esmfold2.sh runs esmfold2_selftest.py. No GPU, no Modal account.
"$PY" bin/genie3_selftest.py || {
  echo "[genie3] REFUSING TO LAUNCH: genie3_selftest.py failed" >&2; exit 1; }

# The problem JSON is rebuilt every launch, so the Asp143 guard and every pinned residue
# identity are re-checked against the target on disk rather than trusted from a stale file.
PROBLEM="$ROOT/runs/genie3-in/problems/tnfa_corrected.json"
"$PY" bin/build_genie3_problem.py "$PROBLEM" "$BINDER_LEN" || {
  echo "[genie3] REFUSING TO LAUNCH: the problem could not be built" >&2; exit 1; }

MODAL_GPU="${MODAL_GPU:-A100-40GB}"
MODAL_TIMEOUT="${MODAL_TIMEOUT:-30}"
N_RES=$(awk '$1=="ATOM" {print substr($0,22,1) substr($0,23,4)}' targets/tnf/tnf_canonical_trimer.pdb | sort -u | wc -l | tr -d ' ')

echo
echo "[genie3] resolved config -- read it before it spends anything"
echo "[genie3]   tag            $TAG"
echo "[genie3]   n_sample       $N_SAMPLE"
echo "[genie3]   cond_strategy  $STRATEGY"
echo "[genie3]   binder length  $BINDER_LEN   (all 8 measured Genie3 binders were 83 aa)"
echo "[genie3]   target         targets/tnf/tnf_canonical_trimer.pdb, ${N_RES} residues, 3 chains"
echo "[genie3]   tokens         ~$((N_RES + BINDER_LEN)) of max_n_token 2048"
echo "[genie3]   chains         4 of max_n_chain 4  <-- at the architectural ceiling"
echo "[genie3]   GPU            $MODAL_GPU"
echo "[genie3]   timeout        ${MODAL_TIMEOUT} min (a ceiling, not a reservation)"
echo "[genie3]   stage          GENERATE ONLY -- no evaluation, no ColabFold, no MSA server"
echo

cd "$ROOT/biomodals"
MODAL_GPU="$MODAL_GPU" MODAL_TIMEOUT="$MODAL_TIMEOUT" \
modal run modal_genie3.py \
  --problem "$PROBLEM" \
  --n-sample "$N_SAMPLE" \
  --cond-strategy "$STRATEGY" \
  --out-dir "${OUT_DIR:-../runs/genie3-out}" \
  --run-name "$TAG"
