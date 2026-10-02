#!/bin/bash
# Per-arm chain: when this arm's run finishes, relaunch it only if it showed SIGNS OF LIFE.
# usage: chain_arm.sh GEN PREFIX OUTSUB INITIAL_RUN INITIAL_PID [LENGTHS PROTOCOL HOTSPOTS BINDER]
set -u
cd "$HOME/code/adaptyv-2026/biomodals" || exit 1
export PATH="$HOME/.local/bin:$PATH"

GEN=$1; PREFIX=$2; OUTSUB=$3; RUN=$4; PID=$5
LENGTHS=${6:-}; PROTO=${7:-}; HOTS=${8:-}; BINDER=${9:-}

OUT="../runs/$OUTSUB"; ABS="$HOME/code/adaptyv-2026/runs/$OUTSUB"
LOG="$ABS/chain.log"; MAX_RUNS=4; STOP_HOUR=7; TARGET_MIN=120; PESSIMISTIC_MIN=20

log(){ echo "[$(date '+%a %-I:%M:%S %p') $PREFIX] $*" >> "$LOG"; }
# NOTE: `grep -c` prints 0 AND exits 1 on no-match, so `grep -c ... || echo 0` yields "0\n0",
# which breaks `[ "$x" -eq 0 ]` and silently defeated the dead-arm guard. Capture, then default.
_count(){ local n; n=$(grep -c "$1" "$2" 2>/dev/null); echo "${n:-0}"; }
accepted(){
  if [ "$GEN" = "boltzgen" ]; then
    find "$ABS/$1/final_ranked_designs" -name "*.cif" 2>/dev/null | wc -l | tr -d ' '
  else
    ls "$ABS/$1"/Accepted/*.pdb 2>/dev/null | wc -l | tr -d ' '
  fi
}
built(){ [ "$GEN" = "boltzgen" ] && accepted "$1" || _count "Trajectory successfully" "$ABS/$1.log"; }
trajs(){ _count "Starting trajectory: " "$ABS/$1.log"; }

log "watching $RUN (pid $PID)"
while kill -0 "$PID" 2>/dev/null; do sleep 60; done
date +%s > "$ABS/$RUN.end"

N=1
while :; do
  [ -f "$ABS/STOP" ] || [ -f "$HOME/code/adaptyv-2026/runs/STOP_ALL" ] && { log "STOP flag — exiting"; break; }
  A=$(accepted "$RUN"); B=$(built "$RUN"); T=$(trajs "$RUN")
  S=$(cat "$ABS/$RUN.start" 2>/dev/null || echo 0); E=$(cat "$ABS/$RUN.end" 2>/dev/null || echo 0)
  MIN=$(( (E-S)/60 ))
  log "$RUN done: accepted=$A built=$B trajectories=$T elapsed=${MIN}min"

  if [ ! -d "$ABS/$RUN" ]; then log "no output dir — crash. STOPPING this arm."; break; fi
  if [ "$A" -eq 0 ] && [ "$B" -eq 0 ]; then
    log "DEAD: built nothing and accepted nothing. This configuration cannot make a structure. STOPPING."
    break
  fi
  log "SIGN OF LIFE (built=$B accepted=$A) — this arm is worth more compute"

  [ "$N" -ge "$MAX_RUNS" ] && { log "run ceiling. STOPPING."; break; }
  H=$(date +%-H); [ "$H" -ge "$STOP_HOUR" ] && { log "past ${STOP_HOUR} am. STOPPING."; break; }

  if [ "$T" -gt 0 ] && [ "$MIN" -gt 0 ]; then PER=$(( MIN/T )); [ "$PER" -lt 1 ] && PER=1; NEXT=$(( TARGET_MIN/PER ));
  else PER=0; NEXT=15; fi
  [ "$NEXT" -lt 8 ] && NEXT=8; CAP=$(( 160 / PESSIMISTIC_MIN )); [ "$NEXT" -gt "$CAP" ] && NEXT=$CAP   # survive even if every trajectory succeeds (slow)

  N=$((N+1)); RUN=$(printf "%s%02d" "$PREFIX" "$N")
  log "launching $RUN (${NEXT} trajectories, measured ${PER} min/traj)"
  date +%s > "$ABS/$RUN.start"
  if [ "$GEN" = "bindcraft" ]; then
    GPU=H100 TIMEOUT=200 modal run modal_bindcraft.py \
      --input-pdb ../targets/egfr/egfr_d3_6aru.pdb --target-chains A \
      --target-hotspot-residues "$HOTS" --lengths "$LENGTHS" \
      --number-of-final-designs 8 --max-trajectories "$NEXT" \
      --design-protocol "$PROTO" --binder-name "$BINDER" \
      --out-dir "$OUT" --run-name "$RUN" > "$ABS/$RUN.log" 2>&1
  elif [ "$GEN" = "germinal" ]; then
    GPU=H100 TIMEOUT=200 uv run --quiet --with PyYAML modal run modal_germinal.py \
      --target-yaml ../targets/egfr/germinal_egfr_d3.yaml --run-type vhh \
      --max-trajectories 20 --max-passing-designs 8 \
      --experiment-name egfr_d3_vhh --run-name "$RUN" \
      --out-dir "$OUT" > "$ABS/$RUN.log" 2>&1
  else
    GPU=L40S TIMEOUT=150 modal run modal_boltzgen.py \
      --input-yaml ../targets/egfr/boltzgen_egfr_d3.yaml \
      --protocol protein-anything --num-designs 60 \
      --out-dir "$OUT" --run-name "$RUN" > "$ABS/$RUN.log" 2>&1
  fi
  date +%s > "$ABS/$RUN.end"
done
log "ARM CHAIN DONE"
