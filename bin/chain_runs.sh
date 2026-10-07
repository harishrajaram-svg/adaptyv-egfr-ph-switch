#!/bin/bash
# Chain bounded BindCraft runs overnight. Stops on genuine zero-yield; retries a crash once.
set -u
cd "$HOME/code/adaptyv-2026/biomodals" || exit 1
export PATH="$HOME/.local/bin:$PATH"

OUT="../runs/egfr-d3-long"
ABS="$HOME/code/adaptyv-2026/runs/egfr-d3-long"
LOG="$ABS/chain.log"
MAX_RUNS=4          # ceiling on NEW runs launched after run01
STOP_HOUR=7         # never launch a new run at/after 07:00
TARGET_MIN=120      # trajectory-time budget per run; TIMEOUT is 165 (45 min margin)

log(){ echo "[$(date '+%a %H:%M:%S')] $*" >> "$LOG"; }
accepted(){ ls "$ABS/$1"/Accepted/*.pdb 2>/dev/null | wc -l | tr -d ' '; }
trajs(){ grep -c "Starting trajectory: " "$ABS/$1.log" 2>/dev/null || echo 0; }
have_output(){ [ -d "$ABS/$1" ]; }

log "chain started; waiting on run01 (pid $1)"
while kill -0 "$1" 2>/dev/null; do sleep 60; done
echo "$(date +%s)" > "$ABS/run01.end"

PREV=run01
RETRIED=0
N=1

while :; do
  [ -f "$ABS/STOP" ] && { log "STOP file present — exiting"; break; }

  A=$(accepted "$PREV"); T=$(trajs "$PREV")
  S=$(cat "$ABS/$PREV.start" 2>/dev/null || echo 0)
  E=$(cat "$ABS/$PREV.end"   2>/dev/null || echo 0)
  MIN=$(( (E - S) / 60 ))
  RAW=$(find "$ABS/$PREV" -name "*.pdb" 2>/dev/null | wc -l | tr -d " ")
  log "$PREV finished: ${A} accepted designs, ${T} trajectories, ${MIN} min (raw pdb files: ${RAW})"

  if ! have_output "$PREV"; then
    if [ "$RETRIED" -eq 0 ]; then
      log "$PREV produced NO output dir — treating as a crash, retrying once"
      RETRIED=1
    else
      log "second run with no output — infrastructure problem, not a design problem. STOPPING."
      break
    fi
  elif [ "$A" -eq 0 ]; then
    log "ZERO accepted designs with output present. 120-160 aa did not fix it."
    log "DECISION FOR HARISH: pivot to nanobody (Germinal) — reviewer email Q2. Not spending more. STOPPING."
    break
  fi

  # size the next run from measured per-trajectory cost
  if [ "$T" -gt 0 ] && [ "$MIN" -gt 0 ]; then
    PER=$(( MIN / T )); [ "$PER" -lt 1 ] && PER=1
    NEXT=$(( TARGET_MIN / PER ))
  else
    PER=0; NEXT=15
  fi
  [ "$NEXT" -lt 8 ]  && NEXT=8
  [ "$NEXT" -gt 30 ] && NEXT=30

  [ "$N" -ge "$MAX_RUNS" ] && { log "run ceiling ($MAX_RUNS) reached. STOPPING."; break; }
  H=$(date +%H); H=${H#0}
  [ "${H:-0}" -ge "$STOP_HOUR" ] && { log "past ${STOP_HOUR}:00 — not launching another. STOPPING."; break; }

  N=$((N+1)); RUN=$(printf "run%02d" "$((N))")
  log "launching $RUN: ${NEXT} trajectories (measured ${PER} min/trajectory)"
  date +%s > "$ABS/$RUN.start"
  GPU=H100 TIMEOUT=165 modal run modal_bindcraft.py \
    --input-pdb ../targets/egfr/egfr_d3_6aru.pdb --target-chains A \
    --target-hotspot-residues "A409,A411,A412" --lengths "120,160" \
    --number-of-final-designs 8 --max-trajectories "$NEXT" \
    --design-protocol "Beta-sheet" --binder-name egfr_d3_h409_long \
    --out-dir "$OUT" --run-name "$RUN" > "$ABS/$RUN.log" 2>&1
  date +%s > "$ABS/$RUN.end"
  PREV=$RUN
done

TOTAL=0
for d in "$ABS"/run*/; do
  r=$(basename "$d"); n=$(accepted "$r"); TOTAL=$((TOTAL+n))
  log "  $r: $n accepted"
done
log "CHAIN DONE — $TOTAL accepted designs total across all runs"
