#!/usr/bin/env bash
# Pull EVERY BindCraft run off the Modal volume, including Rejected/ and MPNN/.
#
# WHY THIS EXISTS. modal_bindcraft.py mounts a Volume over the whole design path and its
# return value already globs `**/*.*`, so nothing is lost REMOTELY. The loss is local: if the
# client dies (a network blip killed one on 2026-10-03) the return value never arrives, and
# nothing in this repo ever read the Volume afterwards. On 2026-10-04 that turned out to have
# hidden **71 distinct BindCraft sequences** -- every design from d3acid, d3acid2, d3acid3,
# d3acidrelax and bs_full -- while `analysis/01-egfr/bindcraft_accepted.csv` held 6 rows and
# was treated as the complete record. BindCraft clears 0.1493 on both species at ~50% against
# BoltzGen's 3.5%, so those were the highest-prior unscreened designs in the project.
#
# Run this after EVERY BindCraft invocation, and before concluding anything about yield.
#   bin/bindcraft-sync.sh [dest]      default dest: runs/bindcraft-recovered
set -euo pipefail
VOL="${VOL:-bindcraft-out}"
DEST="${1:-runs/bindcraft-recovered}"
mkdir -p "$DEST"
RUNS=$(modal volume ls "$VOL" 2>/dev/null | awk 'NF' || true)
[ -z "$RUNS" ] && { echo "no runs on volume $VOL"; exit 0; }
for r in $RUNS; do
  r="${r%/}"
  echo "[sync] $r"
  mkdir -p "$DEST/$r"
  for f in mpnn_design_stats.csv final_design_stats.csv trajectory_stats.csv failure_csv.csv; do
    modal volume get "$VOL" "$r/$f" "$DEST/$r/$f" --force >/dev/null 2>&1 || true
  done
  for sub in Accepted Rejected MPNN Trajectory; do
    n=$(modal volume ls "$VOL" "$r/$sub" 2>/dev/null | wc -l | tr -d ' ')
    [ "$n" = "0" ] && continue
    echo "    $sub: $n entries"
    modal volume get "$VOL" "$r/$sub" "$DEST/$r/" --force >/dev/null 2>&1 || true
  done
done
echo
echo "[sync] distinct sequences across all mpnn_design_stats.csv:"
python3 - "$DEST" <<'PY'
import csv, glob, sys, os
seqs=set()
for f in glob.glob(os.path.join(sys.argv[1], '*', 'mpnn_design_stats.csv')):
    try:
        for r in csv.DictReader(open(f)):
            s=(r.get('Sequence') or '').strip()
            if len(s) >= 30: seqs.add(s)
    except Exception: pass
print(f"  {len(seqs)} distinct designs recovered -> {sys.argv[1]}")
PY
