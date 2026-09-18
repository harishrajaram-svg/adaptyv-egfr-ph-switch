#!/usr/bin/env bash
# Re-run after any `git pull` in biomodals/. Modal's default Python is 3.14 and breaks numba-based builds.
set -euo pipefail
cd "$(dirname "$0")/../biomodals"
BAD=0
for f in modal_boltzgen.py modal_chai1.py; do
  if grep -q "debian_slim()" "$f" 2>/dev/null; then
    echo "UNPINNED (will fail): $f  -> re-apply python_version=\"3.12\""; BAD=1
  else
    echo "ok: $f"
  fi
done
exit $BAD
