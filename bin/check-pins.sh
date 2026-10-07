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

# A second kind of pin, added 2026-10-07: the Mosaic install must carry a revision. Unpinned it
# resolved to whatever that repo's HEAD was on the day a reader ran it, which quietly widened the
# METHODS 13 reproducibility claim. biomodals/ is gitignored, so the PATCH is what a reader applies
# -- both must carry the pin or the fix reaches nobody.
MOSAIC_REV=b94b9d4eb9907a700a6d78ed2d29d3704c5df46c
for f in modal_mosaic.py ../patches/modal_mosaic.patch; do
  if grep -q "escalante-bio/mosaic@${MOSAIC_REV}" "$f" 2>/dev/null; then
    echo "ok: $f  (mosaic @ ${MOSAIC_REV:0:7})"
  else
    echo "UNPINNED (reproducibility): $f  -> install must read escalante-bio/mosaic@${MOSAIC_REV}"; BAD=1
  fi
done
exit $BAD
