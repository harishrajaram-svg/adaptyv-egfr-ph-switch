#!/usr/bin/env bash
# Fetch Anthropic's released binder-design campaign data and derive the TNF-alpha slice.
#
# Why this exists: the campaign shipped 150 TNF-alpha designs with wet-lab binding
# labels, the generator recorded per design, and co-folding metrics from eight models.
# That is the first positive control this project has ever had on its own target --
# 12 measured binders and 138 measured non-binders, same assay vendor (Acro TNA-H4211).
#
#   dataset  huggingface.co/datasets/Anthropic/claude-protein-binder-design
#   licence  CC BY 4.0 (data), MIT (scripts)   verified ungated 2026-10-07
#
# The raw downloads are gitignored (10 MB parquet). The DERIVED TNF-alpha tables are
# committed, so analysis/02-tnf/calibrate_external.py reproduces without the download.
set -euo pipefail

BASE="https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/main/data/tables"
DIR="$(cd "$(dirname "$0")/.." && pwd)/reference/anthropic-campaign"
mkdir -p "$DIR"

fetch() {  # url -> file, skipped when already present and non-empty
  local url="$1" out="$2"
  if [[ -s "$out" ]]; then echo "[have] $(basename "$out")"; return; fi
  echo "[get ] $(basename "$out")"
  curl -fsSL --max-time 300 -o "$out.part" "$url"
  mv "$out.part" "$out"
}

fetch "$BASE/design_summary.csv"                  "$DIR/design_summary.csv"
fetch "$BASE/insilico/cofold_predictions.parquet" "$DIR/cofold_predictions.parquet"

# Derive the TNF-alpha slice as TSV so the analysis stays pure-stdlib (no pyarrow).
# Needs uv only here, once. Fails loudly rather than writing a partial table.
echo "[derive] TNF-alpha cofold slice"
uv run --quiet --with pyarrow python - "$DIR" <<'PY'
import csv, os, sys
import pyarrow.parquet as pq

d = sys.argv[1]
COLS = ["uuid", "full_name", "target", "cofolding_model", "stoichiometry", "seed",
        "ipsae_min", "ipsae_max", "iptm_pae", "plddt_binder", "plddt_target",
        "pae_interface_min", "pae_interface_mean", "sc_dockq", "n_interface_contacts",
        "has_clash"]
rows = [r for r in pq.read_table(os.path.join(d, "cofold_predictions.parquet"),
                                 columns=COLS).to_pylist() if r["target"] == "TNFa"]
if not rows:
    sys.exit("REFUSE: no TNFa rows in cofold_predictions.parquet")
out = os.path.join(d, "tnfa_cofold.tsv")
with open(out + ".part", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
os.replace(out + ".part", out)
print(f"  wrote {len(rows)} rows -> {os.path.basename(out)}")
PY

# Labels + provenance for the same 150 designs, also as a committed TSV.
python3 - "$DIR" <<'PY'
import csv, os, sys
d = sys.argv[1]
KEEP = ["uuid", "full_name", "target", "generator", "sequence_design_method",
        "binder_length", "sequence", "binder_final", "kd_nM_final",
        "mouse_binding_final", "mouse_kd_nM_final", "vendor_agreement",
        "adaptyv_binding", "adaptyv_kd_nM", "adaptyv_expression"]
src = csv.DictReader(open(os.path.join(d, "design_summary.csv")))
rows = [{k: r.get(k, "") for k in KEEP} for r in src if r["target"] == "TNFa"]
if len(rows) != 150:
    sys.exit(f"REFUSE: expected 150 TNFa designs, got {len(rows)}")
out = os.path.join(d, "tnfa_designs.tsv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=KEEP, delimiter="\t")
    w.writeheader(); w.writerows(rows)
print(f"  wrote {len(rows)} rows -> {os.path.basename(out)}")
PY
echo "[done] $DIR"
