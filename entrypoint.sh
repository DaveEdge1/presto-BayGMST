#!/usr/bin/env bash
# BayGMST container entrypoint.
#
# Pipeline:
#   1. (optional) Convert LiPD pickle → PAGES2k-style proxy matrix CSVs
#      [skipped when ptype == ALL_cached_Barboza, since cached RPs are used]
#   2. (optional) Run the R reducer to produce RPind.csv from that matrix
#      [skipped when ptype == ALL_cached_Barboza]
#   3. Run BayGMST_v1.0.R to fit the Stan model and write outputs
#   4. Convert the CSV reconstruction to a 1D NetCDF for downstream viz
#
# All paths are controlled by /app/config/user_config.yml plus the env vars
# baked in by the Dockerfile. CI provides:
#   -v <pkl>:/proxies/lipd_legacy.pkl:ro
#   -v <user_config.yml>:/app/config/user_config.yml:ro
#   -v <output_dir>:/results

set -euo pipefail

CONFIG="${BAYGMST_CONFIG:-/app/config/user_config.yml}"
REFDATA="${BAYGMST_REFDATA:-/app/reference_data}"
OUT="${BAYGMST_OUTPUT:-/results}"
PKL="${LIPD_PICKLE:-/proxies/lipd_legacy.pkl}"

mkdir -p "$OUT/figures" "$OUT/reconstructions"

# Parse ptype out of the user config to decide whether to skip steps 1 and 2.
# Using Python (already installed) so we don't have to ship yq.
PTYPE=$(python3 -c "
import yaml
p = yaml.safe_load(open('$CONFIG')).get('ptype', 'ALL_cached_Barboza')
print(p if isinstance(p, str) else ','.join(str(x) for x in p))
")
echo "[entrypoint] ptype = $PTYPE"

if [ "$PTYPE" != "ALL_cached_Barboza" ]; then
    if [ ! -f "$PKL" ]; then
        echo "[entrypoint] ERROR: ptype=$PTYPE requires a LiPD pickle at $PKL, but none found." >&2
        exit 1
    fi
    echo "[entrypoint] Step 1: LiPD pickle → PAGES2k-style proxy matrix"
    python3 /app/scripts/lipd_to_baygmst.py \
        --pickle "$PKL" \
        --out-matrix   "$REFDATA/PAGES2K_proxy_matrix_screened_1900-2000.csv" \
        --out-metadata "$REFDATA/PAGES2K_proxy_metadata_screened_1900-2000.csv"

    echo "[entrypoint] Step 2: R reducer → RPind.csv"
    # The reducer reads config.yml via here::here(), so we hand it one at
    # /app/config.yml that points at the (mounted) user_config.yml values.
    Rscript /app/scripts/run_reducer.R
fi

echo "[entrypoint] Step 3: Stan / BayGMST fit"
Rscript /app/scripts/run_baygmst.R

echo "[entrypoint] Step 4: CSV → 1D NetCDF for visualization"
python3 /app/scripts/csv_to_netcdf.py \
    --in-csv  "$OUT/reconstructions/gmst_reconstruction_data.csv" \
    --out-nc  "$OUT/reconstructions/gmst_reconstruction.nc"

echo "[entrypoint] Done. Results in $OUT:"
ls -lhR "$OUT" || true
