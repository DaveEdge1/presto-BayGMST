# PReSto BayGMST Template

A [PReSto](https://paleopresto.com) (Paleoclimate Reconstruction Storehouse) template for the **Bayesian Global Mean Surface Temperature** reconstruction of Bagwell et al. (2025), packaged for automated execution via GitHub Actions.

PReSto lowers the barrier to running, reproducing, and customizing paleoclimate reconstructions. This repository wraps the R/Stan code from [paleopresto/BayGMST_R](https://github.com/paleopresto/BayGMST_R) into a workflow that pulls proxy data from [LiPDverse](https://lipdverse.org), fits the Bayesian state-space model, and publishes the reconstruction to GitHub Pages.

## Method

BayGMST is a Bayesian hierarchical model that infers Global Mean Surface Temperature (GMST) anomalies over the past millennium from a *reduced proxy* — a 1D composite derived from many spatially distributed proxy records via PCR / LASSO / sPCR / SIR / SPLS — together with external climate forcings (CO₂, volcanic, solar) and instrumental HadCRUT temperatures. AR(1) structure on both temperature and proxy equations handles serial autocorrelation; missing pre-instrumental temperatures are estimated as model parameters via HMC sampling in Stan (`cmdstanr`).

Proxy observations are drawn from either:
- **Cached Barboza et al. (2019) reduced proxies** baked into the image at `reference_data/barboza_rps/` (set `ptype: ALL_cached_Barboza` in `config/user_config.yml` — fast, deterministic, ignores LiPDverse).
- **LiPDverse queries** issued through PReSto's interactive map (or `query_params.json`), converted to BayGMST's proxy-matrix format by `scripts/lipd_to_baygmst.py`, then reduced by `utils/PAGES2k_reducedProxy_UNSC.R`.

The original BayGMST_R code is at [paleopresto/BayGMST_R](https://github.com/paleopresto/BayGMST_R).

## File Structure

| Path | Purpose |
|------|---------|
| `config/user_config.yml`       | Reconstruction parameters (overwritten per run by PReSto) |
| `query_params.json`            | LiPDverse query filters (committed by PReSto to trigger the workflow) |
| `Dockerfile` + `entrypoint.sh` | R + Stan + Python environment baked for CI; runs the full pipeline |
| `stan/`                        | Stan source for the BayGMST model (`BayGMST_v1.0.stan`) |
| `R_scripts/BayGMST_v1.0.R`     | Upstream Bayesian fit script |
| `utils/PAGES2k_reducedProxy_UNSC.R` | Dimensionality reducer (proxy matrix → 1D reduced proxy) |
| `scripts/run_baygmst.R`        | Container-side wrapper that bridges `user_config.yml` → upstream paths |
| `scripts/run_reducer.R`        | Container-side wrapper around the R reducer |
| `scripts/lipd_to_baygmst.py`   | LiPD legacy pickle → PAGES2k-format proxy matrix CSVs |
| `scripts/csv_to_netcdf.py`     | Reconstruction CSV → 1D CF-NetCDF for downstream viz |
| `reference_data/`              | Forcings, instrumental HadCRUT, cached Barboza RPs, PAGES2k reference matrix |
| `.github/workflows/baygmst.yml`| Two-job CI pipeline (data prep + reconstruct) |
| `.github/workflows/visualize.yml` | Builds a GitHub Pages site from reconstruction artifacts |

## Workflows

### `baygmst.yml` — BayGMST Reconstruction

Two-job pipeline triggered by a push to `query_params.json` or manual dispatch:

1. **prepare-data** — Reads `ptype` from `config/user_config.yml`. If `ptype` is `ALL_cached_Barboza`, this job no-ops (the cached RPs are already baked into the image). Otherwise it acquires proxy data:
   - *Archived*: downloads a pre-built compilation pickle from LiPDverse
   - *Filtered*: runs the `lipdGenerator` Docker container to query LiPDverse and produce `lipd_legacy.pkl`

   Either way, citation metadata is generated (merged into the existing `CITATION.cff`) and the pickle is uploaded as an artifact.

2. **reconstruct** — Builds the BayGMST Docker image (`rocker/r-ver` + cmdstan + R/Python deps), then runs `entrypoint.sh`:
   - LiPD pickle → proxy-matrix CSVs (`scripts/lipd_to_baygmst.py`) — skipped for cached path
   - Proxy matrix → `RPind.csv` via the R reducer — skipped for cached path
   - Stan fit (`scripts/run_baygmst.R` → `R_scripts/BayGMST_v1.0.R`) writes the reconstruction CSV + figures
   - CSV → 1D NetCDF (`scripts/csv_to_netcdf.py`)

   Results are committed to `results/` on the main branch and uploaded as a 90-day artifact.

### `visualize.yml` — Visualization

Triggered after a successful `baygmst.yml` run (or manually). Downloads the reconstruction artifact and publishes the R-generated PNGs + reconstruction CSV/NetCDF to a static GitHub Pages site with the standard PReSto tile UI.

> **Note:** BayGMST is a 1D global mean reconstruction, not a spatial field, so this template does **not** use the `presto-viz` reusable workflow that LMR2 / Holocene DA use. The static GitHub Pages site renders the R-generated time-series figure directly. To wire BayGMST into `presto-viz`, the NetCDF emitted by `scripts/csv_to_netcdf.py` (1D `gmst_*` variables on a `time` axis) would need a compatible visualizer template.

The repo's About-section URL is populated by the PReSto webhook handler using the submitting user's OAuth token (the default workflow `GITHUB_TOKEN` lacks admin scope to update repo settings, so the template doesn't try). For repos created outside PReSto, set the URL once via the ⚙ next to "About".

## How to Use

### Via PReSto (intended path)

1. Submit a reconstruction request through [paleopresto.com](https://paleopresto.com); pick BayGMST as the method and tune the knobs (proxy selection, time window, sampler iterations).
2. PReSto forks this repository, overwrites `config/user_config.yml`, and commits `query_params.json`. The push triggers `baygmst.yml`.
3. Results land in `results/` on your fork's `main` branch within ~30–60 min (depending on `iter_sampling`).
4. The reconstruction is visualized at `https://<your-username>.github.io/<your-fork-name>/`.

### Locally (for development)

```sh
# Build the container (first build takes ~15 min for cmdstan + R packages)
docker build -t baygmst:local .

# Run with the bundled cached Barboza RPs (ptype: ALL_cached_Barboza)
mkdir -p results
docker run --rm \
  -v $(pwd)/config/user_config.yml:/app/config/user_config.yml:ro \
  -v $(pwd)/results:/results \
  baygmst:local

# Or run with a LiPD pickle (after setting ptype != ALL_cached_Barboza in user_config.yml)
docker run --rm \
  -v $(pwd)/lipd_legacy.pkl:/proxies/lipd_legacy.pkl:ro \
  -v $(pwd)/config/user_config.yml:/app/config/user_config.yml:ro \
  -v $(pwd)/results:/results \
  baygmst:local
```

Outputs land in `results/figures/*.png`, `results/reconstructions/gmst_reconstruction_data.csv`, and `results/reconstructions/gmst_reconstruction.nc`.

## Citing

If you use this reconstruction, please cite both the platform and the underlying method (see [`CITATION.cff`](CITATION.cff)).
