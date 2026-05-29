# Container-side wrapper around utils/PAGES2k_reducedProxy_UNSC.R.
#
# Runs the dimensionality-reduction step that turns a PAGES2k-style proxy
# matrix (one column per proxy, one row per year) into the single reduced-
# proxy series RPind.csv that BayGMST_v1.0.R consumes when ptype is not
# "ALL_cached_Barboza".
#
# The upstream reducer uses here::here() relative to a project root that
# contains config.yml + data/. We stage that layout under /app/data and
# write the bridge config.yml the reducer reads.

suppressPackageStartupMessages({
  library(yaml)
  library(here)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

user_config_path <- Sys.getenv("BAYGMST_CONFIG", unset = "/app/config/user_config.yml")
refdata_dir      <- Sys.getenv("BAYGMST_REFDATA", unset = "/app/reference_data")

if (!file.exists(user_config_path)) {
  stop(sprintf("user config not found at %s", user_config_path))
}

user_cfg <- yaml::read_yaml(user_config_path)

setwd("/app")
dir.create("/app/data", showWarnings = FALSE)

# Stage the proxy matrix + metadata under data/ where the reducer expects
# them (here('data', 'PAGES2K_proxy_matrix_screened_1900-2000.csv') etc.).
for (fname in c(
  "PAGES2K_proxy_matrix_screened_1900-2000.csv",
  "PAGES2K_proxy_metadata_screened_1900-2000.csv"
)) {
  src <- file.path(refdata_dir, fname)
  dst <- file.path("/app/data", fname)
  if (!file.exists(src)) {
    stop(sprintf("Required PAGES2k input missing: %s. The LiPD adapter must run before the reducer.", src))
  }
  file.copy(src, dst, overwrite = TRUE)
}

# The reducer also reads cfg$folder_paths$instr_temp_path via here().
file.copy(
  file.path(refdata_dir, "HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv"),
  "/app/data/HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv",
  overwrite = TRUE
)

output_dir <- Sys.getenv("BAYGMST_OUTPUT", unset = "/results")
figures_dir <- file.path(output_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

bridge_cfg <- list(
  ptype     = user_cfg$ptype     %||% "ALL",
  rp_method = toupper(user_cfg$rp_method %||% "PCR"),
  partition_years = list(
    t1 = user_cfg$partition_years$t1 %||% 1,
    t2 = user_cfg$partition_years$t2 %||% 1850,
    t3 = user_cfg$partition_years$t3 %||% 2000
  ),
  folder_paths = list(
    instr_temp_path = "data/HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv",
    # The reducer writes RP_ts.png to figures_dir; without this set, png()
    # gets an empty filename and the upstream script halts.
    figures_dir = figures_dir
  )
)
yaml::write_yaml(bridge_cfg, "/app/config.yml")

message("[run_reducer.R] running upstream reducer...")
source("/app/utils/PAGES2k_reducedProxy_UNSC.R", echo = FALSE)

# Upstream writes RPind.csv to the working dir or data/ — normalize to /app
# so run_baygmst.R can pick it up.
for (cand in c("/app/RPind.csv", "/app/data/RPind.csv", "data/RPind.csv", "RPind.csv")) {
  if (file.exists(cand)) {
    if (cand != "/app/RPind.csv") file.copy(cand, "/app/RPind.csv", overwrite = TRUE)
    message(sprintf("[run_reducer.R] RPind.csv located at %s", cand))
    break
  }
}

if (!file.exists("/app/RPind.csv")) {
  stop("Reducer finished but RPind.csv was not produced.")
}

message("[run_reducer.R] done.")
