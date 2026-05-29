# Container-side wrapper around R_scripts/BayGMST_v1.0.R.
#
# The upstream BayGMST_R script expects:
#   - a config.yml at the working directory
#   - relative paths to data/* and outputs/*
#   - cmdstan installed at cfg$folder_paths$cmdstan_path
#
# In the container, paths and the cmdstan install differ, so this wrapper
# (a) loads the user-supplied user_config.yml, (b) overrides folder_paths
# to point at the baked-in reference_data and the mounted /results dir,
# and (c) writes a synthesized config.yml the upstream script can read.

suppressPackageStartupMessages({
  library(yaml)
  library(cmdstanr)
})

user_config_path <- Sys.getenv("BAYGMST_CONFIG", unset = "/app/config/user_config.yml")
refdata_dir      <- Sys.getenv("BAYGMST_REFDATA", unset = "/app/reference_data")
output_dir       <- Sys.getenv("BAYGMST_OUTPUT", unset = "/results")

if (!file.exists(user_config_path)) {
  stop(sprintf("user config not found at %s", user_config_path))
}

user_cfg <- yaml::read_yaml(user_config_path)

# Tiny helper for "use user value if set, else fallback default".
`%||%` <- function(a, b) if (is.null(a)) b else a

# Build the folder_paths block expected by the upstream R script, plus
# fold in cmdstan's actual install location.
cmdstan_path_val <- cmdstanr::cmdstan_path()

# Normalize rp_method to uppercase at the entry point so all downstream code
# (load_proxies match.arg, RPind.csv staging) sees the canonical form.
# Coerce SPCR -> PCR: sPCR is valid in the reducer but the cached-RP loader
# only accepts LASSO/PCR/SIR/SPLS, and the Barboza directory has no SPCR CSV.
rp_method_norm <- toupper(user_cfg$rp_method %||% "PCR")
if (rp_method_norm == "SPCR") rp_method_norm <- "PCR"

cfg <- list(
  rp_method = rp_method_norm,
  ptype     = user_cfg$ptype     %||% "ALL_cached_Barboza",
  stan_params = list(
    iter_warmup   = user_cfg$stan_params$iter_warmup   %||% 500,
    iter_sampling = user_cfg$stan_params$iter_sampling %||% 1500
  ),
  partition_years = list(
    t1 = user_cfg$partition_years$t1 %||% 1001,
    t2 = user_cfg$partition_years$t2 %||% 1850,
    t3 = user_cfg$partition_years$t3 %||% 2000
  ),
  folder_paths = list(
    cmdstan_path        = cmdstan_path_val,
    stan_code_path      = "/app/stan/BayGMST_v1.0.stan",
    forcings_path       = file.path(refdata_dir, "forcing.csv"),
    instr_temp_path     = file.path(refdata_dir, "HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv"),
    barboza_rps_path    = file.path(refdata_dir, "barboza_rps"),
    output_dir          = output_dir,
    figures_dir         = file.path(output_dir, "figures"),
    reconstruction_dir  = file.path(output_dir, "reconstructions")
  ),
  co2_params = list(
    c0       = user_cfg$co2_params$c0       %||% 280,
    co2_coef = user_cfg$co2_params$co2_coef %||% 5.35
  ),
  vol_params = list(
    vol_coef = user_cfg$vol_params$vol_coef %||% 25.0
  )
)

dir.create(cfg$folder_paths$figures_dir,        recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$folder_paths$reconstruction_dir, recursive = TRUE, showWarnings = FALSE)

# Route a LiPD-derived RPind.csv (produced upstream by run_reducer.R) through
# the *cached-Barboza* code path in BayGMST_v1.0.R. We do this instead of
# letting the upstream script source the reducer inline because the upstream
# reducer starts with `rm(list=ls())`, which would clobber the parent
# script's cfg and crash the run. By copying RPind.csv into the barboza_rps
# directory under the user's chosen rp_method name and pinning ptype to
# ALL_cached_Barboza, the upstream load_proxies() picks up the live RP via
# its safe, well-tested cache loader.
if (file.exists("/app/RPind.csv")) {
  target <- file.path(cfg$folder_paths$barboza_rps_path,
                      sprintf("RP_new_All_%s.csv", cfg$rp_method))
  dir.create(dirname(target), showWarnings = FALSE, recursive = TRUE)
  file.copy("/app/RPind.csv", target, overwrite = TRUE)
  message(sprintf("[run_baygmst.R] LiPD-derived RP staged at %s", target))
  cfg$ptype <- "ALL_cached_Barboza"
}

# The upstream script does `yaml::read_yaml("config.yml")` from the working
# directory, so write the synthesized config there and cd into /app before
# sourcing it.
setwd("/app")
yaml::write_yaml(cfg, "config.yml")

message("[run_baygmst.R] effective ptype=", cfg$ptype, ", rp_method=", cfg$rp_method)
message("[run_baygmst.R] folder_paths:")
str(cfg$folder_paths)

source("/app/R_scripts/BayGMST_v1.0.R", echo = FALSE)
message("[run_baygmst.R] done.")
