# BayGMST through the BayGMST R package (CRAN 0.1.0), replacing the
# config.yml-driven R_scripts/BayGMST_v1.0.R + utils/ reducer pair.
#
# Reads the same user_config.yml and the same inputs, and writes the same
# outputs the rest of the template expects:
#   <out>/reconstructions/gmst_reconstruction_data.csv   (legacy column names,
#       read by scripts/csv_to_netcdf.py and the viz)
#   <out>/reconstructions/fit_post_summaries.csv
#   <out>/reconstructions/reduced_proxy.csv
#   <out>/figures/reconstruction_ts.png   plot(fit): the manuscript figure
#   <out>/figures/trace_plots.png
#
# Proxy input, by ptype:
#   ALL_cached_Barboza  the cached Barboza et al. (2019) reduced proxy
#   anything else       the proxy matrix scripts/lipd_to_baygmst.py wrote,
#                       reduced here with reduce_proxies()

suppressPackageStartupMessages({
  library(yaml)
  library(BayGMST)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

cfg     <- yaml::read_yaml(Sys.getenv("BAYGMST_CONFIG", unset = "/app/config/user_config.yml"))
refdata <- Sys.getenv("BAYGMST_REFDATA", unset = "/app/reference_data")
out     <- Sys.getenv("BAYGMST_OUTPUT", unset = "/results")
dir.create(file.path(out, "reconstructions"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

t1 <- cfg$partition_years$t1 %||% 1001
t2 <- cfg$partition_years$t2 %||% 1850
t3 <- cfg$partition_years$t3 %||% 2000
stopifnot(t1 <= t2, t2 <= t3)
years <- t1:t3

# The package spells supervised PCR "sPCR"; the config has always taken any case.
method <- toupper(cfg$rp_method %||% "PCR")
method <- c(PCR = "PCR", LASSO = "LASSO", SPCR = "sPCR", SPLS = "SPLS", SIR = "SIR")[[method]]
ptype  <- cfg$ptype %||% "ALL"
ptype  <- if (is.character(ptype)) ptype else paste(unlist(ptype), collapse = ",")

# forcing.csv has one row per year from 1 CE and no year column.
forc <- utils::read.csv(file.path(refdata, "forcing.csv"))
forc$year <- seq_len(nrow(forc))
inst <- utils::read.csv(file.path(refdata, "HadCRUT.5.1.0.0.analysis.summary_series.global.annual.csv"))
colnames(inst)[1:2] <- c("year", "T")
inst <- inst[inst$year >= t2 & inst$year <= t3, ]

if (identical(ptype, "ALL_cached_Barboza")) {
  rp_file <- file.path(refdata, "barboza_rps", sprintf("RP_new_All_%s.csv", method))
  rp_df <- utils::read.csv(rp_file)          # columns Year, RP1
  proxy <- as_baygmst_proxy(rp_df)
  rp_out <- data.frame(year = rp_df$Year, RP1 = rp_df$RP1)
  message("[run_package.R] cached Barboza reduced proxy: ", basename(rp_file))
} else {
  pm <- utils::read.csv(file.path(refdata, "PAGES2K_proxy_matrix_screened_1900-2000.csv"), check.names = FALSE)
  pm <- pm[pm$year %in% years, , drop = FALSE]
  message(sprintf("[run_package.R] proxy matrix: %d years x %d records", nrow(pm), ncol(pm) - 1))
  proxy <- reduce_proxies(
    proxy_matrix = pm[, -1, drop = FALSE],
    years        = pm$year,
    temp_calib   = inst$T,
    calib_years  = inst$year,
    method       = method,
    chunk        = cfg$reducer$chunk %||% 250,
    max_na_frac  = cfg$reducer$max_na_frac %||% 0.05
  )
  print(proxy)
  rp_out <- proxy$composite                  # columns year, RP1
}
utils::write.csv(rp_out, file.path(out, "reconstructions", "reduced_proxy.csv"), row.names = FALSE)

iF <- match(years, forc$year)
if (anyNA(iF)) stop("forcing.csv does not cover ", t1, "-", t3)
f <- transform_forcings(
  G = forc$CO2[iF], V = forc$volcanic[iF], S = forc$solar[iF],
  co2_c0   = cfg$co2_params$c0       %||% 280,
  co2_coef = cfg$co2_params$co2_coef %||% 5.35,
  vol_coef = cfg$vol_params$vol_coef %||% 25
)

sp <- cfg$stan_params %||% list()
fit <- fit_baygmst(
  proxy          = proxy,
  instrumental_T = inst$T[match(years, inst$year)],
  forcing_G      = f$G, forcing_V = f$V, forcing_S = f$S,
  years          = years,
  chains          = sp$chains          %||% 4,
  parallel_chains = sp$parallel_chains %||% 2,
  iter_warmup     = sp$iter_warmup     %||% 500,
  iter_sampling   = sp$iter_sampling   %||% 1500,
  seed            = cfg$seed %||% 1
)
print(summary(fit))

# Legacy column names: csv_to_netcdf.py and the viz read these.
rec <- reconstruct(fit)
legacy <- data.frame(
  year = rec$year, T.obs = rec$T_obs, T.mean = rec$T_mean,
  T.lo.68CrI = rec$T_lo_inner, T.hi.68CrI = rec$T_hi_inner,
  T.lolo.95CrI = rec$T_lo_outer, T.hihi.95CrI = rec$T_hi_outer,
  type = rec$type)
utils::write.csv(legacy, file.path(out, "reconstructions", "gmst_reconstruction_data.csv"), row.names = FALSE)

post <- fit$fit$summary(variables = c("alpha0", "alpha1", "phi_R", "phi_T", "beta0",
                                      "betaG", "betaS", "betaV", "sigma_y", "sigma_z"))
utils::write.csv(as.data.frame(post), file.path(out, "reconstructions", "fit_post_summaries.csv"), row.names = FALSE)

ggplot2::ggsave(file.path(out, "figures", "reconstruction_ts.png"), plot(fit),
                width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(file.path(out, "figures", "trace_plots.png"), plot_trace(fit),
                width = 10, height = 7, units = "in", dpi = 200, bg = "white")
message("[run_package.R] done.")
