# ============================================================
# BayGMST_R: Bayesian global temperature reconstruction pipeline 
#
# This script:
#   1. Reads user configuration from config.yml
#   2. Loads reduced proxies, forcings, and instrumental temperatures
#   3. Aligns all inputs onto a common annual time grid
#   4. Applies forcing transformations / normalization
#   5. Prepares the data list required by the Stan model
#   6. Fits the Bayesian hierarchical model with CmdStan
#   7. Saves posterior summaries
#   8. Produces a reconstruction figure and posterior histograms
# ============================================================

library(config)
library(yaml)
library(cmdstanr)
library(car)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

cfg <- yaml::read_yaml("config.yml")
#str(cfg)
#cfg$rp_method

### SET cmdstan PATH
# CmdStan is the command line interface to Stan
set_cmdstan_path(path = cfg$folder_paths$cmdstan_path)

### FUNCTIONS:
# ------------------------------------------------------------
# Helper function: load the reduced-proxy dataset specified
# in the config file. The proxy method determines which CSV
# is loaded from the reduced-proxy directory.
# ------------------------------------------------------------
load_proxies <- function(rp_method = c("LASSO", "PCR", "SIR", "SPLS"),
                         base_dir = cfg$folder_paths$barboza_rps_path) {
  rp_method <- match.arg(toupper(rp_method), c("LASSO", "PCR", "SIR", "SPLS"))
  f <- file.path(base_dir, sprintf("RP_new_All_%s.csv", rp_method))
  stopifnot(file.exists(f))
  read.csv(f)
}

### INPUTS
# ------------------------------------------------------------
# Load all model inputs:
#   - reduced proxies
#   - external forcings
#   - instrumental temperature observations
#
# Also pull the reconstruction / calibration window from the
# config file.
# ------------------------------------------------------------
rp_method        <- cfg$rp_method
Proxies.in       <- load_proxies(rp_method)
Forcings.in      <- read.csv(cfg$folder_paths$forcings_path)
Forcings.in$year <- as.integer(rownames(Forcings.in))
Temperatures.in  <- read.csv(cfg$folder_paths$instr_temp_path)
colnames(Temperatures.in) <- c("year","T","l95","u95")

t1 <- cfg$partition_years$t1
t2 <- cfg$partition_years$t2
t3 <- cfg$partition_years$t3
t4 <- cfg$partition_years$t4

# validate ordering: t1 <= t2 <= t3
# ------------------------------------------------------------
# Validate that the reconstruction window is well defined.
# Required ordering is:
#   t1 = start of full reconstruction window
#   t2 = start of instrumental period
#   t3 = end of analysis window
#   t4 = end of future projections window
# with t1 <= t2 <= t3 <= 4.
# ------------------------------------------------------------
check_year <- function(x, name) {
  if (is.null(x)) {
    stop(sprintf("%s is missing.", name), call. = FALSE)
  }
  if (
    !is.numeric(x) ||
    length(x) != 1 ||
    is.na(x) ||
    !is.finite(x) ||
    x <= 0 ||
    x != as.integer(x)
  ) {
    stop(sprintf("%s must be a single positive integer.", name), call. = FALSE)
  }
}
check_year(t1, "t1")
check_year(t2, "t2")
check_year(t3, "t3")

if (!(t1 <= t2 && t2 <= t3)) {
  stop(sprintf("Invalid partition years: require t1 <= t2 <= t3, got t1=%s, t2=%s, t3=%s", t1, t2, t3))
}
if (!is.null(t4)){
  check_year(t4, "t4")
  if (anyNA(c(t4))) {
    stop("t4 contains NA.")
  }
  if (!is.numeric(t4)) {
    stop("t4 must be numeric.")
  }
  if (!(t3 < t4 && t4 <= 2100)) {
    stop(sprintf("Invalid partition years: require t3 < t4 <= 2100, got t3=%s, t4=%s", t3, t4))
  }
}

instru_year_min <- min(Temperatures.in$year, na.rm = TRUE)
instru_year_max <- max(Temperatures.in$year, na.rm = TRUE)
if (!all(c(t2, t3) >= instru_year_min & c(t2, t3) <= instru_year_max)) {
  stop(
    sprintf(
      "t2 and/or t3 are outside the observed instrumental temperature year range [%s, %s]. t2 = %s, t3 = %s",
      instru_year_min, instru_year_max, t2, t3
    )
  )
}

vol_coef <- cfg$vol_params$vol_coef
co2_coef <- cfg$co2_params$co2_coef
co2_c0   <- cfg$co2_params$c0

### BUILD DATAFRAME
# ------------------------------------------------------------
# Build the main analysis dataframe by aligning all inputs
# onto a common annual sequence from t1 to t3.
#
# The resulting dataframe contains:
#   year = annual time index
#   S    = solar forcing
#   V    = volcanic forcing
#   G    = greenhouse gas forcing
#   T    = instrumental temperature anomaly
#   R    = reduced proxy series
# ------------------------------------------------------------
years <- t1:t3

stopifnot(
  !anyDuplicated(Forcings.in$year),
  !anyDuplicated(Temperatures.in$year),
  !anyDuplicated(Proxies.in$Year)
)

Temps_inst <- subset(Temperatures.in, year >= t2)
iF <- match(years, Forcings.in$year)
iT <- match(years, Temps_inst$year)
iP <- match(years, Proxies.in$Year)

df <- data.frame(
  year = years,
  S = Forcings.in$solar[iF],      # solar forcing
  V = Forcings.in$volcanic[iF],   # volcanism forcing
  G = Forcings.in$CO2[iF],        # greenhouse gas (CO2) forcing
  T = Temps_inst$T[iT],           # will be NA for years not in Temperatures.in (e.g., < t2)
  R = Proxies.in$RP1[iP]          # reduced proxy
)

## quick missingness check by column
colSums(is.na(df))



# REPLACE LATER: add in projections for CO2 and solar
Forcings.projections <- read.csv('/Users/tylerbagwell/Documents/GitHub/BayGMST_R/data/forcings_with_prediction_HanWang.csv')
tail(Forcings.projections)
df_prj <- subset(Forcings.projections, year %in% seq(2001,2100,1))
head(Forcings.projections)



### NORMALIZE
# ------------------------------------------------------------
# Apply the forcing transformations used by the model:
#   - volcanic forcing is transformed to a negative saturating form
#   - CO2 forcing is log-transformed relative to a baseline
#   - solar forcing is centered
# ------------------------------------------------------------
df$V <- -abs(vol_coef)*(1-exp(-df$V))
df$G <- co2_coef*log(df$G/co2_c0)
df$S <- df$S - mean(df$S)

df_prj$solar <- df_prj$solar - mean(df_prj$solar)
df_prj$CO2_RCP_2.6 <- co2_coef*log(df_prj$CO2_RCP_2.6/co2_c0)
df_prj$CO2_RCP_4.5 <- co2_coef*log(df_prj$CO2_RCP_4.5/co2_c0)
df_prj$CO2_RCP_6.0 <- co2_coef*log(df_prj$CO2_RCP_6.0/co2_c0)
df_prj$CO2_RCP_8.5 <- co2_coef*log(df_prj$CO2_RCP_8.5/co2_c0)


### PREPARE DATA FOR STAN (VARIBALES BELOW ARE CONSISTENT WITH STAN CODE)
# ------------------------------------------------------------
# Prepare the observed and missing temperature indices for Stan.
#
# In this setup:
#   y = instrumental temperature series (partially observed)
#   z = reduced proxy series
#
# Stan receives both the observed temperature values and the
# index locations of observed vs. missing entries, so that
# missing historical temperatures can be estimated.
# ------------------------------------------------------------
iter_warmup   <- cfg$stan_params$iter_warmup
iter_sampling <- cfg$stan_params$iter_sampling

if (iter_sampling < 1000) {
  stop("cfg$stan_params$iter_sampling must be at least 1000.")
}

nfold <- 3
nobs_fold <- floor((t3-t2)/nfold)
colors <- c('deepskyblue', 'red', 'lawngreen')#, 'gold2', 'magenta1')


r2_cv_vec <- c()
df_leftout <- matrix(
  NA_real_,
  nrow = (t3 - t2),
  ncol = 8
) 
colnames(df_leftout) <- c('year', 'T.mean', 'T.lolo', 'T.lo', 'T.hi', 'T.hihi', 'fold', 'color')
df_leftout <- as.data.frame(df_leftout)

for (i in 1:nfold){
  y <- df$T
  z <- df$R
  
  # holdout fold
  idx_cv_start <- (t2 - t1 + 1 + (i-1)*nobs_fold)
  idx_cv_end   <- (t2 - t1 + i*nobs_fold)
  idx_cv       <- idx_cv_start:idx_cv_end
  y_cv_true    <- y[idx_cv]
  y[idx_cv]    <- NA
  
  NT      <- length(z)
  idx_obs <- which(!is.na(y))
  idx_mis <- which(is.na(y))
  y_obs   <- as.vector(y[idx_obs])
  
  data_list <- list(
    NT = NT, 
    NT_obs = length(idx_obs), 
    NT_mis = length(idx_mis),
    NT_cv  = length(idx_cv),
    idx_obs = as.integer(idx_obs),
    idx_mis = as.integer(idx_mis),
    idx_cv  = as.integer(idx_cv),
    G = as.vector(df$G),
    S = as.vector(df$S),
    V = as.vector(df$V),
    y_obs = y_obs,
    y_cv_true = y_cv_true,
    z = z
  )
  
  message("Running STAN model now...")
  mod <- cmdstan_model("BayGMST_v1.0_5fcv.stan")
  t <- system.time({
    fit <- mod$sample(
      data = data_list,
      chains = 4,
      parallel_chains = 2,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling
    )
  })
  elapsed_sec <- unname(t["elapsed"])
  elapsed_sec
  message("Done.")
  
  # y_mis draws and summaries
  draws_mean <- fit$draws("y_mis")
  idx_names  <- paste0("y_mis[", seq_along(idx_mis), "]")
  mat        <- posterior::as_draws_matrix(draws_mean)[, idx_names, drop = FALSE]
  y1_post    <- cbind(
    t = idx_mis,
    mean = apply(mat, 2, mean),
    lo = apply(mat, 2, quantile, 0.160),
    hi = apply(mat, 2, quantile, 0.840),
    lolo = apply(mat, 2, quantile, 0.025),
    hihi = apply(mat, 2, quantile, 0.975)
  )
  
  df_pred <- data.frame(
    year = as.numeric(idx_mis + t1),
    T.mean = as.numeric(y1_post[, "mean"]),
    T.lo   = as.numeric(y1_post[, "lo"]),
    T.hi   = as.numeric(y1_post[, "hi"]),
    T.lolo   = as.numeric(y1_post[, "lolo"]),
    T.hihi   = as.numeric(y1_post[, "hihi"])
  )
  
  df_obs <- data.frame(
    year = as.numeric(idx_obs + t1),
    T    = as.numeric(y_obs)
  )
  
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "year"]  <- df$year[idx_cv]
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "T.mean"]  <- df_pred$T.mean[(t2-t1+1):(t2-t1+nobs_fold)]
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "T.lolo"]  <- df_pred$T.lolo[(t2-t1+1):(t2-t1+nobs_fold)]
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "T.lo"]    <- df_pred$T.lo[(t2-t1+1):(t2-t1+nobs_fold)]
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "T.hi"]    <- df_pred$T.hi[(t2-t1+1):(t2-t1+nobs_fold)]
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "T.hihi"]  <- df_pred$T.hihi[(t2-t1+1):(t2-t1+nobs_fold)]
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "color"]   <- colors[i]
  df_leftout[((i - 1) * nobs_fold + 1):(i * nobs_fold), "fold"]    <- i
  
  y.fit <- df_pred$T.mean[(t2-t1+1):(t2-t1+nobs_fold)]
  y.true <- df$T[((t2-t1+1) + (i-1)*nobs_fold):((t2-t1) + (i)*nobs_fold)]
  cor(y.fit, y.true)
  
  ok <- complete.cases(y.fit, y.true)
  ss_res <- sum((y.true[ok] - y.fit[ok])^2)
  ss_tot <- sum((y.true[ok] - mean(y.true[ok]))^2)
  r2 <- 1 - ss_res / ss_tot
  
  r2_cv_vec <- append(r2_cv_vec, r2)
}

r2_cv_vec
median(r2_cv_vec)

head(df_leftout)
tail(df_leftout)

x_lim <- range(c(df_leftout$year), na.rm = TRUE)
y_lim <- range(c(df_leftout$T.lolo, df_leftout$T.hihi), na.rm = TRUE)
r2_cv_label <- sprintf(
  "5-fold CV median R² = %.2f",
  median(r2_cv_vec, na.rm = TRUE)
)
r2_cv_annot <- data.frame(
  x = x_lim[2] - 0.01 * diff(x_lim),
  y = y_lim[1] + 0.06 * diff(y_lim),
  label = r2_cv_label
)

fold_labels <- df_leftout |>
  dplyr::group_by(fold) |>
  dplyr::summarise(
    x = mean(range(year, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::arrange(fold) |>
  dplyr::mutate(
    y = y_lim[2] - 0.04 * diff(y_lim),
    r2 = r2_cv_vec[fold],
    label = sprintf("fold %s\nR² = %.2f", fold, r2)
  )

fold_bounds <- df_leftout |>
  dplyr::group_by(fold) |>
  dplyr::summarise(
    xmin = min(year, na.rm = TRUE),
    xmax = max(year, na.rm = TRUE),
    .groups = "drop"
  )

fold_vlines <- sort(unique(fold_bounds$xmin))
fold_vlines <- fold_vlines[fold_vlines > min(fold_vlines)]

p_ts_cv <- ggplot() +
  geom_vline(
    xintercept = fold_vlines,
    color = "gray35",
    linetype = "dashed",
    linewidth = 0.35,
    alpha = 0.8
  ) +
  geom_text(
    data = fold_labels,
    aes(x = x, y = y, label = label),
    size = 3.2,
    color = "gray20",
    vjust = 1
  ) +
  geom_ribbon(
    data = df_leftout,
    aes(
      x = year,
      ymin = T.lolo,
      ymax = T.hihi,
      group = fold,
      fill = 'gray30'
    ),
    alpha = 0.25
  ) +
  geom_ribbon(
    data = df_leftout,
    aes(
      x = year,
      ymin = T.lo,
      ymax = T.hi,
      group = fold,
      fill = 'gray30'
    ),
    alpha = 0.50
  ) +
  geom_line(
    data = df, aes(x = year, y = T, linetype = "HadCRUT5 (Instrumental Observations)"),
    color = "black", linewidth = 0.8
  ) +
  geom_line(
    data = df_leftout, aes(x = year, y = T.mean, group = fold, color = color), linewidth = 0.8
    ) +
  geom_text(
    data = r2_cv_annot,
    aes(x = x, y = y, label = label),
    hjust = 1,
    vjust = 0,
    size = 3.5,
    color = "gray20"
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_linetype_manual(
    name = "",
    values = c("HadCRUT5 (Instrumental Observations)" = "solid")
  ) +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  labs(x = "year", y = "GMST Anomaly (°C)") +
  theme_light(base_size = 10) +
  theme(legend.position = c(0.01, 0.70),
        legend.justification = c(0, 0),
        legend.background = element_rect(fill = NA, color = NA)
        )

sub_txt <- sprintf(
  "Reconstruction window: (%s, %s);  RP computed via %s;  Fold length: %s years",
  t1, (t2-1), rp_method, nobs_fold
)
p_cv <- p_ts_cv + plot_annotation(
  title = "GMST Reconstructions via 5-folds",
  subtitle = sub_txt
) &
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9)
  )

p_cv

ggsave(paste0(cfg$folder_paths$figures_dir,"/cross_validation_ts.png"), plot = p_cv, width = 7, height = 4, units = "in", dpi = 300, bg = "white")



############## OLD ##############
############## OLD ##############

#fold 1
idx_cv_start <- (t2 - t1 + 1 + 0*nobs_fold)
idx_cv_end   <- (t2 - t1 + 1*nobs_fold)
idx_cv       <- idx_cv_start:idx_cv_end
y_cv_true    <- y[idx_cv]
y[idx_cv]    <- NA

#fold 2
idx_cv_start <- (t2 - t1 + 1 + 1*nobs_fold)
idx_cv_end   <- (t2 - t1 + 2*nobs_fold)
idx_cv       <- idx_cv_start:idx_cv_end
y_cv_true    <- y[idx_cv]
y[idx_cv]    <- NA

#fold 3
idx_cv_start <- (t2 - t1 + 1 + 2*nobs_fold)
idx_cv_end   <- (t2 - t1 + 3*nobs_fold)
idx_cv       <- idx_cv_start:idx_cv_end
y_cv_true    <- y[idx_cv]
y[idx_cv]    <- NA

#fold 4
idx_cv_start <- (t2 - t1 + 1 + 3*nobs_fold)
idx_cv_end   <- (t2 - t1 + 4*nobs_fold)
idx_cv       <- idx_cv_start:idx_cv_end
y_cv_true    <- y[idx_cv]
y[idx_cv]    <- NA

#fold 5
idx_cv_start <- (t2 - t1 + 1 + 4*nobs_fold)
idx_cv_end   <- (t2 - t1 + 5*nobs_fold)
idx_cv       <- idx_cv_start:idx_cv_end
y_cv_true    <- y[idx_cv]
y[idx_cv]    <- NA

#fold 6
idx_cv_start <- (t2 - t1 + 1 + 5*nobs_fold)
idx_cv_end   <- (t2 - t1 + 6*nobs_fold)
idx_cv       <- idx_cv_start:idx_cv_end
y_cv_true    <- y[idx_cv]
y[idx_cv]    <- NA


#
NT      <- length(z)
idx_obs <- which(!is.na(y))
idx_mis <- which(is.na(y))
y_obs   <- as.vector(y[idx_obs])

data_list <- list(
  NT = NT, 
  NT_obs = length(idx_obs), 
  NT_mis = length(idx_mis),
  NT_cv  = length(idx_cv),
  idx_obs = as.integer(idx_obs),
  idx_mis = as.integer(idx_mis),
  idx_cv  = as.integer(idx_cv),
  G = as.vector(df$G),
  S = as.vector(df$S),
  V = as.vector(df$V),
  y_obs = y_obs,
  y_cv_true = y_cv_true,
  z = z
)

### FIT BHM with STAN
# ------------------------------------------------------------
# Fit the Bayesian hierarchical model in Stan.
# The model file is read from the config, compiled, and then
# sampled using HMC through cmdstanr.
# ------------------------------------------------------------





message("Running STAN model now...")
mod <- cmdstan_model("BayGMST_v1.0_5fcv.stan")
t <- system.time({
  fit <- mod$sample(
    data = data_list,
    chains = 4,
    parallel_chains = 2,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling
  )
})
elapsed_sec <- unname(t["elapsed"])
elapsed_sec
message("Done.")


fit$summary(variables = c("r2_cv"))
r2_cv_draws <- fit$draws(variables = "r2_cv")
hist(r2_cv_draws, breaks='scott')


fit$summary(variables = c("mse_cv"))
mse_cv_draws <- fit$draws(variables = "mse_cv")
hist(mse_cv_draws, breaks='scott')


#
draws_mean <- fit$draws("y_mis")
idx_names  <- paste0("y_mis[", seq_along(idx_mis), "]")
mat        <- posterior::as_draws_matrix(draws_mean)[, idx_names, drop = FALSE]
y1_post    <- cbind(
  t = idx_mis,
  mean = apply(mat, 2, mean),
  lo = apply(mat, 2, quantile, 0.160),
  hi = apply(mat, 2, quantile, 0.840),
  lolo = apply(mat, 2, quantile, 0.025),
  hihi = apply(mat, 2, quantile, 0.975)
)

df_pred <- data.frame(
  year = as.numeric(idx_mis + t1),
  T.mean = as.numeric(y1_post[, "mean"]),
  T.lo   = as.numeric(y1_post[, "lo"]),
  T.hi   = as.numeric(y1_post[, "hi"]),
  T.lolo   = as.numeric(y1_post[, "lolo"]),
  T.hihi   = as.numeric(y1_post[, "hihi"])
)

df_obs <- data.frame(
  year = as.numeric(idx_obs + t1),
  T    = as.numeric(y_obs)
)


plot(df_pred$year, df_pred$T.mean, type='l', col='red', lwd=2.0, lty=1, 
     ylim=c(-1.2, +0.70), xlim=c(min(df_pred$year), max(df$year)))
lines(df_pred$year, df_pred$T.lolo,  type='l', col='red', lwd=1.5, lty=5)
lines(df_pred$year, df_pred$T.hihi, type='l', col='red', lwd=1.5, lty=5)
lines(df$year, df$T, Temps_inst$T, type='l', col='black', lwd=1.0, lty=1)


mfold <- 1
y.fit <- df_pred$T.mean[(t2-t1+1):(t2-t1+nobs_fold)]
y.true <- df$T[((t2-t1+1) + (mfold-1)*nobs_fold):((t2-t1) + (mfold)*nobs_fold)]
cor(y.fit, y.true)

ok <- complete.cases(y.fit, y.true)
ss_res <- sum((y.true[ok] - y.fit[ok])^2)
ss_tot <- sum((y.true[ok] - mean(y.true[ok]))^2)
r2 <- 1 - ss_res / ss_tot
r2




