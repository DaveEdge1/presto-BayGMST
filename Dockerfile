# BayGMST_R + cmdstan + a small Python layer for the LiPD → proxy-matrix
# adapter. No pre-built upstream image to extend (BayGMST_R is R/Stan, not
# Python), so we build the full environment here.
#
# The FROM layer + the cmdstan/R-package layer are cached on the runner;
# only the COPY of source code at the end rebuilds when the template's
# scripts/, R_scripts/, stan/, utils/, or config files change.

FROM rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC \
    R_LIBS_USER=/usr/local/lib/R/site-library
# NOTE: do NOT set CMDSTAN here. cmdstanr's .onLoad reads $CMDSTAN at
# package-load time and crashes ("argument is of length zero") if the dir
# exists but is empty — which is the state during the install layer. We
# let cmdstanr fall back to its default ~/.cmdstan/cmdstan-X.Y.Z location;
# scripts/run_baygmst.R retrieves the actual path via cmdstanr::cmdstan_path().

# ── System deps ─────────────────────────────────────────────────────────────
# build-essential / gfortran / cmake : Stan model compilation + RcppParallel
# libcurl/libssl/libxml2/libgit2     : pak / remotes / httr-style R packages
# libnetcdf-dev                      : CSV → 1D NetCDF post-process (ncdf4)
# python3 + pip                      : LiPD pickle adapter (scripts/lipd_to_baygmst.py)
# libxt6                             : silence rocker's font warnings under ggplot2
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran cmake make pkg-config git curl ca-certificates \
        libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
        libfontconfig1-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff5-dev \
        libxt6 libcairo2-dev libnetcdf-dev libhdf5-dev libgeos-dev \
        python3 python3-pip python3-venv \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ── CmdStan ─────────────────────────────────────────────────────────────────
# Install cmdstan via cmdstanr at build time so the toolchain + headers are
# baked into the image. Models still recompile on first sample() call (Stan
# emits per-model .o/.exe), but the cmdstan headers + stanc are pre-built.
RUN R -e "install.packages('cmdstanr', repos = c('https://stan-dev.r-universe.dev', getOption('repos')))" && \
    R -e "cmdstanr::install_cmdstan(cores = 2, overwrite = TRUE)" && \
    R -e "cat('cmdstan path:', cmdstanr::cmdstan_path(), '\n')"

# ── R packages used by BayGMST_v1.0.R + utils/PAGES2k_reducedProxy_UNSC.R ───
# Pinned to CRAN snapshot via rocker/r-ver's RSPM default (reproducible by date).
# Split into two layers to keep the layer size bounded.
RUN R -e "install.packages(c( \
        'yaml','config','here','dplyr','tidyr','ggplot2','patchwork','car', \
        'posterior','reshape2','cowplot','stringr','tibble','readr','jsonlite','ncdf4' \
    ), Ncpus = 2)"

# Dimensionality-reduction packages used by the R reducer. fda + spls + glmnet
# + superpc + dr cover PCR / LASSO / SPLS / SIR. Skip ggmap/maps/maptools —
# they pull in heavy geo deps and the reducer's plots aren't on the critical
# path (we keep them best-effort in the R script).
RUN R -e "install.packages(c( \
        'fda','glmnet','spls','dr','superpc','matrixStats','parallel' \
    ), Ncpus = 2)"

# ── Python deps for the LiPD adapter ────────────────────────────────────────
# scripts/lipd_to_baygmst.py uses the original `lipd` (LiPD-utilities)
# library — same one Holocene DA's da_load_proxies.py uses — because it
# knows how to consume the legacy {'D': {datasetName: {...}}} pickle that
# lipdverse archives. pylipd is the newer RDF-based lib; it doesn't read
# this format. Isolated in a venv so we never fight PEP 668 across Debian
# / Ubuntu base-image versions.
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir --upgrade pip && \
    /opt/venv/bin/pip install --no-cache-dir \
        numpy pandas pyyaml scipy xarray netcdf4 LiPD
ENV PATH=/opt/venv/bin:$PATH

# ── App source ──────────────────────────────────────────────────────────────
# COPY last so cached layers above survive iteration on template code.
# Anything under /app/* is overridden by `docker run -v` mounts in CI for
# config/user_config.yml and the LiPD pickle.
WORKDIR /app
COPY stan/        /app/stan/
COPY R_scripts/   /app/R_scripts/
COPY utils/       /app/utils/
COPY scripts/     /app/scripts/
COPY reference_data/ /app/reference_data/
COPY config/      /app/config/
COPY entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh && touch /app/.here   # marker so here::here() resolves to /app

# Defaults consumed by entrypoint.sh and the R scripts. CI overrides these
# via -e flags / -v mounts; the defaults make `docker run baygmst:local`
# work standalone with the cached Barboza RPs.
ENV BAYGMST_CONFIG=/app/config/user_config.yml \
    BAYGMST_REFDATA=/app/reference_data \
    BAYGMST_OUTPUT=/results \
    LIPD_PICKLE=/proxies/lipd_legacy.pkl

ENTRYPOINT ["/app/entrypoint.sh"]
