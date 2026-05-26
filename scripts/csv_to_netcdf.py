#!/usr/bin/env python3
"""Convert BayGMST's gmst_reconstruction_data.csv to a 1D NetCDF.

The BayGMST R script writes a long-form CSV with both reconstruction and
instrumental rows. We emit a CF-friendly NetCDF with:

  dims:
      time = (n_years,)
  variables:
      time(time)                int     years AD
      gmst_mean(time)           float   posterior mean GMST anomaly (deg C)
      gmst_lo_68(time)          float   16th percentile (68% lower)
      gmst_hi_68(time)          float   84th percentile (68% upper)
      gmst_lo_95(time)          float   2.5th percentile (95% lower)
      gmst_hi_95(time)          float   97.5th percentile (95% upper)
      gmst_obs(time)            float   instrumental observation (NaN outside instrumental period)
      type(time)                S16     'reconstruction' or 'instrumental'

This sidecar format is what downstream visualizers (presto-viz or a static
plot generator) can read; the original CSV is still committed for human
inspection.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def convert(in_csv: Path, out_nc: Path) -> None:
    df = pd.read_csv(in_csv)
    df = df.sort_values("year").reset_index(drop=True)

    def col(name: str, default=np.nan) -> np.ndarray:
        if name in df.columns:
            return df[name].to_numpy(dtype=float)
        return np.full(len(df), default, dtype=float)

    types = df["type"].astype(str).to_numpy() if "type" in df.columns \
            else np.full(len(df), "reconstruction", dtype=object)

    ds = xr.Dataset(
        data_vars={
            "gmst_mean":  ("time", col("T.mean")),
            "gmst_lo_68": ("time", col("T.lo.68CrI")),
            "gmst_hi_68": ("time", col("T.hi.68CrI")),
            "gmst_lo_95": ("time", col("T.lolo.95CrI")),
            "gmst_hi_95": ("time", col("T.hihi.95CrI")),
            "gmst_obs":   ("time", col("T.obs")),
            "type":       ("time", np.array(types, dtype="S16")),
        },
        coords={"time": df["year"].to_numpy(dtype=int)},
        attrs={
            "title": "Bayesian Global Mean Surface Temperature reconstruction (BayGMST)",
            "source": "BayGMST v1.0 (Bagwell et al. 2025)",
            "Conventions": "CF-1.10",
        },
    )

    ds["time"].attrs.update(units="years_AD", standard_name="time")
    for v in ("gmst_mean", "gmst_lo_68", "gmst_hi_68", "gmst_lo_95", "gmst_hi_95", "gmst_obs"):
        ds[v].attrs.update(units="K", long_name="GMST anomaly (deg C)")

    out_nc.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_nc)
    print(f"[csv_to_netcdf] wrote {out_nc} ({len(df)} time steps)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", required=True, type=Path)
    ap.add_argument("--out-nc", required=True, type=Path)
    args = ap.parse_args()
    convert(args.in_csv, args.out_nc)


if __name__ == "__main__":
    main()
