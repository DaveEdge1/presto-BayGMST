#!/usr/bin/env python3
"""LiPD legacy pickle -> PAGES2k-style proxy matrix + metadata CSVs.

BayGMST_R's reducer (utils/PAGES2k_reducedProxy_UNSC.R) consumes two
sibling CSVs in data/:

  PAGES2K_proxy_matrix_screened_1900-2000.csv
      header: year, <pid_1>, <pid_2>, ..., <pid_n>
      one row per year (1..2000 by default), NA for missing

  PAGES2K_proxy_metadata_screened_1900-2000.csv
      header: <pid>, lat, lon, elev, ptype
      one row per proxy column in the matrix
      ptype values like "speleothem.d18O", "tree.TRW", "coral.SrCa"

The lipdverse "legacy" pickle (e.g. Pages2kTemperature2_2_0.pkl) wraps a
dict of LiPD-formatted datasets under a top-level 'D' key. We extract the
flat list of time-series records via the original `lipd` (LiPD-utilities)
library — the same one Holocene DA's da_load_proxies.py uses — and then
bin each record onto the common annual year-AD grid.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _coerce_year(years_raw, ages_raw, year_units: str = "") -> np.ndarray:
    """Return year-AD as a float array.

    Prefer the record's 'year' axis if present; otherwise convert 'age'
    (years BP, 1950 reference) — same convention Holocene DA uses.
    """
    if years_raw is not None and len(years_raw):
        arr = np.asarray(years_raw, dtype=float)
        if "bp" in (year_units or "").lower():
            arr = 1950.0 - arr
        return arr
    if ages_raw is not None and len(ages_raw):
        return 1950.0 - np.asarray(ages_raw, dtype=float)
    return np.array([], dtype=float)


def _aggregate_to_annual(years: np.ndarray, vals: np.ndarray,
                         year_axis: np.ndarray) -> np.ndarray:
    """Bin a record onto the common annual axis.

    LiPD records can be irregular / sub-annual / multi-decadal. Bin by
    floor(year) and average within each bin; missing bins stay NaN.
    """
    mask = np.isfinite(years) & np.isfinite(vals)
    if not mask.any():
        return np.full(year_axis.shape, np.nan)
    df = pd.DataFrame({"year": np.floor(years[mask]).astype(int),
                       "v": vals[mask]})
    agg = df.groupby("year", as_index=True)["v"].mean()
    out = pd.Series(np.nan, index=year_axis)
    common = agg.index.intersection(year_axis)
    out.loc[common] = agg.loc[common]
    return out.to_numpy()


def _safe_float(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


# LiPD uses multi-word CamelCase archive types; PAGES2k (and user_config.yml)
# expects short lowercase abbreviations.  Map the compound forms here.
_ARCHIVE_NORM: dict[str, str] = {
    "glacierice":    "ice",
    "lakesediment":  "lake",
    "marinesediment": "marine",
    "wood":          "tree",
    "molluskshell":  "bivalve",
}


def _extract_ts(pkl_path: Path) -> list[dict]:
    """Load the legacy pickle and flatten to a list of time-series dicts.

    Tries lipd.extractTs first (canonical path used by Holocene DA). If
    that library or call fails — e.g. the pickle is shaped slightly
    differently — falls back to a defensive raw walker that handles the
    {'D': {name: ds}} layout directly.
    """
    with pkl_path.open("rb") as f:
        raw = pickle.load(f)

    # Unwrap the 'D' container if present (Pages2kTemperature, Temp12k, etc.)
    D = raw.get("D", raw) if isinstance(raw, dict) else raw

    try:
        import lipd  # type: ignore
        ts = lipd.extractTs(D)
        if ts:
            print(f"[lipd_to_baygmst] lipd.extractTs returned {len(ts)} records", file=sys.stderr)
            return ts
        print("[lipd_to_baygmst] lipd.extractTs returned 0 records — falling back to raw walk",
              file=sys.stderr)
    except Exception as e:
        print(f"[lipd_to_baygmst] lipd.extractTs failed ({e}); falling back to raw walk",
              file=sys.stderr)

    # Raw walker — handles the legacy {'D': {datasetName: <ds>}} layout where
    # each <ds> is itself a dict with paleoData/geo/etc. Some lipdverse pickles
    # nest one extra layer: {'D': [<ds>, <ds>, ...]} (list of ds). Handle both.
    out: list[dict] = []

    def _flatten_ds(name: str, ds: dict) -> None:
        if not isinstance(ds, dict):
            return
        geo_props = ((ds.get("geo") or {}).get("properties") or {})
        lat  = _safe_float(geo_props.get("latitude"))
        lon  = _safe_float(geo_props.get("longitude"))
        elev = _safe_float(geo_props.get("elevation"))
        _raw_archive = str(ds.get("archiveType", "unknown")).lower()
        archive = _ARCHIVE_NORM.get(_raw_archive, _raw_archive)
        for pd_entry in (ds.get("paleoData") or []):
            for mt in (pd_entry.get("measurementTable") or []):
                cols = mt.get("columns") or []
                year_col = next((c for c in cols
                                 if str(c.get("variableName", "")).lower() == "year"), None)
                age_col  = next((c for c in cols
                                 if str(c.get("variableName", "")).lower() == "age"), None)
                year_vals  = year_col.get("values") if year_col else None
                age_vals   = age_col.get("values")  if age_col  else None
                year_units = (year_col or {}).get("units", "")
                for c in cols:
                    if c is year_col or c is age_col:
                        continue
                    vname = str(c.get("variableName", ""))
                    if not vname:
                        continue
                    out.append({
                        "dataSetName": name,
                        "paleoData_variableName": vname,
                        "paleoData_values": c.get("values", []),
                        "year":  year_vals,
                        "age":   age_vals,
                        "yearUnits":  year_units,
                        "geo_meanLat":  lat,
                        "geo_meanLon":  lon,
                        "geo_meanElev": elev,
                        "archiveType":  archive,
                        "paleoData_proxy": c.get("proxy", c.get("proxyGeneral", vname)),
                    })

    if isinstance(D, dict):
        for name, ds in D.items():
            if isinstance(ds, list):
                # {'D': {name: [ds, ds, ...]}} edge case
                for i, sub in enumerate(ds):
                    _flatten_ds(f"{name}__{i}", sub)
            else:
                _flatten_ds(name, ds)
    elif isinstance(D, list):
        for i, ds in enumerate(D):
            _flatten_ds(str(ds.get("dataSetName", f"ds_{i}")) if isinstance(ds, dict) else f"ds_{i}", ds)

    print(f"[lipd_to_baygmst] raw walker returned {len(out)} records", file=sys.stderr)
    return out


def build_csvs(pkl_path: Path, out_matrix: Path, out_metadata: Path,
               year_start: int = 1, year_end: int = 2000) -> None:
    records = _extract_ts(pkl_path)
    if not records:
        raise SystemExit("No records extracted from LiPD pickle — cannot proceed.")

    year_axis = np.arange(year_start, year_end + 1, dtype=int)
    matrix: dict[str, np.ndarray] = {"year": year_axis}
    meta_rows: list[tuple[str, float, float, float, str]] = []

    seen_ids: set[str] = set()
    dropped = 0

    for rec in records:
        ds_name = str(rec.get("dataSetName", "")) or "unknown"
        var_name = str(rec.get("paleoData_variableName", ""))
        # Skip time/depth/metadata columns that aren't proxy series.
        vn_lower = var_name.lower()
        if vn_lower in ("year", "age", "depth", "depthtop", "depthbottom"):
            continue
        # Skip records with no numeric values.
        raw_vals = rec.get("paleoData_values")
        if raw_vals is None:
            continue
        try:
            vals = np.asarray(raw_vals, dtype=float)
        except (TypeError, ValueError):
            dropped += 1
            continue
        years = _coerce_year(rec.get("year"), rec.get("age"),
                             year_units=str(rec.get("yearUnits", "")))
        if years.size == 0 or vals.size == 0 or years.size != vals.size:
            dropped += 1
            continue

        rid = f"{ds_name}__{var_name}"
        base = rid
        i = 1
        while rid in seen_ids:
            i += 1
            rid = f"{base}__{i}"
        seen_ids.add(rid)

        series = _aggregate_to_annual(years, vals, year_axis)
        if not np.isfinite(series).any():
            dropped += 1
            continue

        _raw_archive = str(rec.get("archiveType", "unknown")).lower()
        archive = _ARCHIVE_NORM.get(_raw_archive, _raw_archive)
        proxy   = str(rec.get("paleoData_proxy", var_name)).lower()
        ptype   = f"{archive}.{proxy}" if proxy else archive

        matrix[rid] = series
        meta_rows.append((
            rid,
            _safe_float(rec.get("geo_meanLat")),
            _safe_float(rec.get("geo_meanLon")),
            _safe_float(rec.get("geo_meanElev")),
            ptype,
        ))

    if len(matrix) <= 1:
        raise SystemExit(
            f"All {len(records)} LiPD records were dropped during alignment "
            f"({dropped} explicitly rejected) — cannot proceed."
        )

    df_matrix = pd.DataFrame(matrix)
    df_meta = pd.DataFrame(meta_rows, columns=["", "lat", "lon", "elev", "ptype"])

    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    df_matrix.to_csv(out_matrix, index=False)
    df_meta.to_csv(out_metadata, index=False)
    print(f"[lipd_to_baygmst] wrote {out_matrix} "
          f"({df_matrix.shape[0]} years x {df_matrix.shape[1]-1} proxies)")
    print(f"[lipd_to_baygmst] wrote {out_metadata} ({len(meta_rows)} records, "
          f"{dropped} dropped during alignment)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pickle",       required=True, type=Path, help="Path to lipd_legacy.pkl")
    ap.add_argument("--out-matrix",   required=True, type=Path)
    ap.add_argument("--out-metadata", required=True, type=Path)
    ap.add_argument("--year-start",   type=int, default=1)
    ap.add_argument("--year-end",     type=int, default=2000)
    args = ap.parse_args()

    build_csvs(args.pickle, args.out_matrix, args.out_metadata,
               year_start=args.year_start, year_end=args.year_end)


if __name__ == "__main__":
    main()
