#!/usr/bin/env python3
"""LiPD pickle -> PAGES2k-style proxy matrix + metadata CSVs.

BayGMST_R's reducer (utils/PAGES2k_reducedProxy_UNSC.R) consumes two
sibling CSVs in data/:

  PAGES2K_proxy_matrix_screened_1900-2000.csv
      header: year, <pid_1>, <pid_2>, ..., <pid_n>
      one row per year (1..2000 by default), NA for missing

  PAGES2K_proxy_metadata_screened_1900-2000.csv
      header: <pid>, lat, lon, elev, ptype
      one row per proxy column in the matrix
      ptype values like "speleothem.d18O", "tree.TRW", "coral.SrCa"

This adapter reads the LiPD "legacy" pickle produced upstream by
davidedge/lipd_webapps:lipdGenerator and emits those two CSVs. The legacy
pickle is a dict keyed by dataset name, with .paleoData entries holding the
measurement tables; pylipd's LiPD class wraps it and exposes
get_timeseries_essentials().

For BayGMST in particular, every record contributes a SINGLE proxy column
(one-record-one-column), unlike CFR which can split a record into multiple
PSM-calibrated series. Mixed-archive selections from PReSto land here as a
heterogenous bag of records; the reducer (PCR/LASSO/SPLS/SIR) handles the
multi-variate compression.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _coerce_year(values) -> np.ndarray:
    """Convert age/year mixed types to a 1D float year-AD array.

    LiPD records sometimes use 'age' (years BP, 1950 reference) and
    sometimes 'year' (year AD). pylipd normalizes 'year' for us when
    possible; if all values look like BP, convert.
    """
    arr = np.asarray(values, dtype=float)
    if arr.size and np.nanmedian(arr) > 1e4:
        # Out of plausible year-AD range; assume BP and convert.
        arr = 1950.0 - arr
    return arr


def _aggregate_to_annual(years: np.ndarray, vals: np.ndarray,
                         year_axis: np.ndarray) -> np.ndarray:
    """Bin a record onto the common annual axis.

    LiPD records can be irregular / sub-annual / multi-decadal. We bin by
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


def _extract_records(pkl_path: Path) -> list[dict]:
    """Pull a flat list of records out of the legacy LiPD pickle.

    Tries pylipd first; falls back to walking the raw dict if pylipd
    isn't available or the pickle is already in legacy-dict form.
    """
    with pkl_path.open("rb") as f:
        raw = pickle.load(f)

    try:
        from pylipd.lipd import LiPD  # type: ignore
        L = LiPD()
        L.load_from_dict(raw if isinstance(raw, dict) else {})
        # essentials = dataframe of (dataSetName, paleoData_variableName,
        # paleoData_values, year, geo_meanLat, geo_meanLon,
        # geo_meanElev, archiveType, paleoData_proxy)
        df = L.get_timeseries_essentials()
        records = []
        for _, row in df.iterrows():
            records.append({
                "id":     str(row.get("dataSetName", "")),
                "var":    str(row.get("paleoData_variableName", "")),
                "lat":    float(row.get("geo_meanLat", np.nan)),
                "lon":    float(row.get("geo_meanLon", np.nan)),
                "elev":   float(row.get("geo_meanElev", np.nan)),
                "ptype":  str(row.get("archiveType", "unknown")).lower()
                          + "." + str(row.get("paleoData_proxy", "")),
                "years":  np.asarray(row.get("year", []), dtype=float),
                "values": np.asarray(row.get("paleoData_values", []), dtype=float),
            })
        return records
    except Exception as e:
        print(f"[lipd_to_baygmst] pylipd path failed ({e}); falling back to raw walk", file=sys.stderr)

    # Fallback: walk the raw dict. Legacy pickle shape:
    #   { datasetName: { 'geo': {...}, 'paleoData': [ { 'measurementTable': [ {'columns': [...]} ] } ], ... } }
    records = []
    if isinstance(raw, dict):
        for name, ds in raw.items():
            geo = ((ds or {}).get("geo") or {}).get("properties") or {}
            lat = float(geo.get("latitude", np.nan))
            lon = float(geo.get("longitude", np.nan))
            elev = float(geo.get("elevation", np.nan))
            archive = str(ds.get("archiveType", "unknown")).lower()
            for pd_entry in (ds.get("paleoData") or []):
                for mt in (pd_entry.get("measurementTable") or []):
                    cols = mt.get("columns") or []
                    year_col = next((c for c in cols
                                     if str(c.get("variableName", "")).lower() in ("year", "age")), None)
                    if year_col is None:
                        continue
                    years = _coerce_year(year_col.get("values", []))
                    for c in cols:
                        if c is year_col:
                            continue
                        vname = str(c.get("variableName", ""))
                        if not vname:
                            continue
                        records.append({
                            "id":     f"{name}__{vname}",
                            "var":    vname,
                            "lat":    lat,
                            "lon":    lon,
                            "elev":   elev,
                            "ptype":  f"{archive}.{c.get('proxy', vname)}",
                            "years":  years,
                            "values": np.asarray(c.get("values", []), dtype=float),
                        })
    return records


def build_csvs(pkl_path: Path, out_matrix: Path, out_metadata: Path,
               year_start: int = 1, year_end: int = 2000) -> None:
    records = _extract_records(pkl_path)
    if not records:
        raise SystemExit("No records extracted from LiPD pickle — cannot proceed.")

    year_axis = np.arange(year_start, year_end + 1, dtype=int)
    matrix = {"year": year_axis}
    meta_rows: list[tuple[str, float, float, float, str]] = []

    seen_ids: set[str] = set()
    for rec in records:
        rid = rec["id"]
        # Disambiguate duplicate IDs (different variables on same dataset).
        base = rid
        i = 1
        while rid in seen_ids:
            i += 1
            rid = f"{base}__{i}"
        seen_ids.add(rid)

        years = _coerce_year(rec["years"])
        vals  = np.asarray(rec["values"], dtype=float)
        if years.size == 0 or vals.size == 0 or years.size != vals.size:
            continue

        series = _aggregate_to_annual(years, vals, year_axis)
        if not np.isfinite(series).any():
            continue

        matrix[rid] = series
        meta_rows.append((rid, rec["lat"], rec["lon"], rec["elev"], rec["ptype"]))

    if len(matrix) <= 1:
        raise SystemExit("All LiPD records dropped during alignment — cannot proceed.")

    df_matrix = pd.DataFrame(matrix)
    df_meta = pd.DataFrame(meta_rows, columns=["", "lat", "lon", "elev", "ptype"])

    out_matrix.parent.mkdir(parents=True, exist_ok=True)
    df_matrix.to_csv(out_matrix, index=False)
    df_meta.to_csv(out_metadata, index=False)
    print(f"[lipd_to_baygmst] wrote {out_matrix} ({df_matrix.shape[0]} years x {df_matrix.shape[1]-1} proxies)")
    print(f"[lipd_to_baygmst] wrote {out_metadata} ({len(meta_rows)} records)")


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
