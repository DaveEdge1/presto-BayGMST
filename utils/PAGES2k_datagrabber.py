import cfr
import xarray as xr
import pandas as pd
import numpy as np
import sys

pdb = cfr.ProxyDatabase().fetch('PAGES2kv2')
df = pdb.to_df()


# Example: df has columns
# ['pid', 'lat', 'lon', 'elev', 'ptype', 'time', 'value']

# Step 1: explode list-columns so each proxy-year is one row
long_df = df[['pid', 'lat', 'lon', 'elev', 'ptype', 'time', 'value']].explode(
    ['time', 'value'],
    ignore_index=True
)

# Step 2: make sure types are numeric
long_df['time'] = pd.to_numeric(long_df['time'], errors='coerce')
long_df['value'] = pd.to_numeric(long_df['value'], errors='coerce')
long_df['year'] = np.floor(long_df['time']).astype('Int64')

# Step 3: create full year index, e.g. 1 to 2000
full_years = pd.Index(range(1, 2001), name='year')

# Step 4: pivot to wide matrix (and account for multiple observations per year by taking the mean)
proxy_matrix = (
    long_df
    .assign(year=np.floor(long_df['time']).astype('Int64'))
    .pivot_table(index='year', columns='pid', values='value', aggfunc='mean')
    .reindex(range(1, 2001))
)

# Step 5: metadata table
proxy_meta = (
    df[['pid', 'lat', 'lon', 'elev', 'ptype']]
    .drop_duplicates('pid')
    .set_index('pid')
    .loc[proxy_matrix.columns]
)

# Step 6: combine metadata and data into one DataFrame (optional)
meta_rows = pd.DataFrame(
    [proxy_meta.loc[proxy_matrix.columns, 'lat'],
     proxy_meta.loc[proxy_matrix.columns, 'lon'],
     proxy_meta.loc[proxy_matrix.columns, 'elev'],
     proxy_meta.loc[proxy_matrix.columns, 'ptype']],
    index=['lat', 'lon', 'elev', 'ptype']
)

combined = pd.concat([meta_rows, proxy_matrix])
# print(combined)
# proxy_matrix.to_csv('data/PAGES2K_proxy_matrix_1-2000.csv', index=True)
# proxy_meta.to_csv('data/PAGES2K_proxy_metadata_1-2000.csv', index=True)
# combined.to_csv('data/PAGES2K_proxy_combined_data_1-2000.csv', index=True)




##############
### --- SCREENING: Screen out proxies that are poor local temperature thermometers (p>0.05) during 1900-2000
import xarray as xr
from scipy.stats import pearsonr
path = "data/HadCRUT.5.1.0.0.analysis.anomalies.ensemble_mean.nc"
ds_temp = xr.open_dataset(path)
print(ds_temp)

# HadCRUT variable
tas = ds_temp["tas_mean"]

# Convert monthly HadCRUT to annual means
tas_ann = tas.groupby("time.year").mean("time", skipna=True)

# HadCRUT longitudes are -180 to 180; many PAGES2k lons are 0 to 360
def lon_to_hadcrut(lon):
    lon = float(lon)
    return ((lon + 180) % 360) - 180

def screen_one_proxy(pid, calib_start=1900, calib_end=2000, min_overlap=10):
    lat = float(proxy_meta.loc[pid, "lat"])
    lon = lon_to_hadcrut(proxy_meta.loc[pid, "lon"])

    proxy_ts = proxy_matrix[pid].loc[calib_start:calib_end]

    local_temp = (
        tas_ann
        .sel(latitude=lat, longitude=lon, method="nearest")
        .to_series()
        .loc[calib_start:calib_end]
    )

    aligned = pd.concat(
        [proxy_ts.rename("proxy"), local_temp.rename("temp")],
        axis=1
    ).dropna()

    n = len(aligned)

    if n < min_overlap:
        return pd.Series({
            "lat": lat,
            "lon": lon,
            "n_overlap": n,
            "r": np.nan,
            "p": np.nan,
            "keep": False,
            "nearest_lat": np.nan,
            "nearest_lon": np.nan,
        })

    r, p = pearsonr(aligned["proxy"], aligned["temp"])

    nearest = tas_ann.sel(latitude=lat, longitude=lon, method="nearest")

    return pd.Series({
        "lat": lat,
        "lon": lon,
        "n_overlap": n,
        "r": r,
        "p": p,
        "keep": p < 0.05,
        "nearest_lat": float(nearest["latitude"]),
        "nearest_lon": float(nearest["longitude"]),
    })

screening = pd.DataFrame({
    pid: screen_one_proxy(pid)
    for pid in proxy_matrix.columns
}).T

keep_pids = screening.index[screening["keep"]]

proxy_matrix_screened = proxy_matrix[keep_pids]
proxy_meta_screened = proxy_meta.loc[keep_pids]

proxy_matrix_screened.to_csv("data/PAGES2K_proxy_matrix_screened_1900-2000.csv")
proxy_meta_screened.to_csv("data/PAGES2K_proxy_metadata_screened_1900-2000.csv")
