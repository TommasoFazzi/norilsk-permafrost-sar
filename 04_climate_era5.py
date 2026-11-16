"""Step 04: ERA5-Land air + soil temperature and snow at each point of interest.

Uses the Open-Meteo archive API (ERA5-Land repackaged, no key). For the
final submission, cite ERA5-Land (Copernicus CDS); values are the same model.
Output per point: daily CSV + annual summary (MAAT, thawing/freezing degree
days, deep soil temperature). Thawing degree days drive the seasonal
subsidence model used later with InSAR (subsidence ~ sqrt(cumulative TDD)).

    python 04_climate_era5.py
"""
import time

import pandas as pd
import requests

from config import DIRS, POINTS

URL = "https://archive-api.open-meteo.com/v1/archive"
START, END = "2015-01-01", pd.Timestamp.today().strftime("%Y-%m-%d")
# snowfall_sum is not provided for models=era5_land (all null): use snow_depth.
DAILY = ["temperature_2m_mean"]
# Soil temperature is kept for reference only: ERA5-Land gives a positive
# mean at 1-2.5 m in continuous permafrost (known bias). NOT USABLE as
# evidence of ground warming; use ESA CCI Permafrost (step 07).
HOURLY = ["soil_temperature_100_to_255cm", "snow_depth"]
# ERA5-Land cells are ~9 km: they give the regional forcing, not
# site-level differences. Three points (lowland, tundra, plateau) are enough.
CLIMATE_POINTS = ["norilsk_centre", "control_NW", "control_SE",
                  "plateau_reference"]
OUT = DIRS["climate"]


class RateLimited(Exception):
    pass


def fetch(lat, lon, tries=3):
    """Long hourly requests count as many calls on the free Open-Meteo API.
    Retry on 429 for the per-minute limit; give up on the hourly one."""
    for k in range(tries):
        r = requests.get(URL, timeout=300, params={
            "latitude": lat, "longitude": lon, "start_date": START,
            "end_date": END, "daily": ",".join(DAILY),
            "hourly": ",".join(HOURLY), "models": "era5_land",
            "timezone": "Asia/Krasnoyarsk"})
        if r.status_code != 429:
            break
        wait = 70 * (k + 1)
        print(f"  rate limit (429), waiting {wait} s...")
        time.sleep(wait)
    else:
        raise RateLimited(r.text[:200])
    r.raise_for_status()
    j = r.json()
    d = pd.DataFrame(j["daily"]).rename(columns={"time": "date"})
    d["date"] = pd.to_datetime(d["date"])
    h = pd.DataFrame(j["hourly"])
    h["date"] = pd.to_datetime(h["time"]).dt.floor("D")
    h = h.drop(columns="time").groupby("date").mean().reset_index()
    df = d.merge(h, on="date", how="left").set_index("date")
    df.attrs["grid"] = (j["latitude"], j["longitude"], j.get("elevation"))
    return df


def thaw_onset(t):
    """First day after 1 April starting 5 consecutive days with Tmean > 0."""
    t = t[t.index.month >= 4]
    warm = (t > 0).astype(int).rolling(5).sum() == 5
    return warm.idxmax() - pd.Timedelta(days=4) if warm.any() else pd.NaT


def annual(df):
    t = df["temperature_2m_mean"]
    g = df.groupby(df.index.year)
    out = pd.DataFrame({
        "maat_c": g["temperature_2m_mean"].mean(),
        "tdd": g["temperature_2m_mean"].apply(lambda s: s.clip(lower=0).sum()),
        "fdd": g["temperature_2m_mean"].apply(lambda s: -s.clip(upper=0).sum()),
        "soil_100_255_mean_c": g["soil_temperature_100_to_255cm"].mean(),
        "soil_100_255_max_c": g["soil_temperature_100_to_255cm"].max(),
        "snow_depth_max_m": g["snow_depth"].max(),
        "days": g.size(),
    })
    out["complete_year"] = out["days"] >= 365
    out["thaw_onset"] = [thaw_onset(t[t.index.year == y]) for y in out.index]
    num = out.columns.drop(["thaw_onset", "complete_year"])
    out[num] = out[num].round(2)
    return out


def main():
    for name in CLIMATE_POINTS:
        path = OUT / f"era5land_daily_{name}.csv"
        if path.exists():
            print(f"{name}: cached")
            continue
        lat, lon = POINTS[name]
        try:
            df = fetch(lat, lon)
        except RateLimited:
            print("Open-Meteo hourly limit reached. Re-run in ~1 hour: "
                  "points already saved are skipped.")
            break
        # cumulative thawing degree days within each year (for InSAR model)
        df["tdd_cum"] = (df["temperature_2m_mean"].clip(lower=0)
                         .groupby(df.index.year).cumsum())
        df.to_csv(path)
        print(f"{name}: grid cell {df.attrs['grid']}, {len(df)} days")
        time.sleep(5)

    rows = []
    for name in CLIMATE_POINTS:
        path = OUT / f"era5land_daily_{name}.csv"
        if path.exists():
            df = pd.read_csv(path, index_col="date", parse_dates=["date"])
            a = annual(df)
            a.insert(0, "point", name)
            rows.append(a)
    if rows:
        summary = pd.concat(rows)
        summary.to_csv(OUT / "era5land_annual.csv")
        print(summary.drop(columns="days").to_string())


if __name__ == "__main__":
    main()
