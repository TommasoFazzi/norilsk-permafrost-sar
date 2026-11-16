"""Step 17: weather-model troposphere between lowland and plateau, the
absolute (lowland vs bedrock) test, and soil moisture vs thaw.

The radar path to a lowland pixel (~80-150 m) crosses the air layer between
the lowland and the plateau reference (~660 m); the path to the reference
does not. Extra zenith delay of that layer, from ERA5 at the Sentinel-1 pass
time (00:40 UTC, hourly values interpolated), at the zone's climate point and
at the plateau point (Open-Meteo era5_seamless, no key):
  hydrostatic  ZHD = 2.2768 mm/hPa * (P_low - P_plateau) * dh_zone/dh_points
  wet          ZWD = 1e-6 * dh_zone * mean(k2' e/T + k3 e/T^2) * 1000  mm
(e from dew point, Magnus; endpoints averaged). More delay looks like motion
away from the satellite, so the apparent vertical motion of the zone relative
to the reference is
  d_app(t) = -(D(t) - D(t0)) / cos^2(inc)          (LOS / cos(inc) as step 09)

Checks
  1. Per pair: empirical height coefficient of step 09 (tropo_k_mm_per_km,
     fitted on 250-1100 m) vs the ERA5 layer gradient
     k_ERA5 = (G_sec - G_ref) / cos(inc), G = D / dh_points. Agreement means
     both see the same stratified troposphere.
  2. Absolute test (hypothesis A, tundra vs bedrock): zone-mean vertical
     series with NO troposphere model (conncomp reference) minus d_app, and
     the height+ramp series, fitted against sqrt(TDD) per zone and summer;
     two-sided p from Student t (n - 2 dof). Also how much d_app explains of
     the lowland common-mode series.
  3. Soil moisture: per-date wetland-minus-moss local series (step 16) vs
     ERA5 0-7 cm soil moisture and sqrt(TDD), pooled over summers with one
     intercept per summer.

Outputs: data/climate/era5_hourly_<point>.csv (cache), outputs/era5_tropo.json,
outputs/fig_era5_tropo.png

    python 17_era5_tropo.py
"""
import importlib
import json
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import requests
from scipy import stats

from common import ZONE_NAMES, save_json
from config import DIRS, OUTPUTS, POINTS, ZONE_CLIMATE

ts = importlib.import_module("09_timeseries")
loc = importlib.import_module("14_local_signal")
an = importlib.import_module("11_analysis")

URL = "https://archive-api.open-meteo.com/v1/archive"
VARS = ["temperature_2m", "dew_point_2m", "surface_pressure", "soil_moisture_0_to_7cm",
        "precipitation"]
PASS_UTC_H = 40 / 60          # Sentinel-1 path 136 descending, ~00:40 UTC
YEARS_FETCH = range(2017, 2027)
YEARS = (2020, 2025)
K2P, K3 = 22.97, 3.75463e5    # K/hPa, K^2/hPa
ZHD_MM_PER_HPA = 2.2768
LOWLAND_ZONES = [2, 4, 5, 6]


def fetch(point):
    """Hourly ERA5 (Jun-Oct of each summer) at one point, cached as CSV."""
    path = DIRS["climate"] / f"era5_hourly_{point}.csv"
    if path.exists():
        return pd.read_csv(path, index_col="time", parse_dates=["time"])
    lat, lon = POINTS[point]
    parts, elev = [], None
    for y in YEARS_FETCH:
        for k in range(3):
            r = requests.get(URL, timeout=300, params={
                "latitude": lat, "longitude": lon, "start_date": f"{y}-06-01",
                "end_date": f"{y}-10-15", "hourly": ",".join(VARS),
                "models": "era5_seamless", "timezone": "GMT"})
            if r.status_code != 429:
                break
            time.sleep(70 * (k + 1))
        r.raise_for_status()
        j = r.json()
        elev = j.get("elevation")
        parts.append(pd.DataFrame(j["hourly"]))
        time.sleep(2)
    df = pd.concat(parts)
    df["time"] = pd.to_datetime(df["time"])
    df["elevation_m"] = elev
    df.set_index("time").to_csv(path)
    return df.set_index("time")


def at_pass(df, dates):
    """Values at the pass time: linear between 00:00 and 01:00 UTC."""
    out = []
    for d in dates:
        t0 = pd.Timestamp(d)
        a, b = df.loc[t0], df.loc[t0 + pd.Timedelta(hours=1)]
        out.append(a * (1 - PASS_UTC_H) + b * PASS_UTC_H)
    return pd.DataFrame(out, index=dates)


def vapour_hpa(td_c):
    return 6.112 * np.exp(17.62 * td_c / (243.12 + td_c))


def wet_refractivity(t_c, td_c):
    T = t_c + 273.15
    e = vapour_hpa(td_c)
    return K2P * e / T + K3 * e / T ** 2


def layer_delay(low, ref, dh_zone):
    """Zenith delay (mm) of the air between a lowland zone and the plateau."""
    dh_pts = ref["elevation_m"] - low["elevation_m"]
    zhd = ZHD_MM_PER_HPA * (low["surface_pressure"] - ref["surface_pressure"]) * dh_zone / dh_pts
    nw = (wet_refractivity(low["temperature_2m"], low["dew_point_2m"])
          + wet_refractivity(ref["temperature_2m"], ref["dew_point_2m"])) / 2
    zwd = 1e-6 * dh_zone * nw * 1000
    return (zhd + zwd).astype(float), float(dh_pts.iloc[0])


def tfit(x, y):
    f = an.fit(np.asarray(x, float), np.asarray(y, float))
    dof = f["n"] - 2
    p = float(2 * stats.t.sf(abs(f["b"] / f["se_b"]), dof)) if dof > 0 and f["se_b"] > 0 else np.nan
    return {"b": round(float(f["b"]), 3), "se": round(float(f["se_b"]), 3), "n": int(f["n"]),
            "t": round(float(f["b"] / f["se_b"]), 2) if f["se_b"] > 0 else None,
            "p": None if np.isnan(p) else round(p, 4)}


def main():
    era = {p: fetch(p) for p in sorted(set(ZONE_CLIMATE.values()))}
    ctx = ts.Context()
    dem = ctx.dem.reshape(ctx.shape)
    ref_h = float(np.nanmean(dem[ctx.ref_patches > 0]))
    rep = {"tropo_k": [], "absolute": {}, "common_mode": {}}
    figdata = {}

    for year in YEARS:
        meta = json.loads((DIRS["ts"] / f"s1_{year}_meta.json").read_text())
        inc = np.radians(meta["incidence_median_deg"])
        dates_all = sorted({d for p in meta["pairs"] for d in p})
        ref_pt = at_pass(era["plateau_reference"], dates_all)

        # 1. empirical k (step 09) vs ERA5 layer gradient, per pair
        low_c = at_pass(era["norilsk_centre"], dates_all)
        D, dh_pts = layer_delay(low_c, ref_pt, ref_pt["elevation_m"].iloc[0] - low_c["elevation_m"].iloc[0])
        G = D / dh_pts                                        # mm per m
        if meta.get("tropo_k_mm_per_km"):
            for (r, s), k in zip(meta["pairs"], meta["tropo_k_mm_per_km"]):
                rep["tropo_k"].append({"pair": [r, s], "k_empirical": k,
                                       "k_era5": round(float((G[s] - G[r]) / np.cos(inc) * 1000), 2)})

        # 2. absolute test: NO-troposphere series (reference component) - ERA5
        items = sorted(ts.find_pairs()[year])
        b_raw, coh, inc_map = ts.load_year(items, ctx, "reference")
        b, c, pairs, dates, _, _ = ts.preprocess(b_raw, coh, items, ctx, "none", False, "reference")
        del b_raw
        _, vert_none, *_ = ts.solve(dates, pairs, b, c, inc_map, ctx)
        with rasterio.open(DIRS["ts"] / f"s1_{year}_vert.tif") as s:
            vert_hr = s.read().astype("float32")
            hr_dates = list(s.descriptions)
        assert hr_dates == dates, "run 09 --conncomp reference first"
        x_by_pt = {p: loc.sqrt_tdd(dates) if p == "norilsk_centre" else
                   np.sqrt(an.tdd_series(p).reindex(pd.to_datetime(dates)).values) for p in era}
        ref_d = ref_pt.loc[dates]
        res, cm_none, cm_app = {}, [], []
        for zc in LOWLAND_ZONES:
            name = ZONE_NAMES[zc]
            m = ctx.z == zc
            pt = ZONE_CLIMATE[name]
            dh_zone = ref_h - float(np.nanmean(dem[m]))
            Dz, _ = layer_delay(at_pass(era[pt], dates), ref_d, dh_zone)
            d_app = -(Dz - Dz.iloc[0]).values / np.cos(inc) ** 2
            y_none = np.array([np.nanmean(v[m]) for v in vert_none])
            y_hr = np.array([np.nanmean(v[m]) for v in vert_hr])
            x = x_by_pt[pt]
            res[name] = {"none": tfit(x, y_none), "era5": tfit(x, y_none - d_app),
                         "height_ramp": tfit(x, y_hr),
                         "end_mm": {"none": round(float(y_none[-1]), 1),
                                    "era5": round(float(y_none[-1] - d_app[-1]), 1),
                                    "height_ramp": round(float(y_hr[-1]), 1)},
                         "rms_about_fit_mm": {
                             k: round(float(np.std(yy - np.polyval(np.polyfit(x, yy, 1), x))), 1)
                             for k, yy in (("none", y_none), ("era5", y_none - d_app), ("height_ramp", y_hr))}}
            if zc in (5, 6):
                cm_none.append(y_none)
                cm_app.append(d_app)
        cm_none, cm_app = np.mean(cm_none, 0), np.mean(cm_app, 0)
        rep["common_mode"][str(year)] = {
            "dates": dates, "tundra_none_mm": np.round(cm_none, 1).tolist(),
            "era5_apparent_mm": np.round(cm_app, 1).tolist(),
            "r": round(float(np.corrcoef(cm_none, cm_app)[0, 1]), 3),
            "regression_slope": round(float(np.polyfit(cm_app, cm_none, 1)[0]), 2)}
        rep["absolute"][str(year)] = res
        figdata[year] = (dates, cm_none, cm_app)

    kk = pd.DataFrame(rep["tropo_k"])
    if len(kk) > 2:
        rep["tropo_k_summary"] = {"n_pairs": len(kk),
                                  "r": round(float(np.corrcoef(kk.k_empirical, kk.k_era5)[0, 1]), 3),
                                  "slope_emp_on_era5": round(float(np.polyfit(kk.k_era5, kk.k_empirical, 1)[0]), 2),
                                  "rms_emp": round(float(kk.k_empirical.std()), 2),
                                  "rms_era5": round(float(kk.k_era5.std()), 2)}

    # 3. soil moisture vs thaw for the wetland local series (step 16)
    rob = OUTPUTS / "robustness.json"
    if rob.exists():
        st = json.loads(rob.read_text())["stefan_v1"]
        rows = []
        for year in YEARS:
            s = st[str(year)]
            sm = np.mean([at_pass(era[p], s["dates"])["soil_moisture_0_to_7cm"].values
                          for p in ("norilsk_centre", "control_NW", "control_SE")], axis=0)
            pr = np.mean([[era[p].loc[pd.Timestamp(d) - pd.Timedelta(days=3):pd.Timestamp(d), "precipitation"].sum()
                           for d in s["dates"]] for p in ("norilsk_centre", "control_NW", "control_SE")], axis=0)
            for d, y, x, m, p3 in zip(s["dates"], s["wet_minus_moss_mm"], s["sqrt_tdd"], sm, pr):
                rows.append({"summer": year, "date": d, "y": y, "sqrt_tdd": x, "sm": m, "precip_3d_mm": p3})
        df = pd.DataFrame(rows)
        X_year = pd.get_dummies(df.summer, prefix="s").astype(float).values
        out = {"r_sqrt_tdd_vs_sm": round(float(df.sqrt_tdd.corr(df.sm)), 3), "models": {}}
        n = len(df)
        for name, cols in (("thaw", ["sqrt_tdd"]), ("moisture", ["sm"]), ("both", ["sqrt_tdd", "sm"])):
            X = np.c_[X_year, df[cols].values]
            coef, *_ = np.linalg.lstsq(X, df.y.values, rcond=None)
            resid = df.y.values - X @ coef
            dof = n - X.shape[1]
            s2 = resid @ resid / dof
            se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
            k = X.shape[1]
            aic = n * np.log(resid @ resid / n) + 2 * k
            out["models"][name] = {c: {"coef": round(float(coef[-len(cols) + i]), 3),
                                       "t": round(float(coef[-len(cols) + i] / se[-len(cols) + i]), 2)}
                                   for i, c in enumerate(cols)}
            out["models"][name].update({"rss": round(float(resid @ resid), 2), "aic": round(float(aic), 2)})
        out["data"] = df.round(3).to_dict("records")
        rep["moisture_vs_thaw"] = out
    save_json(OUTPUTS / "era5_tropo.json", rep)

    print(json.dumps({k: v for k, v in rep.items() if k not in ("tropo_k",)}, indent=1, default=str)[:6000])

    fig, axs = plt.subplots(1, len(figdata), figsize=(6 * len(figdata), 4), squeeze=False)
    for ax, (year, (dates, cmn, cma)) in zip(axs[0], figdata.items()):
        t = pd.to_datetime(dates)
        ax.plot(t, cmn, "o-", label="tundra, no tropo model")
        ax.plot(t, cma, "s--", label="ERA5 apparent (layer delay)")
        ax.plot(t, cmn - cma, "^-", label="tundra - ERA5")
        ax.set_title(f"{year}: tundra controls vs plateau [mm, + up]", fontsize=9)
        ax.legend(fontsize=7)
        ax.tick_params(axis="x", labelrotation=30, labelsize=7)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_era5_tropo.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
