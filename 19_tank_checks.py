"""Step 19: is the tank-5 winter signal of step 18 real? Checks and context.

Runs on a wider crop of the same HyP3 jobs (step 02 download --name
norilsk_tank_nowm --out norilsk_tank_wide --crop-km 8 --crop-point
burst_tanks_centre), which holds tank 5 and the OSM storage tanks of the
Norilsk depots 9-13 km east. The step-18 machinery is reused unchanged
(3x3 window minus winter-coherent targets 0.3-1.5 km away in the same SNAPHU
component, segments between the snowmelt / snow-onset gaps, vertical =
LOS / cos(incidence)).

C1 Unwrapping. Winter 2019-20 pairs: triplet closure V12 + V23 - V13 (LOS mm)
   of the tank window, of each of its 9 pixels and of the null targets; an
   unwrapping error is a whole cycle (27.7 mm LOS). Pixels more than half a
   cycle from the window median in a pair are flagged.
C2 Robust window: median of the window pixels minus median of the ring.
C3 Weighting: SBAS weighted by the window coherence g (w = g^2 / (1 - g^2));
   leave-one-pair-out change Dec 6 -> Mar 11 and -> Apr 4.
C4 DEM error: the winter 2019-20 series fitted with a linear trend plus a
   height-error term B_perp / (R sin(inc)) per date (pair baselines from the
   HyP3 .txt files, made per-date by the same SBAS), then removed.
C5 Tank-only comparison: the same interval changes for the OSM storage tanks
   in the crop and the bright structures of the tank-5 row (amplitude peaks
   >= ROW_PCT percentile within ROW_M), next to the step-18 null (bright
   targets 1.5-5.5 km away with tank 5's winter coherence).
C6 Per-date table, Dec -> Apr of each winter, all of the above.
C7 Snow and wind at tank 5 (ERA5 hourly via Open-Meteo, cached): snowfall,
   snow depth and blowing-snow hours (wind >= BLOW_WIND m/s, air <= BLOW_TMAX
   degC, snow on the ground) between consecutive acquisitions. Dry snow
   delays the radar by ~1 mm LOS per mm of snow water equivalent, so a drift
   at the tank would look like settlement.

Outputs: outputs/tank_checks.json, outputs/tank_winters.csv,
outputs/tank_points_series.csv (segment series of tank 5 and the row
structures; each segment starts at 0), outputs/fig_tank_checks.png,
data/climate/era5_hourly_winter_tank5.csv

    python 19_tank_checks.py
"""
import importlib
import json
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.ndimage import maximum_filter
from shapely.geometry import shape

from common import save_json
from config import DIRS, OUTPUTS, POINTS, S1_WAVELENGTH

ts = importlib.import_module("09_timeseries")
bias = importlib.import_module("15_bias_checks")
t18 = importlib.import_module("18_tank_series")

PROD = DIRS["hyp3"] / "norilsk_tank_wide"
WINTER_2020 = ("20191101", "20200510")
D0, D_MAR, D_APR = "20191206", "20200311", "20200404"
WINTER_STARTS = {"2017-18": "20171204", "2018-19": "20181211", "2019-20": "20191206"}
CYCLE_MM = S1_WAVELENGTH / 2 * 1000        # one phase cycle in LOS mm
ROW_M, ROW_PCT, ROW_MIN_M = 500, 99, 60
ERA5_URL = "https://archive-api.open-meteo.com/v1/archive"
ERA5_VARS = ["temperature_2m", "snowfall", "snow_depth", "wind_speed_10m",
             "wind_direction_10m", "wind_gusts_10m"]
BLOW_WIND, BLOW_TMAX, BLOW_SNOW_M = 7.0, -2.0, 0.05


# ----------------------------------------------------------------- targets
def row_structures(maps, g, tank, px):
    """Bright amplitude peaks of the tank-5 row (probable tanks; OSM has none)."""
    a = 10 * np.log10(maps["amp"])
    d = t18.distance_m(*tank, g, px)
    peak = ((a == maximum_filter(a, size=3)) & (d <= ROW_M) & (d > ROW_MIN_M)
            & (a >= np.nanpercentile(a, ROW_PCT)))
    out = {}
    for r, c in sorted(zip(*np.nonzero(peak)), key=lambda p: d[p]):
        dy, dx = tank[0] - r, c - tank[1]
        bearing = "N NE E SE S SW W NW".split()[int(((np.degrees(np.arctan2(dx, dy)) + 22.5) % 360) // 45)]
        out[f"row_{d[r, c]:.0f}m_{bearing}"] = (int(r), int(c))
    return out


def osm_tanks(g):
    gj = json.loads((DIRS["osm"] / "osm_tanks_plants.geojson").read_text())
    out = {}
    for f in gj["features"]:
        if f["properties"].get("man_made") != "storage_tank":
            continue
        c = shape(f["geometry"]).centroid
        r, cc = t18.rc_of(c.y, c.x, g)
        if 2 <= r < g["height"] - 2 and 2 <= cc < g["width"] - 2:
            out[f"osm_{f['properties']['osm_id'].split('/')[1]}"] = (r, cc)
    return out


# ----------------------------------------------------------------- C1-C4
def window_pixels(pairs, g, point, px, stable):
    """Per pair: LOS mm of each 3x3 pixel in the window's SNAPHU component,
    and the mean and median of the ring's stable targets (NaN if unusable)."""
    win, ring = t18.window_ring(*point, g, px, stable)
    P = np.full((len(pairs), win.size), np.nan)
    ring_mean, ring_med, coh_win = (np.full(len(pairs), np.nan) for _ in range(3))
    for k, (r, s, stem) in enumerate(pairs):
        b, cc, coh = t18.load(stem, g)
        ok = (np.isfinite(b) & (coh >= t18.COH_MIN)).ravel()
        bf, ccf = b.ravel(), cc.ravel()
        coh_win[k] = np.nanmean(coh.ravel()[win])
        w = win[ok[win]]
        if w.size == 0:
            continue
        lab = np.bincount(ccf[w].astype(int)).argmax()
        rg = ring[ok[ring] & (ccf[ring] == lab)]
        if rg.size < t18.MIN_RING_PS:
            continue
        sel = ok[win] & (ccf[win] == lab)
        P[k, sel] = bf[win[sel]]
        ring_mean[k], ring_med[k] = bf[rg].mean(), np.median(bf[rg])
    return P, ring_mean, ring_med, coh_win


def closures(pairs, v):
    """Triplet closures V(a,b) + V(b,c) - V(a,c) (LOS mm)."""
    idx = {(r, s): k for k, (r, s, _) in enumerate(pairs)}
    out = []
    for (a, b), k1 in idx.items():
        for (b2, c), k2 in idx.items():
            if b2 == b and (a, c) in idx:
                x = v[k1] + v[k2] - v[idx[(a, c)]]
                if np.isfinite(x):
                    out.append(float(x))
    return np.array(out)


def segment_series(pairs, v, d0, d1, weights=None):
    """SBAS series (units of v) of the connected set of dates holding d0 and
    d1, optionally weighted: (dates, x, design, used pair indices)."""
    use = [k for k in range(len(pairs)) if np.isfinite(v[k])]
    seg = next((c for c in ts.date_components([pairs[k][:2] for k in use])
                if d0 in c and d1 in c), None)
    if seg is None:
        return None
    ks = [k for k in use if pairs[k][0] in seg]
    A = ts.design(seg, [pairs[k][:2] for k in ks])
    y = v[ks].astype(float)
    if weights is not None:
        sw = np.sqrt(weights[ks])
        A, y = A * sw[:, None], y * sw
    return seg, np.r_[0.0, np.linalg.lstsq(A, y, rcond=None)[0]], ks


def segment_change(pairs, v, d0, d1, weights=None):
    s = segment_series(pairs, v, d0, d1, weights)
    if s is None:
        return np.nan
    seg, x, _ = s
    return float(x[seg.index(d1)] - x[seg.index(d0)])


def height_corrected_change(pairs, v, d0, d1, height_col):
    """Pair baselines are differences of per-date baselines, so a height
    error cannot be told apart from free per-date positions; as in standard
    SBAS, fit the per-date series with a linear trend plus the per-date
    baseline term (-> dh in m), then remove that term."""
    s = segment_series(pairs, v, d0, d1)
    if s is None:
        return np.nan, np.nan
    seg, x, ks = s
    A = ts.design(seg, [pairs[k][:2] for k in ks])
    b = np.r_[0.0, np.linalg.lstsq(A, height_col[ks], rcond=None)[0]]   # mm per m, per date
    t = (pd.to_datetime(seg) - pd.Timestamp(seg[0])).days.values / 365.25
    coef = np.linalg.lstsq(np.c_[np.ones_like(t), t, b], x, rcond=None)[0]
    xc = x - coef[2] * b
    return float(xc[seg.index(d1)] - xc[seg.index(d0)]), float(coef[2])


def unwrapping_checks(win_pairs, P, ring_mean, V_null):
    tank = np.nanmean(P, 1) - ring_mean
    cl_tank = closures(win_pairs, tank)
    cl_px = np.concatenate([closures(win_pairs, P[:, i] - ring_mean) for i in range(P.shape[1])])
    cl_null = np.concatenate([closures(win_pairs, V_null[:, j]) for j in range(V_null.shape[1])])
    dev = np.abs(P - np.nanmedian(P, 1, keepdims=True))
    half = CYCLE_MM / 2
    stat = lambda c: {"n": int(c.size), "rms_mm": round(float(np.sqrt(np.mean(c ** 2))), 2),
                      "frac_over_half_cycle": round(float(np.mean(np.abs(c) > half)), 3)}
    return {"cycle_mm_los": round(CYCLE_MM, 2), "tank_window": stat(cl_tank),
            "tank_pixels": stat(cl_px), "null_windows": stat(cl_null),
            "pixel_offsets_over_half_cycle": int(np.nansum(dev > half)),
            "pixel_values": int(np.isfinite(P).sum()),
            "pairs_with_offset_pixel": [f"{r}_{s}" for (r, s, _), row in zip(win_pairs, dev)
                                        if np.nansum(row > half) > 0]}, cl_tank, cl_null


# ----------------------------------------------------------------- C7
def era5_winter():
    path = DIRS["climate"] / "era5_hourly_winter_tank5.csv"
    if path.exists():
        return pd.read_csv(path, index_col="time", parse_dates=["time"])
    lat, lon = POINTS[t18.TANK_PT]
    parts = []
    for y in (2017, 2018, 2019):
        for k in range(3):
            r = requests.get(ERA5_URL, timeout=300, params={
                "latitude": lat, "longitude": lon, "start_date": f"{y}-11-01",
                "end_date": f"{y + 1}-05-31", "hourly": ",".join(ERA5_VARS),
                "models": "era5_seamless", "wind_speed_unit": "ms", "timezone": "GMT"})
            if r.status_code != 429:
                break
            time.sleep(70 * (k + 1))
        r.raise_for_status()
        parts.append(pd.DataFrame(r.json()["hourly"]))
        time.sleep(2)
    df = pd.concat(parts)
    df["time"] = pd.to_datetime(df["time"])
    df.set_index("time").to_csv(path)
    return df.set_index("time")


def weather_between(era, d0, d1):
    w = era.loc[pd.Timestamp(d0) + pd.Timedelta("40min"):pd.Timestamp(d1) + pd.Timedelta("40min")]
    blow = ((w.wind_speed_10m >= BLOW_WIND) & (w.temperature_2m <= BLOW_TMAX)
            & (w.snow_depth >= BLOW_SNOW_M))
    dirs = np.radians(w.wind_direction_10m[blow])
    prevailing = (round(float(np.degrees(np.arctan2(np.sin(dirs).mean(), np.cos(dirs).mean())) % 360))
                  if blow.any() else None)
    return {"snowfall_cm": round(float(w.snowfall.sum()), 1),
            "snow_depth_change_cm": round(float((w.snow_depth.iloc[-1] - w.snow_depth.iloc[0]) * 100), 1),
            "blowing_snow_h": int(blow.sum()), "max_wind_ms": round(float(w.wind_speed_10m.max()), 1),
            "blowing_from_deg": prevailing}


# ----------------------------------------------------------------- main
def main():
    pairs = t18.product_pairs(PROD)
    if not pairs:
        raise SystemExit(f"No products in {PROD}: run step 02 download with --out norilsk_tank_wide.")
    g = t18.reference_grid(pairs[0][2])
    px = abs(g["transform"].a)
    maps = t18.site_maps(pairs, g)
    bright = maps["amp"] >= np.nanpercentile(maps["amp"], t18.BRIGHT_PCT)
    stable = maps["winter_coh"] >= t18.COH_STABLE
    tank = t18.rc_of(*POINTS[t18.TANK_PT], g)
    row, osm = row_structures(maps, g, tank, px), osm_tanks(g)
    nulls = t18.null_points(maps, bright, g, [tank], tank, px)
    named = {"tank5": tank, **row, **osm}
    last = t18.last_before(pairs)
    pre = [p for p in pairs if p[1] <= last]
    V, C, inc = t18.relative_values(pre, g, list(named.values()) + nulls, px, stable)
    cosi = np.cos(np.radians(inc))
    rep = {"crop": str(PROD.name), "incidence_deg": round(inc, 1), "tank_rc": list(tank),
           "tank_winter_coh": round(float(maps["winter_coh"][tank]), 2),
           "row_structures": {n: {"rc": list(p), "latlon": list(t18.latlon_of(*p, g)),
                                  "winter_coh": round(float(maps["winter_coh"][p]), 2),
                                  "amp_change_2019_2020_db": round(float(maps["2020"][p] - maps["2019"][p]), 1)}
                              for n, p in row.items()},
           "osm_tanks": {n: list(t18.latlon_of(*p, g)) for n, p in osm.items()},
           "n_null": len(nulls)}

    # C1-C4 on winter 2019-20
    wi = [k for k, p in enumerate(pre) if p[0] >= WINTER_2020[0] and p[1] <= WINTER_2020[1]]
    wp = [pre[k] for k in wi]
    P, ring_mean, ring_med, coh_win = window_pixels(wp, g, tank, px, stable)
    rep["c1_unwrapping"], cl_tank, cl_null = unwrapping_checks(wp, P, ring_mean, V[wi][:, len(named):])
    v_mean = np.nanmean(P, 1) - ring_mean
    v_med = np.nanmedian(P, 1) - ring_med
    gam = np.clip(coh_win, 0.05, 0.95)
    w = gam ** 2 / (1 - gam ** 2)
    bp = np.array([bias.pair_baseline(f"{stem}_unw_phase.tif") for *_, stem in wp])
    hcol = bp[:, 0] / (bp[:, 1] * np.sin(np.radians(inc))) * 1000        # LOS mm per m
    vert = lambda x: round(float(x) / cosi, 1)
    rep["c2_c4_change_vertical_mm"] = {}
    for d1 in (D_MAR, D_APR):
        base = segment_change(wp, v_mean, D0, d1)
        med = segment_change(wp, v_med, D0, d1)
        wtd = segment_change(wp, v_mean, D0, d1, weights=w)
        dem, dh = height_corrected_change(wp, v_mean, D0, d1, hcol)
        jack = []
        for k in range(len(wp)):
            if np.isfinite(v_mean[k]):
                vv = v_mean.copy()
                vv[k] = np.nan
                jack.append((f"{wp[k][0]}_{wp[k][1]}", segment_change(wp, vv, D0, d1)))
        jv = np.array([x for _, x in jack])
        jmax = max(jack, key=lambda z: abs(np.nan_to_num(z[1] - base, nan=1e9)))
        rep["c2_c4_change_vertical_mm"][f"{D0}_{d1}"] = {
            "window_mean": vert(base), "window_median": vert(med), "coherence_weighted": vert(wtd),
            "with_height_error_term": vert(dem), "height_error_m": round(dh, 1),
            "jackknife_min": vert(np.nanmin(jv)), "jackknife_max": vert(np.nanmax(jv)),
            "jackknife_n_disconnected": int(np.isnan(jv).sum()),
            "most_influential_pair": jmax[0],
            "without_it": None if np.isnan(jmax[1]) else vert(jmax[1])}

    # C5-C6 interval and per-date changes for all targets
    Vv = V / cosi
    segs = [t18.invert_segments(pre, Vv[:, j]) for j in range(Vv.shape[1])]
    k = len(named)
    intervals = [(d[0], d[-1]) for d, _, _ in segs[0] if len(d) >= t18.MIN_SEG_DATES]
    is_osm = np.array([n.startswith("osm_") for n in named])
    is_row = np.array([n.startswith("row_") for n in named])

    def compare(d0, d1):
        ch = np.array([t18.interval_change(s, d0, d1) for s in segs])
        tank_v, mine, nv = ch[0], ch[:k], ch[k:][np.isfinite(ch[k:])]
        tanks = mine[1:][np.isfinite(mine[1:])]
        out = {"tank5": None if np.isnan(tank_v) else round(float(tank_v), 1),
               "row": {n: (None if np.isnan(x) else round(float(x), 1))
                       for n, x, r in zip(named, mine, is_row) if r},
               "osm_tanks_mean": round(float(np.nanmean(mine[is_osm])), 1) if np.isfinite(mine[is_osm]).any() else None,
               "osm_tanks_sd": round(float(np.nanstd(mine[is_osm])), 1) if np.isfinite(mine[is_osm]).sum() > 1 else None,
               "osm_tanks_n": int(np.isfinite(mine[is_osm]).sum()),
               "tanks_and_row_as_low_or_lower": (int(np.sum(tanks <= tank_v)) if np.isfinite(tank_v) else None),
               "tanks_and_row_n": int(tanks.size),
               "null_mean": round(float(nv.mean()), 1) if nv.size else None,
               "null_sd": round(float(nv.std()), 1) if nv.size else None, "null_n": int(nv.size)}
        if nv.size >= 20 and np.isfinite(tank_v):
            out["z"] = round(float((tank_v - nv.mean()) / nv.std()), 2)
            out["p"] = round(float(np.mean(np.abs(nv - nv.mean()) >= abs(tank_v - nv.mean()))), 3)
        return out

    rep["c5_intervals"] = {f"{a}_{b}": compare(a, b) for a, b in intervals}
    pd.DataFrame([{"point": n, "latitude": t18.latlon_of(*named[n], g)[0],
                   "longitude": t18.latlon_of(*named[n], g)[1], "segment": si,
                   "segment_rms_mm": round(rms, 2), "date": d, "vertical_mm": round(float(x), 2)}
                  for n, s in zip(named, segs[:k]) if not n.startswith("osm_")
                  for si, (dd, xx, rms) in enumerate(s) for d, x in zip(dd, xx)]
                 ).to_csv(OUTPUTS / "tank_points_series.csv", index=False)
    era = era5_winter()
    rows = []
    dates = sorted({d for p in pre for d in p[:2]})
    for wname, d0 in WINTER_STARTS.items():
        end = f"{int(d0[:4]) + 1}0430"
        prev = d0
        for d in [x for x in dates if d0 < x <= end]:
            c = compare(d0, d)
            if c["tank5"] is None and c["null_n"] < 20:
                continue
            rows.append({"winter": wname, "start": d0, "date": d, "tank5_mm": c["tank5"],
                         **{f"{n}_mm": x for n, x in c["row"].items()},
                         "osm_tanks_mean_mm": c["osm_tanks_mean"], "osm_tanks_sd_mm": c["osm_tanks_sd"],
                         "osm_tanks_n": c["osm_tanks_n"], "null_mean_mm": c["null_mean"],
                         "null_sd_mm": c["null_sd"], "null_n": c["null_n"],
                         "z_tank5": c.get("z"), "p_tank5": c.get("p"),
                         **{f"since_prev_{key}": v for key, v in weather_between(era, prev, d).items()}})
            prev = d
    table = pd.DataFrame(rows)
    table.to_csv(OUTPUTS / "tank_winters.csv", index=False)

    # C7 season totals and link between tank steps and blowing snow
    rep["c7_weather"] = {}
    for wname, d0 in WINTER_STARTS.items():
        sub = table[table.winter == wname]
        if sub.empty:
            continue
        rep["c7_weather"][wname] = {"from": d0, "to": sub.date.iloc[-1],
                                    **weather_between(era, d0, sub.date.iloc[-1])}
        steps = sub.tank5_mm.diff().fillna(sub.tank5_mm)
        ok = steps.notna() & sub.since_prev_blowing_snow_h.notna()
        if ok.sum() >= 5:
            rep["c7_weather"][wname]["r_tank_step_vs_blowing_h"] = round(
                float(np.corrcoef(steps[ok], sub.since_prev_blowing_snow_h[ok])[0, 1]), 2)
    save_json(OUTPUTS / "tank_checks.json", rep)
    print(json.dumps(rep, indent=1, default=str))
    print(table.to_string(max_cols=40))
    figure(rep, table, cl_tank, cl_null)


def figure(rep, table, cl_tank, cl_null):
    fig, axs = plt.subplots(1, 3, figsize=(19, 5))
    bins = np.linspace(-1.2, 1.2, 49)
    axs[0].hist(cl_null / CYCLE_MM, bins=bins, density=True, color="grey", alpha=0.6, label="null targets")
    axs[0].hist(cl_tank / CYCLE_MM, bins=bins, density=True, histtype="step", color="red", lw=2, label="tank 5")
    axs[0].set_xlabel("triplet closure, winter 2019-20 [phase cycles]")
    axs[0].set_title("unwrapping check: errors would sit at +-1 cycle", fontsize=9)
    axs[0].legend(fontsize=8)
    for ax, wname in zip(axs[1:], ("2018-19", "2019-20")):
        s = table[table.winter == wname]
        if s.empty:
            ax.set_title(f"winter {wname}: no data", fontsize=9)
            continue
        t = pd.to_datetime(s.date)
        ax.fill_between(t, s.null_mean_mm - 2 * s.null_sd_mm, s.null_mean_mm + 2 * s.null_sd_mm,
                        color="grey", alpha=0.25, label="comparable targets, ±2 sd")
        if s.osm_tanks_mean_mm.notna().any():
            ax.plot(t, s.osm_tanks_mean_mm, "s--", color="C2", ms=4, label="OSM depot tanks, mean")
        for c in [c for c in s.columns if c.startswith("row_") and c.endswith("_mm")]:
            ax.plot(t, s[c], ".:", lw=0.8, label=c.removesuffix("_mm"))
        ax.plot(t, s.tank5_mm, "o-", color="red", lw=2, label="tank 5")
        ax2 = ax.twinx()
        ax2.bar(t, s.since_prev_blowing_snow_h, width=5, color="C0", alpha=0.25, align="edge")
        ax2.set_ylabel("blowing-snow hours since previous date", color="C0", fontsize=8)
        ax.axhline(0, color="k", lw=0.5)
        ax.tick_params(axis="x", labelrotation=30, labelsize=7)
        ax.set_ylabel(f"vertical change since {s.start.iloc[0]} [mm], + up")
        ax.set_title(f"winter {wname}: tank 5 vs other targets", fontsize=9)
        ax.legend(fontsize=6, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_tank_checks.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
