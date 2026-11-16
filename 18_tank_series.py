"""Step 18: TPP-3 tank no. 5 - was the ground under it sinking before the collapse?

Zhang et al. (2022, Remote Sens. 14, 5036) report 60-70 mm of subsidence at
the site between Aug 2018 and the collapse (2020-05-29), about -35 mm/yr,
relative to hills to the SE, and attribute it to permafrost degradation;
Rostekhnadzor (Nov 2020) found no permafrost degradation and blamed the pile
foundation and its operation. Damage comes from differential settlement, so
the test is local: did the tank move relative to the structures around it,
before the collapse, more than comparable places do?

A. Single-burst series (step 02 --name norilsk_tank_nowm --no-water-mask:
   burst 136_290630_IW2, 5x1 looks ~20 m, all seasons 2017-06 .. 2020-10,
   cropped to +-6 km). The HyP3 water mask zeroes the phase on the tank row
   next to the Daldykan river, hence the products without it.
   0. Maps, 3x3 means (each pair is coregistered to its own reference date,
      so a point target such as a tank hops by one pixel between products):
      amplitude normalised by the pair's median (HyP3 *_amp differs between
      pairs with the same reference date), mean before the collapse and dB
      mean per snow-free summer; coherence of the winter pairs (both dates
      Nov-Apr) before the collapse. Stable targets: winter coherence >=
      COH_STABLE (snow-covered tundra decorrelates, structures and rock do
      not; the 5-look estimator floor is ~0.45). Bright: amplitude >=
      BRIGHT_PCT percentile of the crop.
   1. Tank location: within SEARCH_M of the published position, the bright
      pixel whose amplitude drops most from the snow-free summer 2019 to the
      snow-free summer 2020 (the tank shell was destroyed), set against the
      normal change from 2018 to 2019.
   2. Relative series, pairs up to the last acquisition before the collapse
      only: per pair, mean LOS of the point and its 3x3 window minus the mean
      of the stable targets in a RING_M ring, only in the same SNAPHU
      component as the point (whole-cycle offsets between components
      otherwise), >= MIN_RING_PS of them. Structure versus structures, so
      winter pairs count. Vertical = LOS / cos(incidence), + up.
   3. Segments: the tank decorrelates at snowmelt (Apr-Jun) and early winter
      (Nov-Dec), so the pair graph splits into separate connected sets of
      dates; a small SBAS on each, nothing bridges the gaps (no multi-year
      trend). The published point's segments with >= MIN_SEG_DATES dates are
      the observable intervals; the change over each is compared with
      Zhang's rate (ZHANG_RATE x duration).
   4. Null: the same for N_NULL random bright pixels 1.5-5.5 km away whose
      winter coherence is within COH_MATCH of the tank's, so the error bar
      contains atmosphere, snow and unwrapping as they act on comparable
      targets here.
   5. Snow-free check: per summer 2017-2019, change of the relative series
      between the first and last snow-free date, from snow-free pairs only.
   6. Coherence precursor: per snow season (both dates Nov-Apr), mean
      coherence of the 3x3 window minus that of the ring's stable targets.
      A drop in 2019-20 alone would mean the target changed before the
      collapse (settlement, but also roof, fuel level or snow on it).
B. Summers on the 40 m stacks of step 09 (run it with --conncomp reference
   for every year): TPP-3 polygon and tank window vs their ring, seasonal
   sqrt(TDD) amplitude of LOS / cos(inc) (the vertical stack masks slopes
   > 5 deg, and the surface DEM makes tanks look steep), with the
   spatial-shift null of step 16.

Outputs: outputs/tank_series.json, outputs/tank_series.csv,
outputs/fig_tank_series.png

    python 18_tank_series.py
"""
import importlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from rasterio.features import rasterize
from rasterio.warp import Resampling, transform_geom
from scipy.ndimage import binary_dilation, uniform_filter

from common import grid, pair_dates, save_json, to_grid
from config import DIRS, OUTPUTS, POINTS, S1_WAVELENGTH

ts = importlib.import_module("09_timeseries")
loc = importlib.import_module("14_local_signal")
rob = importlib.import_module("16_robustness")

PROD = DIRS["hyp3"] / "norilsk_tank_nowm"
COLLAPSE = "20200529"
ZHANG_RATE = -35.0        # mm/yr at the site, Zhang et al. 2022 (vs SE hills)
TANK_PT = "tank5_zhang2022"
SEARCH_M = 300
N_LOSS, LOSS_SEP_M = 3, 100   # separate amplitude-loss maxima tested as tank candidates
RING_M = (300, 1500)
NULL_KM = (1.5, 5.5)
N_NULL = 200
COH_MIN = 0.3             # per pair, pixels below are not used
COH_STABLE = 0.6          # winter coherence of a stable target (ring)
COH_MATCH = 0.05          # null points: winter coherence within this of the tank's
BRIGHT_PCT = 90           # bright: amplitude percentile within the crop
MIN_RING_PS = 10
SNOW_FREE = ("06-15", "09-30")
WINTER_MONTHS = ("11", "12", "01", "02", "03", "04")
SUMMERS = (2017, 2018, 2019)
WINTERS = ("2017-18", "2018-19", "2019-20")
MIN_SEG_DATES = 5
rng = np.random.default_rng(11)


# ----------------------------------------------------------------- loading
def product_pairs(prod=PROD):
    out = []
    for unw in sorted(prod.rglob("*_unw_phase.tif")):
        r, s = pair_dates(unw.name)
        out.append((r, s, str(unw).removesuffix("_unw_phase.tif")))
    return out


def reference_grid(stem):
    with rasterio.open(f"{stem}_unw_phase.tif") as s:
        return dict(height=s.height, width=s.width, transform=s.transform, crs=s.crs)


def read_layer(path, g, nearest=False):
    with rasterio.open(path) as s:
        a = s.read(1).astype("float32")
        if s.nodata is not None:
            a[a == s.nodata] = np.nan
        if s.shape == (g["height"], g["width"]) and s.transform == g["transform"]:
            return a
        return to_grid(a, s.transform, s.crs, g,
                       Resampling.nearest if nearest else Resampling.bilinear)


def load(stem, g):
    unw = read_layer(f"{stem}_unw_phase.tif", g)
    unw[unw == 0] = np.nan
    cc = read_layer(f"{stem}_conncomp.tif", g, nearest=True)
    coh = read_layer(f"{stem}_corr.tif", g)
    unw[(cc == 0) | ~np.isfinite(cc)] = np.nan
    return -unw * S1_WAVELENGTH / (4 * np.pi) * 1000, cc, coh     # LOS mm, + toward


def snow_free(d):
    md = f"{d[4:6]}-{d[6:]}"
    return SNOW_FREE[0] <= md <= SNOW_FREE[1]


def last_before(pairs):
    return max(d for p in pairs for d in p[:2] if d < COLLAPSE)


# ----------------------------------------------------------------- geometry
def rc_of(lat, lon, g):
    x, y = Transformer.from_crs(4326, g["crs"], always_xy=True).transform(lon, lat)
    c, r = ~g["transform"] * (x, y)
    return int(r), int(c)


def latlon_of(r, c, g):
    x, y = g["transform"] * (c + 0.5, r + 0.5)
    lon, lat = Transformer.from_crs(g["crs"], 4326, always_xy=True).transform(x, y)
    return round(float(lat), 5), round(float(lon), 5)


def distance_m(r, c, g, px):
    rr, cc = np.mgrid[0:g["height"], 0:g["width"]]
    return np.hypot(rr - r, cc - c) * px


def window_ring(r, c, g, px, stable):
    """Flat indices: the 3x3 window on the point; stable targets in the
    RING_M ring."""
    d = distance_m(r, c, g, px)
    rr, cc = np.mgrid[0:g["height"], 0:g["width"]]
    win = (np.abs(rr - r) <= 1) & (np.abs(cc - c) <= 1)
    ring = (d >= RING_M[0]) & (d <= RING_M[1]) & stable
    return np.flatnonzero(win), np.flatnonzero(ring)


# ----------------------------------------------------------------- steps 0-1
def site_maps(pairs, g):
    """3x3 means: normalised amplitude before the collapse, dB amplitude of
    the snow-free pairs of each summer, coherence of the winter pairs before
    the collapse."""
    last = last_before(pairs)
    winter = lambda d: d[4:6] in WINTER_MONTHS
    acc = {}

    def add(key, a):
        s, n = acc.get(key, (0.0, 0))
        acc[key] = (s + a, n + 1)

    for r, s, stem in pairs:
        pre = s <= last
        year = r[:4] if r[:4] == s[:4] and snow_free(r) and snow_free(s) else None
        if pre and winter(r) and winter(s):
            add("winter_coh", uniform_filter(np.nan_to_num(read_layer(f"{stem}_corr.tif", g)), 3))
        if not (pre or year):
            continue
        a = read_layer(f"{stem}_amp.tif", g).astype("float64")
        a = uniform_filter(np.nan_to_num(a / np.nanmedian(a)), 3)
        a[a <= 0] = np.nan
        if pre:
            add("amp", a)
        if year:
            add(year, 10 * np.log10(a))
    return {k: s / n for k, (s, n) in acc.items()}


def locate_tank(maps, bright, g, px):
    """Up to N_LOSS separate (> LOSS_SEP_M apart) amplitude-loss maxima among
    bright pixels within SEARCH_M of the published position."""
    r0, c0 = rc_of(*POINTS[TANK_PT], g)
    change = maps["2020"] - maps["2019"]
    normal = maps["2019"] - maps["2018"]
    ok = bright & np.isfinite(change)
    score = np.where((distance_m(r0, c0, g, px) <= SEARCH_M) & ok, change, np.inf)
    found, info = [], {"published": list(POINTS[TANK_PT]),
                       "published_winter_coh": round(float(maps["winter_coh"][r0, c0]), 2),
                       "bright_px_change_2019_2020_p5_p95_db":
                           np.round(np.percentile(change[ok], [5, 95]), 1).tolist()}
    for k in range(N_LOSS):
        r, c = np.unravel_index(np.argmin(score), score.shape)
        if not np.isfinite(score[r, c]):
            break
        found.append((int(r), int(c)))
        info[f"loss_{k + 1}"] = {
            "latlon": list(latlon_of(r, c, g)), "row_col": [int(r), int(c)],
            "offset_m": round(float(np.hypot(r - r0, c - c0) * px)),
            "amp_change_2019_2020_db": round(float(change[r, c]), 1),
            "amp_change_2018_2019_db": round(float(normal[r, c]), 1),
            "percentile_among_bright_px": round(float(np.mean(change[ok] < change[r, c]) * 100), 2),
            "winter_coh": round(float(maps["winter_coh"][r, c]), 2)}
        score[distance_m(r, c, g, px) <= LOSS_SEP_M] = np.inf
    return found, (r0, c0), info, change


# ----------------------------------------------------------------- steps 2-5
def relative_values(pairs, g, points, px, stable):
    """V[k, j]: point-minus-ring LOS mm for pair k and point j (NaN if
    unusable); C[k, j]: point-minus-ring coherence."""
    sets = [window_ring(r, c, g, px, stable) for r, c in points]
    V = np.full((len(pairs), len(points)), np.nan, "float32")
    C = np.full_like(V, np.nan)
    inc = []
    for k, (r, s, stem) in enumerate(pairs):
        b, cc, coh = load(stem, g)
        ok = np.isfinite(b) & (coh >= COH_MIN)
        bf, ccf, okf, cohf = b.ravel(), cc.ravel(), ok.ravel(), coh.ravel()
        for j, (win, ring) in enumerate(sets):
            C[k, j] = np.nanmean(cohf[win]) - np.nanmean(cohf[ring])
            w = win[okf[win]]
            if w.size == 0:
                continue
            lab = np.bincount(ccf[w].astype(int)).argmax()
            rg = ring[okf[ring] & (ccf[ring] == lab)]
            if rg.size < MIN_RING_PS:
                continue
            w = w[ccf[w] == lab]
            V[k, j] = bf[w].mean() - bf[rg].mean()
        if k % 10 == 0:
            with rasterio.open(f"{stem}_lv_theta.tif") as t:
                inc.append(90 - np.degrees(np.nanmedian(t.read(1))))
    return V, C, float(np.median(inc))


def invert_segments(pairs, v):
    """SBAS separately on each connected set (>= 3 dates) of the pairs with a
    finite value: [(dates, values, rms)], each segment starting at 0. Values
    of different segments are not comparable (the gap is not observed)."""
    use = [(p[0], p[1], x) for p, x in zip(pairs, v) if np.isfinite(x)]
    segs = []
    for dates in ts.date_components([(r, s) for r, s, _ in use]):
        if len(dates) < 3:
            continue
        sub = [(r, s, x) for r, s, x in use if r in dates]
        A = ts.design(dates, [(r, s) for r, s, _ in sub])
        obs = np.array([x for *_, x in sub])
        sol, *_ = np.linalg.lstsq(A, obs, rcond=None)
        rms = float(np.sqrt(np.mean((A @ sol - obs) ** 2)))
        segs.append((dates, np.r_[0.0, sol], rms))
    return sorted(segs, key=lambda s: s[0][0])


def interval_change(segs, d0, d1):
    """x(d1) - x(d0) if both dates are in the same segment, else NaN."""
    for dates, x, _ in segs:
        if d0 in dates and d1 in dates:
            return float(x[dates.index(d1)] - x[dates.index(d0)])
    return np.nan


def summer_change(pairs, v, year):
    """Change of the relative series over one snow-free season, from snow-free
    pairs inside that summer only (no snow, no winter bridging)."""
    sub = [(p, x) for p, x in zip(pairs, v)
           if p[0][:4] == p[1][:4] == str(year) and snow_free(p[0]) and snow_free(p[1])]
    segs = invert_segments([p for p, _ in sub], np.array([x for _, x in sub]))
    if not segs:
        return np.nan
    dates, x, _ = max(segs, key=lambda s: len(s[0]))
    return float(x[-1] - x[0])


def winter_of(r, s):
    """'2018-19' if both dates are in the same snow season, else None."""
    if not (r[4:6] in WINTER_MONTHS and s[4:6] in WINTER_MONTHS):
        return None
    y = [int(d[:4]) - (d[4:6] < "11") for d in (r, s)]
    return f"{y[0]}-{str(y[0] + 1)[2:]}" if y[0] == y[1] else None


def point_stats(pairs, v, c):
    """Segments and fixed metrics of one point; v: vertical relative values
    per pair, c: relative coherence per pair."""
    out = {"segments": invert_segments(pairs, v)}
    for y in SUMMERS:
        out[f"summer_{y}_mm"] = summer_change(pairs, v, y)
    season = np.array([winter_of(r, s) for r, s, _ in pairs])
    for w in WINTERS:
        out[f"winter_coh_{w}"] = float(np.nanmean(c[season == w]))
    return out


def add_intervals(stats, intervals):
    """Change over each (d0, d1) interval, for every point that observed it."""
    for s in stats:
        for d0, d1 in intervals:
            s[f"change_{d0}_{d1}_mm"] = interval_change(s["segments"], d0, d1)


def null_points(maps, bright, g, points, pub, px):
    """Bright pixels NULL_KM from the published point whose winter coherence
    is within COH_MATCH of the range of the tested points."""
    d = distance_m(*pub, g, px) / 1000
    wc = maps["winter_coh"]
    lo, hi = min(wc[p] for p in points), max(wc[p] for p in points)
    like = (wc >= lo - COH_MATCH) & (wc <= hi + COH_MATCH)
    cand = np.flatnonzero(((d >= NULL_KM[0]) & (d <= NULL_KM[1]) & bright & like).ravel())
    pick = rng.choice(cand, size=min(N_NULL, cand.size), replace=False)
    return [tuple(int(i) for i in np.unravel_index(k, wc.shape)) for k in pick]


def against_null(obs, nulls, metrics):
    out = {}
    for m in metrics:
        nv = np.array([n[m] for n in nulls])
        nv = nv[np.isfinite(nv)]
        if nv.size < 20 or not np.isfinite(obs[m]):
            out[m] = None
            continue
        dev = np.abs(nv - nv.mean())
        nd = 3 if m.startswith("winter_coh") else 1
        out[m] = {"value": round(obs[m], nd), "null_mean": round(float(nv.mean()), nd),
                  "null_std": round(float(nv.std()), nd), "n_null": int(nv.size),
                  "z": round(float((obs[m] - nv.mean()) / nv.std()), 2),
                  "p_two_sided": round(float(np.mean(dev >= abs(obs[m] - nv.mean()))), 3)}
    return out


def series_part(rep):
    pairs = product_pairs()
    if not pairs:
        print(f"No products in {PROD}: run step 02 with --name norilsk_tank_nowm (part A skipped).")
        return None
    g = reference_grid(pairs[0][2])
    px = abs(g["transform"].a)
    maps = site_maps(pairs, g)
    bright = maps["amp"] >= np.nanpercentile(maps["amp"], BRIGHT_PCT)
    stable = maps["winter_coh"] >= COH_STABLE
    found, pub, rep["location"], change = locate_tank(maps, bright, g, px)
    rep["location"]["stable_px_fraction"] = round(float(stable.mean()), 3)
    print("tank location:", rep["location"])

    last = last_before(pairs)
    pre = [p for p in pairs if p[1] <= last]
    tested = found + [pub]
    names = [f"loss_{k + 1}" for k in range(len(found))] + ["published"]
    nulls = null_points(maps, bright, g, tested, pub, px)
    V, C, inc = relative_values(pre, g, tested + nulls, px, stable)
    V /= np.cos(np.radians(inc))
    stats = [point_stats(pre, V[:, j], C[:, j]) for j in range(V.shape[1])]
    # Observable intervals: the published point's segments (gaps = melt / early
    # winter, where the tank decorrelates and nothing bridges them).
    intervals = [(d[0], d[-1]) for d, _, _ in stats[len(found)]["segments"]
                 if len(d) >= MIN_SEG_DATES]
    add_intervals(stats, intervals)
    metrics = ([f"change_{a}_{b}_mm" for a, b in intervals]
               + [f"summer_{y}_mm" for y in SUMMERS] + [f"winter_coh_{w}" for w in WINTERS])
    null_stats = stats[len(tested):]
    rep["incidence_deg"] = round(inc, 1)
    rep["last_date_before"] = last
    rep["n_pairs_before"] = len(pre)
    rep["n_null_points"] = len(nulls)
    rep["intervals"] = [{"from": a, "to": b, "days": (pd.Timestamp(b) - pd.Timestamp(a)).days,
                         "zhang_rate_expected_mm": round(ZHANG_RATE * (pd.Timestamp(b) - pd.Timestamp(a)).days / 365.25, 1)}
                        for a, b in intervals]
    rep["points"] = {}
    for name, s in zip(names, stats):
        rep["points"][name] = {
            "segments": [{"from": d[0], "to": d[-1], "n_dates": len(d), "rms_mm": round(rms, 1)}
                         for d, _, rms in s["segments"]],
            "vs_null": against_null(s, null_stats, metrics)}
    pd.DataFrame([{"point": n, "segment": k, "date": d, "minus_ring_vert_mm": round(float(v), 1)}
                  for n, s in zip(names, stats)
                  for k, (dates, x, _) in enumerate(s["segments"]) for d, v in zip(dates, x)]
                 ).to_csv(OUTPUTS / "tank_series.csv", index=False)
    return dict(stats=dict(zip(names, stats)), null_stats=null_stats, intervals=intervals,
                change=change, stable=stable, found=found, pub=pub, px=px)


# ----------------------------------------------------------------- part B
def summer_site_tests():
    """TPP-3 polygon and tank window vs ring on the 40 m summer stacks."""
    prof = grid()
    gj = json.loads((DIRS["osm"] / "osm_landuse.geojson").read_text())
    f = [f for f in gj["features"] if f["properties"]["osm_id"] == "way/525598517"][0]
    plant = rasterize([(transform_geom("EPSG:4326", prof["crs"], f["geometry"]), 1)],
                      out_shape=(prof["height"], prof["width"]), transform=prof["transform"],
                      fill=0, dtype="uint8").astype(bool)
    r, c = rc_of(*POINTS[TANK_PT], prof)
    tank = np.zeros_like(plant)
    tank[r - 1:r + 2, c - 1:c + 2] = True
    rin, rout = int(RING_M[0] / loc.PX), int(RING_M[1] / loc.PX)
    shifts = rob.random_shifts(400)
    out = {}
    # LOS / cos(inc), not the vertical stack: Copernicus DEM is a surface model,
    # so tanks and buildings look like steep slopes and are masked there.
    for f in sorted(DIRS["ts"].glob("s1_*_los.tif")):
        year = f.stem.split("_")[1]
        inc = json.loads((DIRS["ts"] / f"s1_{year}_meta.json").read_text())["incidence_median_deg"]
        with rasterio.open(f) as s:
            st = s.read().astype("float32") / np.cos(np.radians(inc))
            dates = list(s.descriptions)
        seas = loc.seasonal_map(st, loc.sqrt_tdd(dates))
        out[year] = {"dates": [dates[0], dates[-1]]}
        for name, m in (("tpp3_plant", plant), ("tank_window", tank)):
            ring = binary_dilation(m, iterations=rout) & ~binary_dilation(m, iterations=rin)
            iy, ix = np.nonzero(m)
            ry, rx = np.nonzero(ring)
            d = lambda dy, dx: (np.nanmean(rob.take(seas, iy, ix, dy, dx))
                                - np.nanmean(rob.take(seas, ry, rx, dy, dx)))
            obs = d(0, 0)
            null = np.array([d(dy, dx) for dy, dx in shifts])
            null = null[np.isfinite(null)]
            out[year][name] = {
                "inside_minus_ring_mm": None if not np.isfinite(obs) else round(float(obs), 1),
                "coverage": round(float(np.isfinite(seas[m]).mean()), 2),
                "null_std": round(float(null.std()), 1),
                "z": None if not np.isfinite(obs) else round(float((obs - null.mean()) / null.std()), 2)}
    return out


# ----------------------------------------------------------------- main
def main():
    rep = {}
    a = series_part(rep)
    rep["summers_40m"] = summer_site_tests()
    save_json(OUTPUTS / "tank_series.json", rep)
    print(json.dumps(rep, indent=1, default=str))
    if a:
        figure(rep, a)


def figure(rep, a):
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    colors = {name: f"C{i}" for i, name in enumerate(a["stats"])}
    for name, s in a["stats"].items():
        for k, (dates, x, _) in enumerate(s["segments"]):
            axs[0].plot(pd.to_datetime(dates), x, ".-", ms=3, lw=0.7, color=colors[name],
                        label=name if k == 0 else None)
    for k, (d0, d1) in enumerate(a["intervals"]):
        t = pd.to_datetime([d0, d1])
        axs[0].plot(t, [0, ZHANG_RATE * (t[1] - t[0]).days / 365.25], "k:", lw=1,
                    label="Zhang 2022 rate, -35 mm/yr" if k == 0 else None)
    axs[0].axvline(pd.Timestamp(COLLAPSE), color="red", lw=1, label="collapse 2020-05-29")
    axs[0].set_ylabel("minus stable targets 0.3-1.5 km, vertical [mm], + up")
    axs[0].legend(fontsize=7)
    axs[0].set_title("relative series before the collapse (each segment starts at 0; "
                     "gaps = snowmelt / early winter, not observable)", fontsize=8)
    keys = [f"change_{d0}_{d1}_mm" for d0, d1 in a["intervals"]]
    nv = [[v for s in a["null_stats"] if np.isfinite(v := s[m])] for m in keys]
    axs[1].boxplot(nv, whis=(5, 95), showfliers=False)
    for j, (name, s) in enumerate(a["stats"].items()):
        axs[1].plot(np.arange(1, len(keys) + 1) + 0.08 * (j - 1.5), [s[m] for m in keys], "o",
                    color=colors[name], ms=5, label=name)
    zr = [ZHANG_RATE * (pd.Timestamp(d1) - pd.Timestamp(d0)).days / 365.25 for d0, d1 in a["intervals"]]
    axs[1].plot(np.arange(1, len(keys) + 1), zr, "k_", ms=20, mew=2, label="Zhang rate")
    axs[1].set_xticks(np.arange(1, len(keys) + 1),
                      [f"{d0[2:4]}/{d0[4:6]}-{d1[2:4]}/{d1[4:6]}" for d0, d1 in a["intervals"]], fontsize=8)
    axs[1].axhline(0, color="grey", lw=0.5)
    axs[1].set_ylabel("change over the interval [mm], + up")
    axs[1].set_title("tank vs comparable targets 1.5-5.5 km (box 25-75%, whiskers 5-95%)", fontsize=9)
    axs[1].legend(fontsize=7)
    r0, c0 = a["pub"]
    h = int(1000 / a["px"])
    win = slice(r0 - h, r0 + h + 1), slice(c0 - h, c0 + h + 1)
    im = axs[2].imshow(a["change"][win], cmap="RdBu", vmin=-12, vmax=12)
    axs[2].contour(a["stable"][win], levels=[0.5], colors="black", linewidths=0.3)
    for name, (r, c) in zip(a["stats"], a["found"] + [a["pub"]]):
        axs[2].plot(c - c0 + h, r - r0 + h, "x" if name == "published" else "+",
                    color=colors[name], ms=10, mew=2, label=name)
    axs[2].set_title("amplitude change, summer 2019 -> 2020 [dB] (2 km; black = stable targets)",
                     fontsize=9)
    axs[2].legend(fontsize=7)
    fig.colorbar(im, ax=axs[2], shrink=0.7)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_tank_series.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
