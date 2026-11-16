"""Step 14: LOCAL seasonal deformation, and is it the same in both summers?

Why: zone means at 10-50 km scale are dominated by turbulent troposphere,
which is independent from one acquisition to the next. Ground deformation
that depends on the ground (ice content, wetness, foundations) should repeat
in the same place every summer. So:

  1. per summer, per pixel: fit disp(t) = a + b*sqrt(TDD(t)) on the vertical
     stack of step 09; seasonal map = b * (sqrt(TDD_end) - sqrt(TDD_start)), mm
  2. local map = seasonal map - its Gaussian-smoothed version (sigma HP_KM):
     what the pixel does relative to its own ~2 km neighbourhood. This removes
     the long-wavelength atmosphere and any residual ramp.
  3. reproducibility: Pearson r between the 2020 and 2025 local maps
     (decimated every DEC px ~ independent samples). Null: the same r with the
     2025 map shifted by >= 6 km in many directions (spatial structure kept,
     co-location broken).
  4. persistent map: mean of the two summers where both have the same sign;
     hotspots = pixels below HOT_MM in BOTH summers, grouped into clusters.
  5. per OSM facility (industrial / quarry / landfill, named where possible):
     mean seasonal displacement inside the polygon minus a 0.3-1.5 km ring
     around it, per summer. 'persistent' = same sign and > 2 sigma in both.
  6. by land cover (ESA WorldCover): mean local signal per class and summer.

Outputs: data/ts/local_<year>.tif, data/ts/local_persistent.tif,
outputs/local_meta.json, outputs/local_landcover.csv,
outputs/local_facilities.csv, outputs/local_hotspots.csv,
outputs/fig_local.png

    python 14_local_signal.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from rasterio.warp import transform, transform_geom
from scipy.ndimage import binary_dilation, gaussian_filter, label

from common import grid, read, save_json, write_stack, zones
from config import DIRS, OUTPUTS

YEARS = (2020, 2025)
HP_KM = 2.0
DEC = 5                  # decimation for correlation (150 m)
HOT_MM = -10.0
PX = 30.0
RING = (300, 1500)       # m, facility comparison ring
WC = {10: "tree cover", 20: "shrubland", 30: "grassland", 50: "built-up",
      60: "bare / sparse", 80: "water edge", 90: "herbaceous wetland",
      100: "moss and lichen"}


def nangauss(a, sigma_px):
    m = np.isfinite(a)
    num = gaussian_filter(np.where(m, a, 0.0), sigma_px)
    den = gaussian_filter(m.astype("float32"), sigma_px)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0.2, num / den, np.nan)


def sqrt_tdd(dates):
    s = pd.read_csv(DIRS["climate"] / "era5land_daily_norilsk_centre.csv",
                    index_col="date", parse_dates=["date"])["tdd_cum"]
    return np.sqrt(np.array([s[pd.Timestamp(d)] for d in dates]))


def seasonal_map(stack, x):
    """Per-pixel OLS slope vs x, times the summer's x range (mm). Pixels with
    more than one missing date are NaN."""
    m = np.isfinite(stack)
    X = np.where(m, x[:, None, None], np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        xm, ym = np.nanmean(X, 0), np.nanmean(stack, 0)
        b = np.nansum((X - xm) * (stack - ym), 0) / np.nansum((X - xm) ** 2, 0)
    b[m.sum(0) < len(x) - 1] = np.nan
    return (b * (x[-1] - x[0])).astype("float32")


def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    sub = np.zeros_like(m)
    sub[::DEC, ::DEC] = True
    m &= sub
    return float(np.corrcoef(a[m], b[m])[0, 1]), int(m.sum())


def shift_null(a, b, min_km=6, n=24):
    """r between a and b shifted by >= min_km in n directions."""
    d = int(min_km * 1000 / PX)
    rs = []
    for k in range(n):
        ang = 2 * np.pi * k / n
        dy, dx = int(round(d * np.sin(ang))), int(round(d * np.cos(ang)))
        bs = np.full_like(b, np.nan)
        H, W = b.shape
        ys, xs = slice(max(dy, 0), H + min(dy, 0)), slice(max(dx, 0), W + min(dx, 0))
        yd, xd = slice(max(-dy, 0), H + min(-dy, 0)), slice(max(-dx, 0), W + min(-dx, 0))
        bs[ys, xs] = b[yd, xd]
        rs.append(corr(a, bs)[0])
    return float(np.mean(rs)), float(np.std(rs)), float(np.max(np.abs(rs)))


def se_mean(v):
    v = v[np.isfinite(v)]
    return float(np.std(v) / np.sqrt(max(v.size / 25, 1))) if v.size else np.nan


def facilities(seas, prof):
    gj = json.loads((DIRS["osm"] / "osm_landuse.geojson").read_text())
    rows = []
    ring_out = int(RING[1] / PX)
    ring_in = int(RING[0] / PX)
    for f in gj["features"]:
        geom = transform_geom("EPSG:4326", prof["crs"], f["geometry"])
        m = rasterize([(geom, 1)], out_shape=(prof["height"], prof["width"]),
                      transform=prof["transform"], fill=0, dtype="uint8").astype(bool)
        if m.sum() * PX * PX / 1e6 < 0.05:
            continue
        near = binary_dilation(m, iterations=ring_in)
        ring = binary_dilation(m, iterations=ring_out) & ~near
        r = {"osm_id": f["properties"]["osm_id"], "landuse": f["properties"]["landuse"],
             "name": f["properties"].get("name"), "area_km2": round(m.sum() * PX * PX / 1e6, 2)}
        ys, xs = np.nonzero(m)
        lon, lat = transform(prof["crs"], "EPSG:4326",
                             [prof["transform"].c + (xs.mean() + .5) * PX],
                             [prof["transform"].f - (ys.mean() + .5) * PX])
        r["lon"], r["lat"] = round(lon[0], 4), round(lat[0], 4)
        for y, s in seas.items():
            inside, outside = s[m], s[ring]
            cov = np.isfinite(inside).mean()
            if cov < 0.2 or np.isfinite(outside).sum() < 50:
                r[f"d_{y}"], r[f"se_{y}"], r[f"cov_{y}"] = np.nan, np.nan, round(float(cov), 2)
                continue
            r[f"d_{y}"] = round(float(np.nanmean(inside) - np.nanmean(outside)), 2)
            r[f"se_{y}"] = round(float(np.hypot(se_mean(inside), se_mean(outside))), 2)
            r[f"cov_{y}"] = round(float(cov), 2)
        rows.append(r)
    df = pd.DataFrame(rows)
    ys = [str(y) for y in seas]
    ok = df[[f"d_{y}" for y in ys]].notna().all(axis=1)
    sig = lambda y: (df[f"d_{y}"].abs() > 2 * df[f"se_{y}"])
    same = np.sign(df[f"d_{ys[0]}"]) == np.sign(df[f"d_{ys[1]}"])
    df["persistent"] = ok & same & sig(ys[0]) & sig(ys[1])
    df["d_mean"] = df[[f"d_{y}" for y in ys]].mean(axis=1).round(2)
    return df.sort_values("d_mean")


def hotspots(pers, both, prof, lc):
    lab, n = label(both)
    rows = []
    for i in range(1, n + 1):
        m = lab == i
        if m.sum() < 20:          # < 0.018 km2
            continue
        ys, xs = np.nonzero(m)
        lon, lat = transform(prof["crs"], "EPSG:4326",
                             [prof["transform"].c + (xs.mean() + .5) * PX],
                             [prof["transform"].f - (ys.mean() + .5) * PX])
        cls = np.bincount(lc[m].astype(int), minlength=101).argmax()
        rows.append({"lon": round(lon[0], 4), "lat": round(lat[0], 4),
                     "area_km2": round(m.sum() * PX * PX / 1e6, 3),
                     "mean_mm": round(float(np.nanmean(pers[m])), 1),
                     "landcover": WC.get(int(cls), str(cls))})
    return pd.DataFrame(rows).sort_values("area_km2", ascending=False)


def main():
    prof = grid()
    z, _, _ = zones()
    lc = read(DIRS["landcover"] / "worldcover_utm30.tif")[0]
    seas, local, series = {}, {}, {}
    for y in YEARS:
        with rasterio.open(DIRS["ts"] / f"s1_{y}_vert.tif") as s:
            st = s.read().astype("float32")
            dates = list(s.descriptions)
        x = sqrt_tdd(dates)
        seas[y] = seasonal_map(st, x)
        local[y] = seas[y] - nangauss(seas[y], HP_KM * 1000 / PX)
        # per-date local series by land cover (each date high-passed)
        hp = np.array([b - nangauss(b, HP_KM * 1000 / PX) for b in st])
        series[y] = {"dates": dates, "sqrt_tdd": x.round(2).tolist(),
                     "by_class": {WC[c]: [round(float(np.nanmean(h[lc == c])), 2) for h in hp]
                                  for c in WC if (lc == c).sum() > 2000}}
        write_stack(DIRS["ts"] / f"local_{y}.tif", [local[y]], prof, [f"local seasonal {y}"])

    a, b = local[YEARS[0]], local[YEARS[1]]
    r, n = corr(a, b)
    null_mean, null_std, null_max = shift_null(a, b)
    r_raw, _ = corr(seas[YEARS[0]], seas[YEARS[1]])
    dem = read(DIRS["dem"] / "dem_utm30.tif")[0]
    hd = dem - nangauss(dem, HP_KM * 1000 / PX)
    r_dem = [corr(local[y], hd)[0] for y in YEARS]

    same = np.isfinite(a) & np.isfinite(b) & (np.sign(a) == np.sign(b))
    pers = np.where(same, (a + b) / 2, np.nan).astype("float32")
    both = np.isfinite(a) & np.isfinite(b) & (a < HOT_MM) & (b < HOT_MM)
    valid = np.isfinite(a) & np.isfinite(b)
    p_a, p_b = (a[valid] < HOT_MM).mean(), (b[valid] < HOT_MM).mean()
    write_stack(DIRS["ts"] / "local_persistent.tif", [pers], prof, ["persistent local seasonal"])

    lcrows = []
    for c, name in WC.items():
        m = valid & (lc == c)
        if m.sum() < 2000:
            continue
        lcrows.append({"class": name, "code": c, "area_km2": round(m.sum() * PX * PX / 1e6, 1),
                       **{f"local_{y}_mm": round(float(np.nanmean(local[y][m])), 2) for y in YEARS},
                       **{f"se_{y}": round(se_mean(local[y][m]), 2) for y in YEARS},
                       "frac_hotspot": round(float(both[m].mean()), 4)})
    lcdf = pd.DataFrame(lcrows)
    lcdf.to_csv(OUTPUTS / "local_landcover.csv", index=False)
    fac = facilities(seas, prof)
    fac.to_csv(OUTPUTS / "local_facilities.csv", index=False)
    hot = hotspots(pers, both, prof, lc)
    hot.to_csv(OUTPUTS / "local_hotspots.csv", index=False)

    meta = {"hp_km": HP_KM, "hot_mm": HOT_MM, "decimation_px": DEC,
            "r_local": round(r, 3), "n_independent": n,
            "r_null_shift_mean": round(null_mean, 3), "r_null_shift_std": round(null_std, 3),
            "r_null_shift_max_abs": round(null_max, 3),
            "r_raw_seasonal": round(r_raw, 3),
            "r_local_vs_dem_relief": [round(v, 3) for v in r_dem],
            "hotspot_px_both": int(both.sum()),
            "hotspot_px_expected_if_independent": int(p_a * p_b * valid.sum()),
            "hotspot_area_km2": round(both.sum() * PX * PX / 1e6, 1),
            "valid_area_km2": round(valid.sum() * PX * PX / 1e6, 1),
            "n_facilities": int(len(fac)), "n_facilities_persistent": int(fac.persistent.sum()),
            "series": series}
    save_json(OUTPUTS / "local_meta.json", meta)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5.5))
    for ax, y in zip(axs[:2], YEARS):
        im = ax.imshow(local[y], cmap="RdBu", vmin=-20, vmax=20)
        ax.set_title(f"local seasonal {y} [mm] (HP {HP_KM:g} km)")
        ax.set_xticks([]); ax.set_yticks([])
    axs[2].hexbin(a[::DEC, ::DEC].ravel(), b[::DEC, ::DEC].ravel(), gridsize=60,
                  extent=(-40, 40, -40, 40), bins="log", cmap="Greys")
    axs[2].set_xlabel("2020 local [mm]"); axs[2].set_ylabel("2025 local [mm]")
    axs[2].set_title(f"r = {r:.2f} (shift null {null_mean:.2f} ± {null_std:.2f})")
    fig.colorbar(im, ax=axs[:2], shrink=0.7)
    fig.savefig(OUTPUTS / "fig_local.png", dpi=110)
    plt.close(fig)

    print(f"local r(2020, 2025) = {r:.3f} (n={n}); shift null {null_mean:.3f} ± {null_std:.3f} "
          f"(max |r| {null_max:.3f}); raw seasonal r = {r_raw:.3f}; r vs relief {r_dem}")
    print(f"hotspots (< {HOT_MM} mm both summers): {both.sum()} px, expected {int(p_a * p_b * valid.sum())}")
    print(lcdf.to_string(index=False))
    print(fac.head(12).to_string(index=False))
    print(f"persistent facilities: {fac.persistent.sum()} / {len(fac)}")
    print(hot.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
