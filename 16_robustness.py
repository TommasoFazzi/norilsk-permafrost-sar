"""Step 16: how robust are the local (step 14) and zone (step 11) results?

Variants of the Sentinel-1 series (same pairs and inversion as step 09):
  v0       conncomp any,       troposphere height + ramp   (results v0)
  v1       conncomp reference, troposphere height + ramp   (main from v1)
  v1_h     conncomp reference, troposphere height only
  v1_none  conncomp reference, no troposphere model
For each: wetland local signal (and wetland minus moss/lichen), replication r
of the local maps between summers, paired zone slopes (industrial - tundra
controls, tundra NW - SE).

Null tests on v1. The maps are shifted by SHIFT_KM relative to the masks:
this keeps the spatial structure of the noise (atmosphere, residual ramps)
and breaks the co-location, so the spread of the shifted statistic is the
honest error bar, unlike the per-pixel standard errors.
  - land cover: mean local signal inside the shifted wetland mask
  - zones: paired slope (industrial - controls) with the industrial mask
    shifted over lowland ground
  - OSM facilities: number of 'persistent' sites (same sign, > 2 se in both
    summers, step 14 rule) at shifted positions; mean inside-minus-ring
  - replication r for high-pass scales 1, 2 and 4 km
Stefan check: per date, wetland minus moss/lichen of the high-passed map
(both relative to their 2 km surroundings) vs sqrt(TDD); slopes of the two
summers compared.
Independent sensor: the same wetland contrast in the NISAR L-band 2026
vertical maps (step 10), vs the value the Sentinel-1 Stefan slope predicts.

Outputs: outputs/robustness.json, outputs/fig_robustness.png

    python 16_robustness.py [--shifts 400]
"""
import argparse
import importlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import rasterize
from scipy import stats
from rasterio.warp import transform_geom
from scipy.ndimage import binary_dilation

from common import read, save_json
from config import DIRS, OUTPUTS

ts = importlib.import_module("09_timeseries")
loc = importlib.import_module("14_local_signal")
an = importlib.import_module("11_analysis")

YEARS = (2020, 2025)
WET, MOSS = 90, 100
VARIANTS = {"v0": ("any", "height_ramp"), "v1": ("reference", "height_ramp"),
            "v1_h": ("reference", "height"), "v1_none": ("reference", "none")}
SHIFT_KM = (3.0, 20.0)
LOWLAND_MAX_M = 250.0
rng = np.random.default_rng(7)


def local(a, hp_km=loc.HP_KM):
    return a - loc.nangauss(a, hp_km * 1000 / loc.PX)


def build_variants(ctx):
    """{variant: {year: (dates, vert stack)}}; each conncomp mode is loaded once."""
    out = {v: {} for v in VARIANTS}
    for year in YEARS:
        items = sorted(ts.find_pairs()[year])
        for cc in ("any", "reference"):
            b_raw, coh, inc = ts.load_year(items, ctx, cc)
            for v, (vcc, tropo) in VARIANTS.items():
                if vcc != cc:
                    continue
                b, c, pairs, dates, _, _ = ts.preprocess(b_raw, coh, items, ctx, tropo, False, cc)
                _, vert, rms, *_ = ts.solve(dates, pairs, b, c, inc, ctx)
                out[v][year] = (dates, vert)
                print(f"  {v} {year}: {len(dates)} dates, {len(pairs)} pairs, rms {rms:.1f} mm")
            del b_raw
    return out


def seasonal(dates, vert):
    return loc.seasonal_map(vert, loc.sqrt_tdd(dates))


def contrast(lm, lc):
    w = lm[(lc == WET) & np.isfinite(lm)]
    m = lm[(lc == MOSS) & np.isfinite(lm)]
    return {"wetland": round(float(w.mean()), 2), "moss": round(float(m.mean()), 2),
            "wetland_minus_moss": round(float(w.mean() - m.mean()), 2)}


def zone_pair(dates, vert, mask_a, mask_b):
    """Paired slope of (zone a - zone b) vs sqrt(TDD) and end-of-summer difference."""
    x = loc.sqrt_tdd(dates)
    da = np.array([np.nanmean(v[mask_a]) for v in vert])
    db = np.array([np.nanmean(v[mask_b]) for v in vert])
    f = an.fit(x, da - db)
    return {"diff_b": round(float(f["b"]), 3), "se": round(float(f["se_b"]), 3),
            "end_diff_mm": round(float((da - db)[-1]), 1)}


def random_shifts(n):
    d = rng.uniform(*SHIFT_KM, n) * 1000 / loc.PX
    ang = rng.uniform(0, 2 * np.pi, n)
    return np.c_[np.round(d * np.sin(ang)), np.round(d * np.cos(ang))].astype(int)


def take(a, rows, cols, dy, dx):
    """Values of map a at (rows + dy, cols + dx), out-of-grid dropped."""
    r, c = rows + dy, cols + dx
    ok = (r >= 0) & (r < a.shape[0]) & (c >= 0) & (c < a.shape[1])
    return a[r[ok], c[ok]]


def mask_null(stat_map, mask, shifts, min_frac=0.3):
    rows, cols = np.nonzero(mask)
    obs = float(np.nanmean(stat_map[mask]))
    null = []
    for dy, dx in shifts:
        v = take(stat_map, rows, cols, dy, dx)
        if np.isfinite(v).sum() >= min_frac * np.isfinite(stat_map[mask]).sum():
            null.append(np.nanmean(v))
    null = np.array(null)
    return {"obs": round(obs, 2), "null_mean": round(float(null.mean()), 2),
            "null_std": round(float(null.std()), 2), "z_emp": round(float((obs - null.mean()) / null.std()), 1),
            "p_two_sided": float(np.mean(np.abs(null - null.mean()) >= abs(obs - null.mean()))),
            "n_shifts": int(null.size)}


def zone_null(dates, vert, z, lowland, shifts):
    """Paired slope industrial - controls with the industrial mask shifted."""
    ctrl = (z == 5) | (z == 6)
    ind = z == 2
    x = loc.sqrt_tdd(dates)
    dc = np.array([np.nanmean(v[ctrl]) for v in vert])
    rows, cols = np.nonzero(ind)

    def slope(r, c):
        da = np.array([np.nanmean(v[r, c]) for v in vert])
        return an.fit(x, da - dc)["b"]

    obs = slope(rows, cols)
    allowed = lowland & ~ctrl
    null = []
    for dy, dx in shifts:
        r, c = rows + dy, cols + dx
        ok = (r >= 0) & (r < z.shape[0]) & (c >= 0) & (c < z.shape[1])
        r, c = r[ok], c[ok]
        keep = allowed[r, c]
        if keep.sum() < 0.3 * rows.size:
            continue
        null.append(slope(r[keep], c[keep]))
    null = np.array([v for v in null if np.isfinite(v)])
    return {"obs_diff_b": round(float(obs), 3), "null_std": round(float(null.std()), 3),
            "z_emp": round(float((obs - null.mean()) / null.std()), 2), "n_shifts": int(null.size)}


def facility_sets(prof):
    gj = json.loads((DIRS["osm"] / "osm_landuse.geojson").read_text())
    H, W = prof["height"], prof["width"]
    rin, rout = int(loc.RING[0] / loc.PX), int(loc.RING[1] / loc.PX)
    sets = []
    for f in gj["features"]:
        geom = transform_geom("EPSG:4326", prof["crs"], f["geometry"])
        m = rasterize([(geom, 1)], out_shape=(H, W), transform=prof["transform"],
                      fill=0, dtype="uint8").astype(bool)
        if m.sum() * loc.PX ** 2 / 1e6 < 0.05:
            continue
        ys, xs = np.nonzero(m)
        y0, y1 = max(ys.min() - rout - 1, 0), min(ys.max() + rout + 2, H)
        x0, x1 = max(xs.min() - rout - 1, 0), min(xs.max() + rout + 2, W)
        mc = m[y0:y1, x0:x1]
        near = binary_dilation(mc, iterations=rin)
        ring = binary_dilation(mc, iterations=rout) & ~near
        iy, ix = np.nonzero(mc)
        ry, rx = np.nonzero(ring)
        sets.append((iy + y0, ix + x0, ry + y0, rx + x0))
    return sets


def facility_stats(seas, sets, dy=0, dx=0):
    """Count of persistent sites and mean inside-minus-ring (both summers)."""
    n_pers, ds = 0, []
    for iy, ix, ry, rx in sets:
        rows = []
        for s in seas:
            vin, vout = take(s, iy, ix, dy, dx), take(s, ry, rx, dy, dx)
            if vin.size == 0 or np.isfinite(vin).mean() < 0.2 or np.isfinite(vout).sum() < 50:
                rows = None
                break
            d = np.nanmean(vin) - np.nanmean(vout)
            se = np.hypot(loc.se_mean(vin), loc.se_mean(vout))
            rows.append((d, se))
        if rows is None:
            continue
        (d1, s1), (d2, s2) = rows
        ds.append((d1 + d2) / 2)
        if np.sign(d1) == np.sign(d2) and abs(d1) > 2 * s1 and abs(d2) > 2 * s2:
            n_pers += 1
    return n_pers, (float(np.mean(ds)) if ds else np.nan), len(ds)


def stefan(dates, vert, lc):
    """Per-date wetland minus moss/lichen of the high-passed map vs sqrt(TDD)."""
    x = loc.sqrt_tdd(dates)
    y = []
    for v in vert:
        h = local(v)
        y.append(np.nanmean(h[lc == WET]) - np.nanmean(h[lc == MOSS]))
    y = np.array(y) - y[0]
    f = an.fit(x, y)
    return {"dates": dates, "sqrt_tdd": np.round(x, 2).tolist(), "wet_minus_moss_mm": np.round(y, 2).tolist(),
            "b": round(float(f["b"]), 4), "se_b": round(float(f["se_b"]), 4), "r2": round(float(f["r2"]), 3)}


def nisar_wetland(lc, b_stefan, shifts):
    """Independent sensor: NISAR L-band 2026 vertical (asc + desc, step 10).
    Wetland minus moss/lichen of the high-passed 12-day maps, its shift null,
    and the value predicted by the Sentinel-1 Stefan slope over the window."""
    meta_p = DIRS["ts"] / "nisar_meta.json"
    if not meta_p.exists():
        return None
    meta = json.loads(meta_p.read_text())
    tdd = an.tdd_series("norilsk_centre")
    dsq = lambda p: np.sqrt(tdd[pd.Timestamp(p[1])]) - np.sqrt(tdd[pd.Timestamp(p[0])])
    out = {}
    for win in ("early", "late"):
        f = DIRS["ts"] / f"nisar_{win}_vertical.tif"
        if not f.exists() or win not in meta:
            continue
        h = local(read(f)[0])
        dx = np.mean([dsq(meta[win][k].split("-")) for k in ("asc", "desc")])
        out[win] = {**contrast(h, lc), "null": mask_null(h, lc == WET, shifts),
                    "expected_wet_minus_moss_mm": round(float(b_stefan * dx), 2)}
    return out


def all_summers(lc, min_dates=4):
    """Every summer stack in data/ts (run 09 --conncomp reference for all
    years): replication r of the local maps for every pair of summers (with
    the shift null), wetland contrast and Stefan slope per summer."""
    stacks = {}
    for f in sorted(DIRS["ts"].glob("s1_*_vert.tif")):
        with rasterio.open(f) as s:
            dates = list(s.descriptions)
            if len(dates) >= min_dates:
                stacks[int(f.stem.split("_")[1])] = (dates, s.read().astype("float32"))
    locs = {y: local(seasonal(d, st)) for y, (d, st) in stacks.items()}
    years = sorted(locs)
    pairs = {}
    for i, a in enumerate(years):
        for b in years[i + 1:]:
            r, _ = loc.corr(locs[a], locs[b])
            nm, ns, _ = loc.shift_null(locs[a], locs[b])
            pairs[f"{a}-{b}"] = {"r": round(r, 3), "null_mean": round(nm, 3),
                                 "null_std": round(ns, 3), "z": round((r - nm) / ns, 1)}
    tdd = an.tdd_series("norilsk_centre")
    per = {}
    for y, (d, st) in stacks.items():
        s = stefan(d, st, lc)
        per[str(y)] = {**contrast(locs[y], lc), "dates": [d[0], d[-1]], "n_dates": len(d),
                       "stefan_b": s["b"], "stefan_se": s["se_b"], "stefan_r2": s["r2"],
                       "tdd_at_last_date": round(float(tdd[pd.Timestamp(d[-1])]), 0),
                       "tdd_year": round(float(tdd[f"{y}-01-01":f"{y}-12-31"].max()), 0)}
    bs = np.array([per[str(y)]["stefan_b"] for y in years])
    ses = np.array([per[str(y)]["stefan_se"] for y in years])
    w = 1 / ses ** 2
    b_mean = float(np.sum(w * bs) / w.sum())
    chi2 = float(np.sum(w * (bs - b_mean) ** 2))
    warm = np.array([per[str(y)]["tdd_year"] for y in years])
    return {"years": years, "pairs": pairs, "per_summer": per,
            "stefan_weighted_mean": round(b_mean, 4), "stefan_weighted_se": round(float(w.sum() ** -0.5), 4),
            "stefan_chi2": round(chi2, 2), "stefan_dof": len(years) - 1,
            "stefan_p_constant": round(float(stats.chi2.sf(chi2, len(years) - 1)), 3),
            "r_stefan_b_vs_summer_tdd": round(float(np.corrcoef(bs, warm)[0, 1]), 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shifts", type=int, default=400)
    ap.add_argument("--all-summers-only", action="store_true",
                    help="only the all-summers replication (outputs/robustness_all_summers.json)")
    a = ap.parse_args()
    lc = read(DIRS["landcover"] / "worldcover_utm30.tif")[0]
    if a.all_summers_only:
        res = all_summers(lc)
        save_json(OUTPUTS / "robustness_all_summers.json", res)
        print(json.dumps({k: v for k, v in res.items() if k != "per_summer"}, indent=1))
        for y, v in res["per_summer"].items():
            print(y, {k: v[k] for k in ("wetland_minus_moss", "stefan_b", "stefan_se", "stefan_r2",
                                        "n_dates", "tdd_year")})
        return
    ctx = ts.Context()
    dem = ctx.dem.reshape(ctx.shape)
    lowland = (dem > 0) & (dem < LOWLAND_MAX_M)
    z = ctx.z
    stacks = build_variants(ctx)
    rep = {"variants": {}}

    # ---- variants
    seas_v1 = {}
    for v, per_year in stacks.items():
        r = {}
        loc_maps = {}
        for year, (dates, vert) in per_year.items():
            s = seasonal(dates, vert)
            loc_maps[year] = local(s)
            if v == "v1":
                seas_v1[year] = s
            r[str(year)] = {"dates": dates, **contrast(loc_maps[year], lc),
                            "industrial_minus_controls": zone_pair(dates, vert, z == 2, (z == 5) | (z == 6)),
                            "controlNW_minus_SE": zone_pair(dates, vert, z == 5, z == 6)}
        r["r_local"] = round(loc.corr(loc_maps[YEARS[0]], loc_maps[YEARS[1]])[0], 3)
        rep["variants"][v] = r
        print(f"{v}: r_local {r['r_local']}, wetland-moss "
              f"{[r[str(y)]['wetland_minus_moss'] for y in YEARS]} mm, ind-ctrl b "
              f"{[r[str(y)]['industrial_minus_controls']['diff_b'] for y in YEARS]}, NW-SE end "
              f"{[r[str(y)]['controlNW_minus_SE']['end_diff_mm'] for y in YEARS]} mm")

    # ---- null tests on v1
    shifts = random_shifts(a.shifts)
    v1 = stacks["v1"]
    locs = {y: local(seas_v1[y]) for y in YEARS}
    nulls = {"wetland": {}, "zones": {}, "hp_scales": {}}
    for y in YEARS:
        nulls["wetland"][str(y)] = mask_null(locs[y], lc == WET, shifts)
        nulls["zones"][str(y)] = zone_null(*v1[y], z, lowland, shifts[:200])
    nulls["wetland"]["mean_of_summers"] = mask_null((locs[YEARS[0]] + locs[YEARS[1]]) / 2, lc == WET, shifts)
    for hp in (1.0, 2.0, 4.0):
        la, lb = local(seas_v1[YEARS[0]], hp), local(seas_v1[YEARS[1]], hp)
        r, n = loc.corr(la, lb)
        null = loc.shift_null(la, lb)
        nulls["hp_scales"][f"{hp:g}km"] = {"r": round(r, 3), "null_mean": round(null[0], 3),
                                          "null_std": round(null[1], 3)}
    sets = facility_sets(ctx.prof)
    seas_list = [seas_v1[y] for y in YEARS]
    obs_n, obs_d, n_fac = facility_stats(seas_list, sets)
    null_n, null_d = [], []
    for dy, dx in shifts[:200]:
        k, d, _ = facility_stats(seas_list, sets, dy, dx)
        null_n.append(k)
        null_d.append(d)
    null_n, null_d = np.array(null_n), np.array(null_d)
    nulls["facilities"] = {"n_sites": n_fac, "persistent_obs": obs_n,
                           "persistent_null_mean": round(float(null_n.mean()), 1),
                           "persistent_null_p95": float(np.percentile(null_n, 95)),
                           "p_persistent": float(np.mean(null_n >= obs_n)),
                           "mean_inside_minus_ring_obs": round(obs_d, 2),
                           "mean_inside_minus_ring_null_std": round(float(np.nanstd(null_d)), 2),
                           "p_mean_two_sided": float(np.mean(np.abs(null_d - np.nanmean(null_d))
                                                             >= abs(obs_d - np.nanmean(null_d))))}
    rep["nulls_v1"] = nulls
    rep["stefan_v1"] = {str(y): stefan(*v1[y], lc) for y in YEARS}
    s1, s2 = (rep["stefan_v1"][str(y)] for y in YEARS)
    rep["stefan_v1"]["b_difference_z"] = round((s1["b"] - s2["b"]) / np.hypot(s1["se_b"], s2["se_b"]), 2)
    rep["nisar_wetland"] = nisar_wetland(lc, (s1["b"] + s2["b"]) / 2, shifts)
    save_json(OUTPUTS / "robustness.json", rep)

    print("\nnull tests (v1):")
    print(json.dumps(nulls, indent=1))
    print("NISAR wetland:", json.dumps(rep["nisar_wetland"], indent=1))
    print("stefan:", {k: (v["b"], v["se_b"], v["r2"]) if isinstance(v, dict) else v
                      for k, v in rep["stefan_v1"].items()})

    # ---- figure: Stefan series
    fig, ax = plt.subplots(figsize=(6, 4))
    for y in YEARS:
        s = rep["stefan_v1"][str(y)]
        ax.plot(s["sqrt_tdd"], s["wet_minus_moss_mm"], "o-", label=f"{y} (b={s['b']:.2f}±{s['se_b']:.2f})")
    ax.set_xlabel("sqrt(cumulative TDD) [sqrt(degC day)]")
    ax.set_ylabel("wetland minus moss/lichen, local [mm]")
    ax.axhline(0, color="grey", lw=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_robustness.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
