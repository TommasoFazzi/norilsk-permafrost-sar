"""Step 13: build the data files of the website (docs/, served by GitHub Pages).

Reads the pipeline outputs (data/ts, outputs/, outputs/web, data/climate) and
writes docs/data/:
  overlay_<name>.webp  displacement layers reprojected to Web Mercator,
                       diverging red (subsidence) / blue (uplift), alpha grows
                       with |value| so the basemap stays readable near zero
  zones.geojson        analysis zones
  site.json            everything the page shows: fits, comparisons, zone
                       time series, climate, NISAR stats, data timeline,
                       reference stability, verdict, overlay bounds
The page (docs/index.html) contains no numbers of its own: rerun this step
after ./run_v0.sh and the site updates.

    python 13_build_site.py
"""
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from PIL import Image
from rasterio.warp import Resampling, calculate_default_transform, reproject, transform_bounds

from common import ZONE_NAMES
from config import DIRS, OUTPUTS, POINTS, ZONE_CLIMATE

DOCS = Path(__file__).parent / "docs"
D = DOCS / "data"
MAX_PX = 1400
RED, MID, BLUE = np.array([230, 103, 103]), np.array([56, 56, 53]), np.array([57, 135, 229])


def colorize(a, lim):
    """Diverging RGBA: red < 0 (subsidence), blue > 0; alpha ~ |a|/lim."""
    t = np.clip(a / lim, -1, 1)
    rgb = np.where(t[..., None] < 0, MID + (RED - MID) * (-t[..., None]),
                   MID + (BLUE - MID) * t[..., None])
    alpha = np.clip(np.abs(t) * 1.4, 0.0, 1.0) * 235
    alpha[~np.isfinite(a)] = 0
    rgb[~np.isfinite(a)] = 0
    return np.dstack([rgb, alpha]).astype("uint8")


def overlay(src, name, lim, band=None, out_dir=None, prefix="data/"):
    """Reproject one band to EPSG:3857, colorize, save WebP in out_dir
    (default docs/data); return its path (prefix + name) and lon/lat corners."""
    with rasterio.open(src) as s:
        b = band or s.count
        tr, w, h = calculate_default_transform(s.crs, "EPSG:3857", s.width, s.height, *s.bounds)
        out = np.full((h, w), np.nan, "float32")
        reproject(rasterio.band(s, b), out, dst_transform=tr, dst_crs="EPSG:3857",
                  src_nodata=s.nodata if s.nodata is not None else np.nan,
                  dst_nodata=np.nan, resampling=Resampling.bilinear)
        bounds = transform_bounds("EPSG:3857", "EPSG:4326",
                                  tr.c, tr.f + tr.e * h, tr.c + tr.a * w, tr.f)
    img = Image.fromarray(colorize(out, lim))
    if img.width > MAX_PX:   # web weight: ~1400 px is plenty for a 50 km AOI
        img = img.resize((MAX_PX, round(img.height * MAX_PX / img.width)), Image.LANCZOS)
    img.save((out_dir or D) / f"overlay_{name}.webp", quality=82, method=6)
    w_, s_, e_, n_ = bounds
    return {"file": f"{prefix}overlay_{name}.webp", "lim_mm": lim,
            "coords": [[w_, n_], [e_, n_], [e_, s_], [w_, s_]]}


def climate_block():
    out = {}
    for name, pt in ZONE_CLIMATE.items():
        f = DIRS["climate"] / f"era5land_daily_{pt}.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, index_col="date", parse_dates=["date"])
        out[pt] = {}
        for y in (2020, 2025, 2026):
            s = df.loc[f"{y}-05-01":f"{y}-10-15", "tdd_cum"].dropna()
            out[pt][str(y)] = [[d.strftime("%Y-%m-%d"), round(float(v), 1)] for d, v in s.iloc[::3].items()]
    ann = pd.read_csv(DIRS["climate"] / "era5land_annual.csv")
    ann = ann[ann.point == "norilsk_centre"][["date", "maat_c", "tdd", "complete_year"]]
    out["annual_norilsk"] = ann.rename(columns={"date": "year"}).to_dict("records")
    return out


FACILITY_EN = {
    "Надеждинский металлургический завод": "Nadezhda Metallurgical Plant",
    "Норильская ТЭЦ-3": "TPP-3 power plant (2020 spill)",
    "Норильская ТЭЦ-2": "TPP-2 power plant",
    "Кайерканский угольный разрез": "Kayerkan open-pit coal mine",
    "Норильскгеология": "Norilskgeologiya base",
    "ВС-1, 2, 4": "Ventilation shafts VS-1, 2, 4",
    "ВС-5, ВС-6": "Ventilation shafts VS-5, 6",
    "Вентиляционный ствол рудника Маяк": "Mayak mine ventilation shaft",
    "ЮЗВС": "South-west ventilation shaft",
    "Талнахская обогатительная фабрика": "Talnakh concentrator",
    "Медный завод": "Copper plant", "Никелевый завод": "Nickel plant",
    "Норильская нефтебаза": "Norilsk oil depot",
}
LU_EN = {"industrial": "industrial site", "landfill": "landfill / tailings", "quarry": "quarry"}


def local_block():
    """Step 14 outputs: reproducibility, land cover, facilities, hotspots."""
    f = OUTPUTS / "local_meta.json"
    if not f.exists():
        return None
    meta = json.loads(f.read_text())
    lc = pd.read_csv(OUTPUTS / "local_landcover.csv")
    fac = pd.read_csv(OUTPUTS / "local_facilities.csv")
    fac["label"] = [FACILITY_EN.get(n, n) if isinstance(n, str) else f"Unnamed {LU_EN[lu]}"
                    for n, lu in zip(fac["name"], fac["landuse"])]
    hot = pd.read_csv(OUTPUTS / "local_hotspots.csv")
    clean = lambda df: df.round(3).replace({np.nan: None}).to_dict("records")
    by_lu = fac.groupby("landuse")[["d_2020", "d_2025"]].mean().round(2)
    return {**meta, "landcover": clean(lc), "facilities": clean(fac[fac.persistent]),
            "facilities_all_n": int(len(fac)),
            "facility_r": round(float(fac[["d_2020", "d_2025"]].corr().iloc[0, 1]), 3),
            "facilities_by_landuse": by_lu.reset_index().to_dict("records"),
            "hotspots": clean(hot.head(25))}


ABS_SE_USABLE = 0.25   # 2*se_b below 0.5 mm/sqrt(degC d): absolute slope is informative


def verdict(fits, comps):
    """Rule-based reading (documented on the page).

    B anthropogenic: PAIRED slope of (industrial - mean of controls) below zero
      by > 2 sigma (common-mode atmosphere cancels in the difference).
    A climatic: both controls' absolute slopes below zero by > 2 sigma, only
      if the absolute slope is informative (se_b <= ABS_SE_USABLE); else
      'untestable' (absolute signal dominated by troposphere).
    C stable: controls and industry absolute slopes within 2 sigma of zero AND
      informative; else 'untestable'."""
    if fits.empty:
        return {"status": "pending"}
    res = {"status": "ready", "summers": {}}
    for yr, g in fits.groupby("summer"):
        g = g.set_index("zone")
        c = comps[comps.summer == yr].set_index(["zone_a", "zone_b"])
        get = lambda zone: (g.at[zone, "b"], g.at[zone, "se_b"]) if zone in g.index else (np.nan, np.nan)
        pz = lambda a, b: float(c.loc[(a, b), "z"]) if (a, b) in c.index else None
        ctrl = [get("control_NW"), get("control_SE")]
        ind = get("industrial")
        usable = all(np.isfinite(se) and se <= ABS_SE_USABLE for _, se in ctrl + [ind])
        z_ind = pz("industrial", "controls")
        nw_se = c.loc[("control_NW", "control_SE"), "seasonal_diff_mm"] if ("control_NW", "control_SE") in c.index else np.nan
        ind_ctl = c.loc[("industrial", "controls"), "seasonal_diff_mm"] if ("industrial", "controls") in c.index else np.nan
        climatic = ("untestable" if not usable else
                    "yes" if all(b / se < -2 for b, se in ctrl) else "no")
        stable = ("untestable" if not usable else
                  "yes" if all(abs(b / se) < 2 for b, se in ctrl + [ind]) else "no")
        anthro = "yes" if (z_ind is not None and z_ind < -2) else "no"
        res["summers"][str(yr)] = {
            "climatic": climatic, "anthropogenic": anthro, "structural": stable,
            "z_industrial_vs_controls": z_ind,
            "seasonal_industrial_minus_controls_mm": None if not np.isfinite(ind_ctl) else round(float(ind_ctl), 1),
            "seasonal_controlNW_minus_SE_mm": None if not np.isfinite(nw_se) else round(float(nw_se), 1),
            "controls_disagree_as_much": bool(np.isfinite(nw_se) and np.isfinite(ind_ctl) and abs(nw_se) >= abs(ind_ctl)),
            "abs_slope_informative": usable,
            "max_abs_se_b": round(float(max(se for _, se in ctrl + [ind])), 3)}
    for h in ("climatic", "anthropogenic", "structural"):
        vals = [v[h] for v in res["summers"].values()]
        res[h] = {"yes": vals.count("yes"), "no": vals.count("no"), "untestable": vals.count("untestable")}
    res["n_summers"] = len(res["summers"])
    return res


def main():
    D.mkdir(parents=True, exist_ok=True)
    site = {"generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "points": POINTS, "zones": ZONE_NAMES, "overlays": {}}

    # zones
    zg = OUTPUTS / "web" / "zones.geojson"
    if zg.exists():
        shutil.copy(zg, D / "zones.geojson")

    # S1 seasonal overlays (common colour limit so summers are comparable)
    s1 = sorted(DIRS["ts"].glob("s1_*_vert.tif"))
    if s1:
        lasts = [rasterio.open(f).read(rasterio.open(f).count) for f in s1]
        lim = float(np.nanpercentile(np.abs(np.concatenate([a[np.isfinite(a)] for a in lasts])), 97))
        lim = max(5.0, round(lim / 5) * 5)
        for f in s1:
            yr = f.stem.split("_")[1]
            with rasterio.open(f) as s:
                d0, d1 = s.descriptions[0], s.descriptions[-1]
            site["overlays"][f"s1_{yr}"] = {**overlay(f, f"s1_{yr}", lim),
                                            "label": f"Sentinel-1 summer {yr}",
                                            "detail": f"vertical, {d0[:4]}-{d0[4:6]}-{d0[6:]} to {d1[:4]}-{d1[4:6]}-{d1[6:]}"}
            meta = json.loads((DIRS["ts"] / f"s1_{yr}_meta.json").read_text())
            site.setdefault("s1_meta", {})[yr] = meta
    # NISAR overlays
    for win in ("early", "late"):
        for comp in ("vertical", "east"):
            f = DIRS["ts"] / f"nisar_{win}_{comp}.tif"
            if f.exists():
                site["overlays"][f"nisar_{win}_{comp}"] = {
                    **overlay(f, f"nisar_{win}_{comp}", 20.0),
                    "label": f"NISAR {win} summer 2026, {comp}",
                    "detail": "12-day L-band pair, asc 111 + desc 6"}
    for name, lim, label, detail in [
            ("local_persistent", 15.0, "Repeats in both summers",
             "local seasonal motion vs ~2 km surroundings, mean of 2020 and 2025 where both agree in sign"),
            ("local_2020", 15.0, "Local seasonal 2020", "vs ~2 km surroundings, sqrt(TDD) fit"),
            ("local_2025", 15.0, "Local seasonal 2025", "vs ~2 km surroundings, sqrt(TDD) fit")]:
        f = DIRS["ts"] / f"{name}.tif"
        if f.exists():
            site["overlays"][name] = {**overlay(f, name, lim, band=1), "label": label, "detail": detail}
    site["local"] = local_block()
    nm = DIRS["ts"] / "nisar_meta.json"
    site["nisar"] = json.loads(nm.read_text()) if nm.exists() else None

    # analysis tables
    ts_csv = OUTPUTS / "web" / "zone_timeseries.csv"
    fits_csv = OUTPUTS / "zone_fits.csv"
    comp_csv = OUTPUTS / "zone_comparisons.csv"
    gap_csv = OUTPUTS / "gap_prediction.csv"
    ts = pd.read_csv(ts_csv) if ts_csv.exists() else pd.DataFrame()
    fits = pd.read_csv(fits_csv) if fits_csv.exists() else pd.DataFrame()
    comps = pd.read_csv(comp_csv) if comp_csv.exists() else pd.DataFrame()
    site["timeseries"] = ts.round(3).to_dict("records")
    site["fits"] = fits.round(4).replace({np.nan: None}).to_dict("records")
    site["comparisons"] = comps.round(3).replace({np.nan: None}).to_dict("records")
    site["gap"] = (pd.read_csv(gap_csv).round(3).replace({np.nan: None}).to_dict("records")
                   if gap_csv.exists() else [])
    site["verdict"] = verdict(fits, comps)
    L = site["local"]
    if L and site["verdict"].get("status") == "ready":
        wet = next((r for r in L["landcover"] if r["code"] == 90), None)
        sig = lambda r, y: r[f"local_{y}_mm"] < -2 * r[f"se_{y}"]
        site["verdict"]["climatic_local"] = bool(
            wet and L["r_local"] > L["r_null_shift_mean"] + 5 * L["r_null_shift_std"]
            and all(sig(wet, y) for y in (2020, 2025)))
    site["climate"] = climate_block()
    site["zone_climate"] = ZONE_CLIMATE

    # data availability timeline
    sc = pd.read_csv(DIRS["s1"] / "s1_bursts.csv")
    per_year = sc.drop_duplicates("date").assign(y=lambda d: d.date.str[:4]).groupby("y").size()
    site["timeline"] = {str(y): {"s1_dates": int(per_year.get(str(y), 0)),
                                 "processed": y in (2020, 2025), "nisar": y == 2026}
                        for y in range(2017, 2027)}

    # images
    for src in sorted((OUTPUTS / "web").glob("s2_spill_*.png")):
        img = Image.open(src).convert("RGB")
        img.thumbnail((1600, 3200))
        img.save(D / src.with_suffix(".jpg").name, quality=82, optimize=True)

    (D / "site.json").write_text(json.dumps(site, default=str))
    print(f"site.json: {len(site['overlays'])} overlays, {len(site['fits'])} fits, "
          f"verdict {site['verdict'].get('status')}")


if __name__ == "__main__":
    main()
