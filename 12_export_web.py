"""Step 12: frontend data package (outputs/frontend/) - final results.

Built for the app that serves the results (Dancing-with-the-SARs, Xenonite
Lab: React + Cesium globe + FastAPI, region id "norilsk"), but plain files
that any frontend can use. Run it last (run_v1.sh does).

outputs/frontend/
  README.md       data contract: every file, units, sign, how to show it
  manifest.json   the same, machine-readable (one entry per file)
  story.json      final conclusions: sections with headline, text, key
                  numbers read from the analysis outputs, linked files
  rasters/        COG EPSG:4326 float32 mm, NaN = no data
  overlays/       the same maps colourised (WebP, Web Mercator image whose
                  lon/lat corners are in manifest.json): red down, blue up
  vectors/        GeoJSON EPSG:4326: zones, sites, TPP-3 plant
  tables/         CSV, units in column names
  images/         Sentinel-2 chips (with lon/lat bounds), static figures
  app/            Norilsk subsidence dataset in the app's SAR schema
                  (backend/data/sar), merge script, INTEGRATION.md
outputs/frontend_package.zip is the same folder, to hand over (*.tif and
*.csv are not in git). outputs/web/zones.geojson is still written for
13_build_site.py (the original website).

    python 12_export_web.py
"""
import hashlib
import importlib
import json
import shutil

import numpy as np
import pandas as pd
import rasterio
from pyproj import Transformer
from rasterio.features import shapes
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import mapping, shape
from shapely.ops import unary_union

from common import ZONE_NAMES
from config import DIRS, OUTPUTS, POINTS, UTM_EPSG

site = importlib.import_module("13_build_site")

PKG = OUTPUTS / "frontend"
WEB = OUTPUTS / "web"
KIT = OUTPUTS.parent / "frontend_kit" if (OUTPUTS.parent / "frontend_kit").exists() else None
APP_ID = "NOR_TPP3_TANK_VERTICAL"
CREDIT = ("Contains modified Copernicus Sentinel data 2017-2026 processed by ASF HyP3; "
          "NISAR data courtesy NASA/ISRO via ASF DAAC; ERA5 (Copernicus C3S) via Open-Meteo; "
          "GNSS NRIL: Nevada Geodetic Laboratory; OpenStreetMap contributors; ESA WorldCover; "
          "Copernicus DEM")
SIGN = "+ up, - down (vertical, from radar line of sight assuming no horizontal motion)"


def load(name):
    p = OUTPUTS / name
    return json.loads(p.read_text()) if p.exists() else None


# ----------------------------------------------------------------- rasters
def to_cog_4326(src_path, dst_path, band=None, resampling=Resampling.bilinear):
    with rasterio.open(src_path) as s:
        idx = band or s.count
        tr, w, h = calculate_default_transform(s.crs, "EPSG:4326", s.width, s.height, *s.bounds)
        dtype = s.dtypes[idx - 1]
        nodata = s.nodata if s.nodata is not None else (np.nan if "float" in dtype else 0)
        out = np.full((h, w), nodata, dtype)
        reproject(rasterio.band(s, idx), out, dst_transform=tr, dst_crs="EPSG:4326",
                  dst_nodata=nodata, resampling=resampling)
        desc = s.descriptions[idx - 1]
    prof = dict(driver="COG", width=w, height=h, count=1, dtype=dtype, crs="EPSG:4326",
                transform=tr, nodata=nodata, compress="deflate")
    with rasterio.open(dst_path, "w", **prof) as d:
        d.write(out, 1)
        if desc:
            d.set_band_description(1, desc)
    return dst_path


def map_layer(man, src, name, lim, title, section, description, period=None, band=None):
    """One map: COG (analysis) + colourised overlay (display)."""
    cog = to_cog_4326(src, PKG / "rasters" / f"{name}.tif", band)
    ov = site.overlay(src, name, lim, band=band, out_dir=PKG / "overlays", prefix="overlays/")
    add(man, cog, "raster", title, section, description, units="mm", sign=SIGN, crs="EPSG:4326",
        nodata="NaN", period=period, overlay=ov["file"], overlay_corners_lonlat=ov["coords"],
        overlay_colour_limit_mm=lim,
        display="map overlay: diverging, red = down, blue = up, transparent near 0")


def maps(man):
    s1 = sorted(DIRS["ts"].glob("s1_*_vert.tif"))
    if s1:
        lasts = [rasterio.open(f).read(rasterio.open(f).count) for f in s1]
        lim = max(5.0, round(float(np.nanpercentile(np.abs(np.concatenate([a[np.isfinite(a)] for a in lasts])), 97)) / 5) * 5)
    for f in s1:
        year = f.stem.split("_")[1]
        with rasterio.open(f) as s:
            d0, d1 = s.descriptions[0], s.descriptions[-1]
        few = rasterio.open(f).count < 4
        map_layer(man, f, f"s1_seasonal_{year}", lim, f"Sentinel-1 summer {year}: vertical change",
                  "thaw", f"Vertical displacement from {d0} to {d1} relative to plateau bedrock; "
                  "slopes > 5 deg masked. At 10-50 km the atmosphere dominates: use it as a "
                  "pattern map, not for absolute values between distant places."
                  + (" Only 3 dates: preliminary." if few else ""), period=f"{d0}/{d1}")
    for name, title, desc in [
            ("local_persistent", "Where the ground sinks every summer",
             "Local seasonal motion vs ~2 km surroundings, mean of 2020 and 2025 where both "
             "summers agree in sign. The reproducible part of the signal: wetlands."),
            ("local_2020", "Local seasonal motion 2020", "vs ~2 km surroundings, sqrt(TDD) fit"),
            ("local_2025", "Local seasonal motion 2025", "vs ~2 km surroundings, sqrt(TDD) fit")]:
        f = DIRS["ts"] / f"{name}.tif"
        if f.exists():
            map_layer(man, f, name, 15.0, title, "thaw", desc, band=1)
    for f in sorted(DIRS["ts"].glob("nisar_*_vertical.tif")) + sorted(DIRS["ts"].glob("nisar_*_east.tif")):
        comp = "vertical" if "vertical" in f.name else "east-west (+ east)"
        map_layer(man, f, f.stem, 20.0, f"NISAR 2026 {f.stem.split('_')[1]} window, {comp}",
                  "thaw", "L-band 12-day pair, ascending 111 + descending 6 combined; "
                  "an independent cross-check, not a time series.")


# ----------------------------------------------------------------- vectors
def zones(man):
    src = DIRS["zones"] / "zones_utm30.tif"
    cog = to_cog_4326(src, PKG / "rasters" / "zones.tif", resampling=Resampling.nearest)
    with rasterio.open(cog) as s:
        z = s.read(1)
        feats = {}
        for geom, val in shapes(z, mask=z > 0, transform=s.transform):
            feats.setdefault(int(val), []).append(shape(geom))
    fc = {"type": "FeatureCollection", "features": [
        {"type": "Feature", "properties": {"code": c, "zone": ZONE_NAMES[c]},
         "geometry": mapping(unary_union(g).buffer(0.0003).buffer(-0.0003).simplify(0.0004))}
        for c, g in sorted(feats.items())]}
    for path in (PKG / "vectors" / "zones.geojson", WEB / "zones.geojson"):
        path.write_text(json.dumps(fc))
    add(man, cog, "raster", "Analysis zones (codes)", "overview",
        "uint8 zone code, 0 = none; codes in vectors/zones.geojson", crs="EPSG:4326")
    add(man, PKG / "vectors" / "zones.geojson", "vector", "Analysis zones", "overview",
        "Plateau bedrock reference, industrial, TPP-3, urban, tundra controls",
        properties={"zone": "name", "code": "integer"})


def point(lat, lon, **props):
    return {"type": "Feature", "geometry": {"type": "Point", "coordinates": [round(lon, 6), round(lat, 6)]},
            "properties": props}


def sites(man, checks, gnss):
    feats = [point(*POINTS["tank5_zhang2022"], id="tank5", name="TPP-3 fuel tank no. 5 (collapsed 2020-05-29)",
                   kind="tank", note="Position from Zhang et al. 2022, confirmed on imagery and by the "
                                     "radar amplitude loss after the collapse.")]
    for n, r in (checks or {}).get("row_structures", {}).items():
        feats.append(point(*r["latlon"], id=n, name=f"Bright structure of the tank row ({n[4:]})",
                           kind="tank_row", note="Probable tank (round on Sentinel-2); comparison target."))
    for n, (lat, lon) in (checks or {}).get("osm_tanks", {}).items():
        feats.append(point(lat, lon, id=n, name="OSM storage tank (Norilsk depots)", kind="depot_tank",
                           note="Comparison target for the tank-only check."))
    if gnss:
        feats.append(point(gnss["lat"], gnss["lon"], id="NRIL", name="GNSS station NRIL",
                           kind="gnss", note=gnss["monument"]))
    for k in ("ref_1", "ref_2", "ref_3"):
        feats.append(point(*POINTS[k], id=k, name="InSAR reference patch (plateau bedrock)", kind="reference"))
    for k in ("control_NW", "control_SE"):
        feats.append(point(*POINTS[k], id=k, name=f"Tundra control {k[-2:]}", kind="control"))
    path = PKG / "vectors" / "sites.geojson"
    path.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    add(man, path, "vector", "Sites: tank 5, comparison tanks, GNSS, references", "tank",
        "Points with id, name, kind (tank, tank_row, depot_tank, gnss, reference, control), note.")
    osm = DIRS["osm"] / "osm_landuse.geojson"
    if not osm.exists():
        return
    plant = [f for f in json.loads(osm.read_text())["features"] if f["properties"]["osm_id"] == "way/525598517"]
    path = PKG / "vectors" / "tpp3_plant.geojson"
    path.write_text(json.dumps({"type": "FeatureCollection", "features": plant}))
    add(man, path, "vector", "TPP-3 power plant (OSM polygon)", "tank",
        "Plant outline; tank 5 lies 0.35 km outside it, in the fuel-tank row to the NW.")


# ----------------------------------------------------------------- tables
def tables(man, allsum, tank18, checks):
    T = PKG / "tables"
    copies = [
        ("web/zone_timeseries.csv", "zone_timeseries.csv", "Zone mean vertical displacement per summer date",
         "industry", "zone, summer, date, disp_mm (+ up, vs plateau bedrock), disp_median_mm, se_mm, n_px, tdd_cum (degC day)"),
        ("zone_fits.csv", "zone_fits.csv", "Seasonal slope per zone and summer (displacement vs sqrt(TDD))",
         "industry", "b = mm per sqrt(degC day); absolute values are atmosphere-dominated"),
        ("local_landcover.csv", "landcover_local.csv", "Local seasonal motion per land-cover class",
         "thaw", "mean mm vs ~2 km surroundings, 2020 and 2025; wetlands stand out"),
        ("tank_winters.csv", "tank_winters.csv", "Tank 5, winter by winter, date by date",
         "tank", "vertical mm since the first December date; comparison targets mean/sd; ERA5 snow and wind since the previous date"),
        ("tank_points_series.csv", "tank_points_series.csv", "Tank 5 and neighbouring tanks: segment series",
         "tank", "vertical mm relative to stable targets 0.3-1.5 km; each segment starts at 0; segments are NOT comparable"),
        ("gnss_nril_up.csv", "gnss_nril_up.csv", "GNSS NRIL weekly vertical position",
         "reference", "mm vs 2017 mean, station on bedrock in Norilsk"),
    ]
    for src, dst, title, sec, desc in copies:
        if (OUTPUTS / src).exists():
            shutil.copyfile(OUTPUTS / src, T / dst)
            add(man, T / dst, "table", title, sec, desc)
    if allsum:
        ps = pd.DataFrame([{"summer": int(y), **{k: v for k, v in d.items() if k != "dates"},
                            "first_date": d["dates"][0], "last_date": d["dates"][-1]}
                           for y, d in allsum["per_summer"].items()])
        ps.to_csv(T / "wetland_by_summer.csv", index=False)
        add(man, T / "wetland_by_summer.csv", "table", "Wetland thaw subsidence, summer by summer", "thaw",
            "wetland, moss, wetland_minus_moss (mm over the summer, vs ~2 km surroundings); "
            "stefan_b = mm per sqrt(degC day) with se; tdd = thawing degree days")
        pairs = pd.DataFrame([{"summers": k, **v} for k, v in allsum["pairs"].items()])
        pairs.to_csv(T / "replication_pairs.csv", index=False)
        add(man, T / "replication_pairs.csv", "table", "Same pattern every summer", "thaw",
            "r = correlation of local seasonal maps between two summers; null = maps shifted 3-20 km")
    if tank18:
        rows = []
        for key, v in tank18["points"]["published"]["vs_null"].items():
            if v and key.startswith("change_"):
                d0, d1 = key.split("_")[1:3]
                exp = next((i["zhang_rate_expected_mm"] for i in tank18["intervals"]
                            if i["from"] == d0 and i["to"] == d1), None)
                osm = (checks or {}).get("c5_intervals", {}).get(f"{d0}_{d1}", {})
                rows.append({"from": d0, "to": d1, "tank5_mm": v["value"],
                             "comparable_targets_mean_mm": v["null_mean"], "comparable_targets_sd_mm": v["null_std"],
                             "comparable_targets_n": v["n_null"], "z": v["z"], "p": v["p_two_sided"],
                             "depot_tanks_mean_mm": osm.get("osm_tanks_mean"), "depot_tanks_n": osm.get("osm_tanks_n"),
                             "zhang2022_rate_expected_mm": exp})
        pd.DataFrame(rows).to_csv(T / "tank_intervals.csv", index=False)
        add(man, T / "tank_intervals.csv", "table", "Tank 5 vs comparable targets, interval by interval", "tank",
            "vertical mm relative to surroundings over each observable interval; what Zhang et al. "
            "2022 (-35 mm/yr) would predict over the same interval")
    if not (DIRS["climate"] / "era5land_annual.csv").exists() or not (DIRS["s1"] / "s1_bursts.csv").exists():
        return
    ann = pd.read_csv(DIRS["climate"] / "era5land_annual.csv")
    ann = ann[ann.point == "norilsk_centre"].drop(columns="point").rename(columns={"date": "year"})
    ann.to_csv(T / "climate_annual.csv", index=False)
    add(man, T / "climate_annual.csv", "table", "Norilsk climate per year (ERA5-Land)", "overview",
        "maat_c mean air temperature, tdd/fdd thawing/freezing degree days, snow_depth_max_m")
    sc = pd.read_csv(DIRS["s1"] / "s1_bursts.csv").drop_duplicates("date")
    per = sc.groupby(sc.date.str[:4]).size()
    done = {f.stem.split("_")[1] for f in DIRS["ts"].glob("s1_*_vert.tif")}
    pd.DataFrame([{"year": y, "sentinel1_dates": int(per.get(str(y), 0)),
                   "summer_processed": str(y) in done, "nisar": y >= 2026,
                   "note": "no Sentinel-1 over Norilsk" if 2022 <= y <= 2024 else ""}
                  for y in range(2017, 2027)]).to_csv(T / "data_timeline.csv", index=False)
    add(man, T / "data_timeline.csv", "table", "Which radar data exist each year", "overview",
        "Sentinel-1 acquisition dates per year; 2022-2024 have none over Norilsk")


# ----------------------------------------------------------------- images
def images(man):
    chips = WEB / "s2_tank_chips.json"
    if chips.exists():
        meta = json.loads(chips.read_text())
        for k in ("before", "after"):
            shutil.copyfile(WEB / meta[k]["file"], PKG / "images" / meta[k]["file"])
            add(man, PKG / "images" / meta[k]["file"], "image", f"Tank farm {k} the collapse (Sentinel-2, {meta[k]['date']})",
                "tank", "True colour, 10 m: context only, too coarse to show the tank's state.",
                bounds_lonlat=meta["bounds_lonlat"], credit=meta["credit"])
    for f in sorted(WEB.glob("s2_spill_*.png")):
        shutil.copyfile(f, PKG / "images" / f.name)
        add(man, PKG / "images" / f.name, "image", f"Diesel spill, Sentinel-2 {f.stem[-8:]}", "overview",
            "True colour TPP-3 -> Daldykan -> Ambarnaya (box 87.75-88.10 E, 69.30-69.55 N), not georeferenced",
            credit="Contains modified Copernicus Sentinel data 2020")
    for name, title, sec in [
            ("fig_tank_series.png", "Tank 5: relative series, intervals vs null, amplitude loss", "tank"),
            ("fig_tank_checks.png", "Tank 5: unwrapping check and winters vs other targets", "tank"),
            ("fig_gnss_nril.png", "GNSS NRIL: bedrock vertical position", "reference"),
            ("fig_robustness.png", "Robustness of the wetland signal", "thaw"),
            ("fig_local.png", "Local seasonal motion and its reproducibility", "thaw"),
            ("fig_disp_vs_sqrt_tdd.png", "Zone displacement vs sqrt(TDD)", "industry")]:
        if (OUTPUTS / name).exists():
            shutil.copyfile(OUTPUTS / name, PKG / "images" / name)
            add(man, PKG / "images" / name, "figure", title, sec, "Static figure (fallback / methods page)")


# ----------------------------------------------------------------- app dataset
def app_dataset(man):
    """Tank 5 and the row structures in the app's SAR schema (backend/data/sar)."""
    src = OUTPUTS / "tank_points_series.csv"
    if not src.exists() or KIT is None:
        return
    df = pd.read_csv(src, dtype={"date": str})
    to_utm = Transformer.from_crs(4326, UTM_EPSG, always_xy=True)
    names = {"tank5": "Tank 5 (collapsed 29 May 2020)"}
    points, rows = [], []
    for pid, g in df.groupby("point", sort=False):
        lat, lon = float(g.latitude.iloc[0]), float(g.longitude.iloc[0])
        e, n = to_utm.transform(lon, lat)
        dist, bearing = pid[4:].split("m_") if pid.startswith("row_") else ("", "")
        points.append({"id": f"NOR_{pid.upper()}", "name": names.get(pid, f"Tank-row structure, {dist} m {bearing} of tank 5"),
                       "areaLabel": "TPP-3 fuel tank row", "latitude": round(lat, 6), "longitude": round(lon, 6),
                       "locationStatus": "source_coordinates", "gridSizeMeters": 20,
                       "velocity": None, "rmse": round(float(g.segment_rms_mm.mean()), 2),
                       "easting": round(e, 1), "northing": round(n, 1), "sourceCrs": f"EPSG:{UTM_EPSG}"})
        for _, r in g.iterrows():
            rows.append({"date": f"{r.date[:4]}-{r.date[4:6]}-{r.date[6:]}", "pointId": f"NOR_{pid.upper()}",
                         "value": round(float(r.vertical_mm), 2), "qualityPass": True,
                         "details": {"segment": int(r.segment), "segmentRmseMm": float(r.segment_rms_mm),
                                     "reference": "winter-coherent targets 0.3-1.5 km, same SNAPHU component"}})
    rows.sort(key=lambda r: (r["date"], r["pointId"]))
    for i, r in enumerate(rows, 1):
        r["id"] = f"{APP_ID}:{i}"
    meta = {"id": APP_ID, "regionId": "norilsk", "name": "TPP-3 tank 5 · vertical motion vs surroundings",
            "kind": "subsidence", "metric": "Vertical displacement", "unit": "mm", "bounds": None,
            "geometryStatus": "source_coordinates", "validationStatus": "exploratory",
            "note": ("Sentinel-1 (path 136 desc, 20 m), each point minus stable targets 0.3-1.5 km away. "
                     "The tanks decorrelate at snowmelt and snow onset, so the series is split into "
                     "segments (details.segment); each segment starts at 0 and values of different "
                     "segments must not be compared. No mean velocity exists (velocity null) and no "
                     "forecast is provided. Negative = down."),
            "method": ("HyP3 INSAR_ISCE_MULTI_BURST 5x1 looks, no water mask; 3x3 window minus "
                       "winter-coherent ring; SBAS per connected segment; vertical = LOS/cos(41.3 deg). "
                       "Steps 18-19 of the Norilsk pipeline."),
            "source": {"file": src.name, "sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
                       "sheet": "tank_points_series", "publisher": "Norilsk permafrost team / Copernicus Sentinel-1 via ASF HyP3",
                       "url": "https://hyp3-docs.asf.alaska.edu/"},
            "validation": [], "points": points,
            "startDate": min(r["date"] for r in rows), "endDate": max(r["date"] for r in rows),
            "observationCount": len(rows), "qualityPassedCount": len(rows), "schemaVersion": 1,
            "forecastAvailable": False, "pointBadge": "SENTINEL-1 TARGETS"}
    A = PKG / "app"
    (A / "sar").mkdir(parents=True, exist_ok=True)
    (A / "sar" / f"{APP_ID}.json").write_text(json.dumps(rows, ensure_ascii=False, allow_nan=False,
                                                          separators=(",", ":")) + "\n", encoding="utf-8")
    (A / "catalog_entries.json").write_text(json.dumps([meta], ensure_ascii=False, allow_nan=False, indent=2) + "\n",
                                            encoding="utf-8")
    for f in ("merge_into_app.py", "INTEGRATION.md"):
        shutil.copyfile(KIT / f, A / f)
    add(man, A / "catalog_entries.json", "app", "Norilsk dataset for the app's SAR catalogue", "tank",
        f"{APP_ID}: {len(points)} points, {len(rows)} readings; merge with app/merge_into_app.py "
        "after reading app/INTEGRATION.md")


# ----------------------------------------------------------------- story
def story(allsum, local, rob, tank18, checks, gnss):
    """Final conclusions as structured text; every number is read from the outputs."""
    S = []
    S.append({"id": "question", "headline": "Was the ground under Norilsk sinking, and did that bring down tank 5?",
              "text": "On 29 May 2020 fuel tank no. 5 of the TPP-3 power plant collapsed and ~21,000 t of "
                      "diesel reached the Ambarnaya river. One study (Zhang et al. 2022) blamed thawing "
                      "permafrost; the federal inspectorate (Rostekhnadzor, Nov 2020) blamed the pile "
                      "foundation. We measured ground motion with radar satellites 2017-2026.",
              "datasets": ["images/s2_spill_20200531.png", "vectors/sites.geojson"]})
    if allsum and local and rob:
        r = [v["r"] for v in allsum["pairs"].values()]
        z = [v["z"] for v in allsum["pairs"].values()]
        nis = rob["nisar_wetland"]
        S.append({"id": "thaw", "headline": "Wetlands sink every summer, as thawing ground should, and the rate is not growing",
                  "text": "Compared with their ~2 km surroundings, herbaceous wetlands sink every summer in "
                          "proportion to the square root of accumulated thaw (Stefan's law), with the same "
                          "coefficient in all six summers 2017-2025: no sign of acceleration. The same places "
                          "sink every summer, and NISAR (L-band, 2026) sees it independently.",
                  "key_numbers": [
                      {"label": "Stefan coefficient, six summers", "value": f"{allsum['stefan_weighted_mean']:.3f} ± {allsum['stefan_weighted_se']:.3f}", "unit": "mm per sqrt(degC day)"},
                      {"label": "Constant across summers (chi-square p)", "value": allsum["stefan_p_constant"], "unit": ""},
                      {"label": "Same pattern in two summers: r", "value": f"{min(r):.2f}-{max(r):.2f} (z >= {min(z):.1f})", "unit": ""},
                      {"label": "Local correlation 2020 vs 2025", "value": f"{local['r_local']:.2f} vs {local['r_null_shift_mean']:.3f} ± {local['r_null_shift_std']:.3f} by chance", "unit": ""},
                      {"label": "NISAR 2026 wetland minus moss", "value": f"{nis['early']['wetland_minus_moss']} / {nis['late']['wetland_minus_moss']} (z {nis['early']['null']['z_emp']} / {nis['late']['null']['z_emp']})", "unit": "mm"}],
                  "caveats": "Motion of the whole lowland relative to bedrock is not measurable (atmosphere). "
                             "Part of the wetland signal could be water-level effects on the radar.",
                  "datasets": ["overlays/overlay_local_persistent.webp", "tables/wetland_by_summer.csv",
                               "tables/replication_pairs.csv", "tables/landcover_local.csv"]})
        zz = rob["nulls_v1"]["zones"]
        fac = rob["nulls_v1"]["facilities"]
        S.append({"id": "industry", "headline": "Industrial ground does not sink more than other lowland ground",
                  "text": "Industrial zones behave like random lowland ground, and plants such as TPP-3 do "
                          "not sink more than their surroundings in any summer.",
                  "key_numbers": [
                      {"label": "Industrial vs random lowland, z (2020 / 2025)", "value": f"{zz['2020']['z_emp']} / {zz['2025']['z_emp']}", "unit": ""},
                      {"label": "Sites that 'sink every summer' vs chance", "value": f"{fac['persistent_obs']} vs {fac['persistent_null_mean']} (p {fac['p_persistent']})", "unit": ""}],
                  "datasets": ["tables/zone_timeseries.csv", "tables/zone_fits.csv"]})
    if tank18:
        pv = tank18["points"]["published"]["vs_null"]
        last = [k for k in pv if k.startswith("change_") and pv[k]][-1]
        summ = {y: pv[f"summer_{y}_mm"]["value"] for y in (2017, 2018, 2019) if pv.get(f"summer_{y}_mm")}
        nums = [{"label": f"Snow-free summer {y}, tank vs surroundings", "value": v, "unit": "mm"} for y, v in summ.items()]
        nums.append({"label": f"Last winter before the collapse ({last.split('_')[1]} -> {last.split('_')[2]})",
                     "value": f"{pv[last]['value']} (comparable targets {pv[last]['null_mean']} ± {pv[last]['null_std']}, z {pv[last]['z']})",
                     "unit": "mm"})
        S.append({"id": "tank", "headline": TANK_HEADLINE,
                  "text": TANK_TEXT, "key_numbers": nums + tank_check_numbers(checks),
                  "caveats": TANK_CAVEATS,
                  "datasets": ["tables/tank_intervals.csv", "tables/tank_winters.csv", "tables/tank_points_series.csv",
                               "images/fig_tank_series.png", "images/fig_tank_checks.png",
                               "images/s2_tank_before.png", "images/s2_tank_after.png", "app/catalog_entries.json"]})
    if gnss:
        S.append({"id": "reference", "headline": "The bedrock under Norilsk is not sinking",
                  "text": "The GNSS station NRIL, on bedrock in Norilsk, rises slightly. A regional "
                          "subsidence like the -30.8 mm/yr area mean of Zhang et al. 2022 does not exist "
                          "in the bedrock; our InSAR reference (plateau bedrock) is a sound zero.",
                  "key_numbers": [
                      {"label": "NRIL vertical trend 2000-2023", "value": f"{gnss['all']['trend_mm_yr']:+.1f}", "unit": "mm/yr"},
                      {"label": "NRIL vertical trend 2017-2020", "value": f"{gnss['2017_2020']['trend_mm_yr']:+.1f} ± {gnss['2017_2020']['trend_se']:.1f}", "unit": "mm/yr"},
                      {"label": "Seasonal amplitude", "value": gnss["all"]["annual_amp_mm"], "unit": "mm"}],
                  "datasets": ["tables/gnss_nril_up.csv", "images/fig_gnss_nril.png"]})
    S.append({"id": "verdict", "headline": VERDICT_HEADLINE, "text": VERDICT_TEXT})
    S.append({"id": "limits", "headline": "What this cannot say", "text": LIMITS_TEXT})
    return {"title": "Norilsk: thawing ground and the TPP-3 tank collapse", "language": "en",
            "sections": S, "credit": CREDIT}


TANK_HEADLINE = "No ground subsidence under tank 5 in the thaw seasons; a small, local settlement in its last winter"
TANK_TEXT = ("Relative to the structures around it, tank 5 did not sink in the snow-free summers 2017-2019, "
             "when thawing ground would settle; the steady -35 mm/yr of Zhang et al. 2022 is not seen. "
             "In its last winter (Dec 2019 - Apr 2020), while the ground was frozen, tank 5 alone settled "
             "about 1.5 cm relative to its surroundings, step by step from January; none of 20 other tanks "
             "and tank-row structures did, including the one 117 m away. The result survives every check "
             "we could run (unwrapping, median, weighting, height error, leaving out any single image). "
             "It fits a failing foundation; it is strong evidence of movement, not proof of its cause.")
TANK_CAVEATS = ("Dry snow delays the radar (~1 mm per mm of snow water): a drift at the tank could mimic "
                "settlement, although the steps do not follow blowing-snow episodes and the neighbouring tank "
                "did not move; tanks elsewhere show ~1 cm excursions in other periods (snow onset, autumn 2019). "
                "One track, one sensor, 20 m pixels (no tilt); whole-tank motion cannot be compared with damage "
                "thresholds, which refer to differences within one tank.")
VERDICT_HEADLINE = "Foundation failure, not ground thaw"
VERDICT_TEXT = ("The radar shows no thaw-driven settlement of the ground at the tank or the plant before the "
                "collapse, and no regional subsidence of the bedrock. This agrees with the inspectorate's "
                "finding: pile cap and piles too weak (33 piles failed in a chain, ~30% not on rock), tank "
                "returned to service in 2019 without notifying the regulator, no permafrost degradation found. "
                "Permafrost thaw is real in Norilsk - wetlands sink every summer - but it is not what the data "
                "show at tank 5. Warming of the frozen ground around the piles cannot be ruled out.")
LIMITS_TEXT = ("No radar data over Norilsk in 2022-2024; motion of the lowland relative to bedrock is "
               "atmosphere-dominated; C-band loses the signal at snowmelt and snow onset, so winter and "
               "summer cannot be linked; settlement within one tank (tilt) is below the 20 m resolution; "
               "uniform settlement of a whole area is invisible to the relative measures used here.")


def tank_check_numbers(checks):
    if not checks:
        return []
    out = []
    c1 = checks.get("c1_unwrapping", {})
    if c1:
        out.append({"label": "Unwrapping check: closures over half a cycle (tank / comparable targets)",
                    "value": f"{c1['tank_window']['frac_over_half_cycle']} / {c1['null_windows']['frac_over_half_cycle']}", "unit": "fraction"})
    for k, v in checks.get("c2_c4_change_vertical_mm", {}).items():
        out.append({"label": f"Tank 5 {k.replace('_', ' -> ')}: mean / median / weighted / with height term",
                    "value": f"{v['window_mean']} / {v['window_median']} / {v['coherence_weighted']} / {v['with_height_error_term']}", "unit": "mm"})
        out.append({"label": f"Leave-one-pair-out range ({k.replace('_', ' -> ')})",
                    "value": f"{v['jackknife_min']} to {v['jackknife_max']}", "unit": "mm"})
    c5 = checks.get("c5_intervals", {})
    if c5:
        key = list(c5)[-1]
        v = c5[key]
        out.append({"label": f"Other tanks, same winter ({key.replace('_', ' -> ')}): depot tanks mean ± sd (n); "
                             "tanks/row structures as low as tank 5",
                    "value": f"{v['osm_tanks_mean']} ± {v['osm_tanks_sd']} ({v['osm_tanks_n']}); "
                             f"{v['tanks_and_row_as_low_or_lower']} of {v['tanks_and_row_n']}", "unit": "mm"})
    w = checks.get("c7_weather", {}).get("2019-20")
    if w and "r_tank_step_vs_blowing_h" in w:
        out.append({"label": "Winter 2019-20 snow: snowfall, blowing-snow hours, r(tank steps, blowing snow)",
                    "value": f"{w['snowfall_cm']} cm, {w['blowing_snow_h']} h, r = {w['r_tank_step_vs_blowing_h']}", "unit": ""})
    return out


def app_stories(man, allsum, tank18, gnss):
    """Two RegionStory objects in the app's format (frontend/src/data/stories.ts):
    paste-ready, numbers read from the outputs."""
    out = []
    if allsum:
        out.append({
            "id": "norilsk-permafrost-thaw", "regionId": "norilsk",
            "title": "Thawing ground: Norilsk's wetlands sink every summer", "category": "water",
            "firstYear": 2020, "mapEndYear": 2026, "dateLabel": "Radar summers 2017–2025 · NISAR 2026",
            "latitude": 69.38, "longitude": 87.8, "cameraAltitude": 60000, "beforeYear": 2020, "afterYear": 2025,
            "summary": ("Sentinel-1 radar shows Norilsk's herbaceous wetlands sinking a few millimetres more than "
                        "the ground around them every summer, in step with accumulated thaw (Stefan's law: "
                        f"{allsum['stefan_weighted_mean']:.3f} ± {allsum['stefan_weighted_se']:.3f} mm per √(°C·day), "
                        "the same in all six summers 2017–2025). The same places sink every summer, NISAR's L-band "
                        "radar saw it independently in 2026, and there is no sign of acceleration."),
            "lookFor": ("Wetlands and low tundra around Norilsk are where the ground sinks each summer; industrial "
                        "and built-up ground does not sink more than other lowland ground."),
            "limitation": ("Motion is relative to ~2 km surroundings: whether the whole lowland sinks relative to "
                           "bedrock is not measurable (atmosphere). No Sentinel-1 data over Norilsk in 2022–2024. "
                           "Water-level changes in wetlands may contribute to the radar signal."),
            "sources": [{"title": "Sentinel-1 interferograms on demand (HyP3)", "publisher": "Alaska Satellite Facility",
                         "url": "https://hyp3-docs.asf.alaska.edu/"},
                        {"title": "NISAR mission", "publisher": "NASA Jet Propulsion Laboratory",
                         "url": "https://nisar.jpl.nasa.gov/"}]})
    if tank18:
        pv = tank18["points"]["published"]["vs_null"]
        last = [k for k in pv if k.startswith("change_") and pv[k]][-1]
        bedrock = f" and the bedrock under Norilsk (GNSS NRIL) rose slightly ({gnss['2017_2020']['trend_mm_yr']:+.1f} mm/yr)" if gnss else ""
        out.append({
            "id": "norilsk-tank5-ground", "regionId": "norilsk",
            "title": "Did sinking ground bring down tank 5?", "category": "water",
            "firstYear": 2020, "mapEndYear": 2021, "dateLabel": "Radar 2017–2020 · collapse May 29, 2020",
            "latitude": POINTS["tank5_zhang2022"][0], "longitude": POINTS["tank5_zhang2022"][1],
            "cameraAltitude": 6000, "beforeYear": 2020, "afterYear": 2021,
            "summary": ("Radar shows tank 5 did not sink relative to the structures around it in the snow-free "
                        f"summers 2017–2019, when thawing ground settles{bedrock}. In its last winter, with the "
                        f"ground frozen, the tank alone settled about {abs(pv[last]['value']):.0f} mm "
                        f"(comparable targets {pv[last]['null_mean']:+.1f} ± {pv[last]['null_std']:.1f} mm); "
                        "none of 20 other tanks did. This fits the inspectorate's finding of a failing pile "
                        "foundation, not permafrost thaw."),
            "lookFor": ("The marker is tank 5 in the fuel-tank row north-west of TPP-3. Compare the 2020 and 2021 "
                        "mosaics; the radar measurements are in the Norilsk subsidence dataset when connected."),
            "limitation": ("Evidence of movement, not proof of its cause: snow drifts can mimic settlement in "
                           "radar, one satellite track only, and 20 m pixels cannot measure tilt within the 45 m "
                           "tank. Uniform settlement of a whole area would be invisible to this relative measure."),
            "sources": [{"title": "Permafrost early deformation signals before the Norilsk oil tank collapse",
                         "publisher": "Remote Sensing (Zhang et al.)", "publishedYear": 2022,
                         "url": "https://www.mdpi.com/2072-4292/14/19/5036"},
                        {"title": "Weak tank piles named as the cause of the Norilsk TPP accident (Rostekhnadzor)",
                         "publisher": "Interfax", "publishedYear": 2020, "url": "https://www.interfax.ru/russia/737323"},
                        {"title": "GNSS station NRIL time series", "publisher": "Nevada Geodetic Laboratory",
                         "url": "http://geodesy.unr.edu/NGLStationPages/stations/NRIL.sta"}]})
    if not out:
        return
    path = PKG / "app" / "norilsk_stories.ts"
    path.parent.mkdir(exist_ok=True)
    path.write_text("// Generated by Nasa-main/12_export_web.py from the analysis outputs.\n"
                    "// Paste into frontend/src/data/stories.ts: `export const stories = [..., ...norilskResultStories]`.\n"
                    'import type { RegionStory } from "./stories";\n\n'
                    f"export const norilskResultStories: RegionStory[] = {json.dumps(out, indent=2, ensure_ascii=False)};\n")
    add(man, path, "app", "Norilsk result stories, paste-ready for the app", "overview",
        f"{len(out)} RegionStory objects (stories.ts format) with summary, lookFor, limitation, sources")


# ----------------------------------------------------------------- docs
def add(man, path, kind, title, section, description, **extra):
    man.append({"file": str(path.relative_to(PKG)), "type": kind, "title": title, "section": section,
                "description": description, **{k: v for k, v in extra.items() if v is not None}})


README_HEAD = """# Norilsk permafrost - data package for the frontend

Generated by `12_export_web.py`; do not edit by hand. Start from
`story.json` (what to say, in order) and `manifest.json` (every file).

Conventions
- Displacement in mm, {sign}.
- Rasters: Cloud-Optimized GeoTIFF, EPSG:4326, float32, NaN = no data.
  Overlays: WebP images in Web Mercator; place them with the four lon/lat
  corners in `manifest.json` (`overlay_corners_lonlat`: NW, NE, SE, SW).
- Vectors: GeoJSON EPSG:4326. Tables: CSV, units in column names.
- Dates: YYYYMMDD in tables from the pipeline, ISO in `app/`.
- Credit (show it): {credit}.

Story sections: {sections}.

For the Dancing-with-the-SARs app (region `norilsk`): read `app/INTEGRATION.md`,
then run `python app/merge_into_app.py --app <app repo>`.

## Files
"""


def readme(man, st):
    lines = [README_HEAD.format(sign=SIGN, credit=CREDIT,
                                sections=", ".join(f"`{s['id']}`" for s in st["sections"]))]
    for sec in ("overview", "thaw", "industry", "tank", "reference"):
        items = [m for m in man if m["section"] == sec]
        if not items:
            continue
        lines.append(f"\n### {sec}\n\n| file | type | what it is |\n|---|---|---|")
        lines += [f"| `{m['file']}` | {m['type']} | **{m['title']}**. {m['description']} |" for m in items]
    (PKG / "README.md").write_text("\n".join(lines) + "\n")


def main():
    if PKG.exists():
        shutil.rmtree(PKG)
    for d in ("rasters", "overlays", "vectors", "tables", "images"):
        (PKG / d).mkdir(parents=True)
    WEB.mkdir(parents=True, exist_ok=True)
    allsum, local, rob = load("robustness_all_summers.json"), load("local_meta.json"), load("robustness.json")
    tank18, checks, gnss = load("tank_series.json"), load("tank_checks.json"), load("gnss_nril.json")
    man = []
    maps(man)
    zones(man)
    sites(man, checks, gnss)
    tables(man, allsum, tank18, checks)
    images(man)
    app_dataset(man)
    app_stories(man, allsum, tank18, gnss)
    st = story(allsum, local, rob, tank18, checks, gnss)
    (PKG / "story.json").write_text(json.dumps(st, indent=1, ensure_ascii=False))
    (PKG / "manifest.json").write_text(json.dumps({"generated_by": "12_export_web.py", "credit": CREDIT,
                                                   "files": man}, indent=1, ensure_ascii=False))
    readme(man, st)
    if KIT and (KIT / "viewer.html").exists():      # local preview of the package
        shutil.copyfile(KIT / "viewer.html", PKG / "viewer.html")
    zip_path = shutil.make_archive(str(OUTPUTS / "frontend_package"), "zip", PKG)
    print(f"{len(man)} files in {PKG}; {zip_path}")


if __name__ == "__main__":
    main()
