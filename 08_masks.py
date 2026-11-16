"""Step 08: analysis zones, v3 (needs 05 + 06).

Writes data/zones/zones_utm30.tif (codes below), ref_patches_utm30.tif
(1..N, one id per reference patch), zones_points.json (one representative
point per zone) and outputs/masks_v3.png.

Zone codes
  1 reference    plateau bedrock: WorldCover bare (60), elevation > 300 m,
                 slope < 10 deg, > MINE_BUFFER_KM from the convex hulls of
                 the mine-shaft clusters (OSM man_made=mineshaft) and
                 > 3 km from OSM industrial/quarry/landfill, east of
                 REF_MIN_LON; the N_REF largest patches (>= 0.1 km2) at
                 least REF_SEP_KM apart are kept, so they can be checked
                 against each other
  2 industrial   OSM landuse industrial/quarry/landfill (unbuffered),
                 minus the accident zone
  3 accident     TPP-3 site: ACCIDENT_RADIUS_M around ACCIDENT
  4 urban        WorldCover built-up (50) outside zones 2-3
  5 control_NW   tundra rule v2 (below), patches west of NW_MAX_LON
  6 control_SE   tundra rule v2, patches east of SE_MIN_LON and south of
                 SE_MAX_LAT
Tundra rule v2 (approved): WorldCover shrub/grass/moss-lichen (20/30/100),
no wetland/water/trees; elevation and slope within the urban-terrace range
(built-up p5-p95 elevation, p90 slope); >= HAND_MIN m above the nearest
water body > 1 km2 (not floodplain); > 2 km from built-up or OSM land use
(buffered); patches >= 0.5 km2.

    python 08_masks.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import ListedColormap
from pyproj import Transformer
from scipy import ndimage as nd
from shapely.geometry import MultiPoint, Point
from shapely.ops import transform as sh_transform
from shapely.ops import unary_union

from config import DATA, DIRS, OUTPUTS, POINTS, UTM_EPSG

RES = 30.0
PX_KM2 = (RES / 1000) ** 2
BUILT, BARE, WATER = 50, 60, 80
TUNDRA_CLASSES = [20, 30, 100]
HAND_MIN = 5.0
BIG_WATER_KM2 = 1.0
TUNDRA_BUILT_DIST_KM = 2.0
MINE_BUFFER_KM = 3.0
SHAFT_LINK_KM = 2.5
REF_OSM_DIST_KM = 3.0
N_REF = 3
REF_MIN_LON = 88.60
REF_SEP_KM = 2.0
ACCIDENT = (69.3227, 87.9499)   # OSM way/525598517 "Норильская ТЭЦ-3" centroid
ACCIDENT_RADIUS_M = 1000
NW_MAX_LON = 88.05
SE_MIN_LON = 88.30
SE_MAX_LAT = 69.42
ZONES = {1: "reference", 2: "industrial", 3: "accident", 4: "urban",
         5: "control_NW", 6: "control_SE"}
COLORS = {1: "red", 2: "orange", 3: "magenta", 4: "yellow", 5: "green", 6: "cyan"}
ZDIR = DATA / "zones"
OUT = OUTPUTS


def read(path):
    with rasterio.open(path) as s:
        return s.read(1), s.profile


def write(path, arr, prof):
    with rasterio.open(path, "w", **dict(prof, dtype="uint8", nodata=0)) as d:
        d.write(arr.astype("uint8"), 1)


def min_patch(mask, km2):
    lab, n = nd.label(mask)
    if n == 0:
        return mask
    size = nd.sum(mask, lab, range(1, n + 1)) * PX_KM2
    return np.isin(lab, 1 + np.where(size >= km2)[0])


def grids(prof):
    h, w = prof["height"], prof["width"]
    cols, rows = np.meshgrid(np.arange(w) + 0.5, np.arange(h) + 0.5)
    x, y = prof["transform"] * (cols, rows)
    lon, lat = Transformer.from_crs(UTM_EPSG, 4326, always_xy=True).transform(x, y)
    return x, y, lon, lat


def km_from(mask):
    return nd.distance_transform_edt(~mask) * RES / 1000


def centre_point(mask, lon, lat):
    """Pixel closest to the centroid of the largest patch -> (lat, lon)."""
    lab, n = nd.label(mask)
    if n > 1:
        mask = lab == 1 + int(np.argmax(nd.sum(mask, lab, range(1, n + 1))))
    rr, cc = np.where(mask)
    k = np.argmin((rr - rr.mean()) ** 2 + (cc - cc.mean()) ** 2)
    return round(float(lat[rr[k], cc[k]]), 4), round(float(lon[rr[k], cc[k]]), 4)


def main():
    dem, prof = read(DIRS["dem"] / "dem_utm30.tif")
    slope, _ = read(DIRS["dem"] / "slope_deg_utm30.tif")
    lc, _ = read(DIRS["landcover"] / "worldcover_utm30.tif")
    osm_buf = read(DIRS["osm"] / "osm_landuse_mask.tif")[0].astype(bool)
    osm_core = read(DIRS["osm"] / "osm_landuse_core.tif")[0].astype(bool)
    valid = dem > 0
    x, y, lon, lat = grids(prof)
    to_utm = Transformer.from_crs(4326, UTM_EPSG, always_xy=True)

    # ---- accident zone
    ax_, ay_ = to_utm.transform(ACCIDENT[1], ACCIDENT[0])
    accident = np.hypot(x - ax_, y - ay_) <= ACCIDENT_RADIUS_M

    # ---- mine exclusion: one convex hull per shaft cluster (shafts closer
    # than 2*SHAFT_LINK_KM join), + MINE_BUFFER_KM. Norilsk and Talnakh mines
    # form separate clusters.
    shafts = json.loads((DIRS["osm"] / "osm_mineshafts.geojson").read_text())["features"]
    sp = [sh_transform(to_utm.transform, Point(f["geometry"]["coordinates"]))
          for f in shafts]
    blobs = unary_union([q.buffer(SHAFT_LINK_KM * 1000) for q in sp])
    blobs = list(getattr(blobs, "geoms", [blobs]))
    hulls = [MultiPoint([q for q in sp if b.contains(q)]).convex_hull for b in blobs]
    hull_buf = unary_union([h.buffer(MINE_BUFFER_KM * 1000) for h in hulls])
    from rasterio.features import rasterize
    near_mine = rasterize([(hull_buf, 1)], out_shape=dem.shape,
                          transform=prof["transform"], dtype="uint8").astype(bool)

    # ---- reference patches
    ref_all = (valid & (lc == BARE) & (dem > 300) & (slope < 10) & ~near_mine
               & (km_from(osm_buf) > REF_OSM_DIST_KM))
    ref_all = min_patch(ref_all, 0.1)
    lab, n = nd.label(ref_all)
    if n == 0:
        raise SystemExit("No reference patch found: relax the rule.")
    idx = range(1, n + 1)
    area = nd.sum(ref_all, lab, idx)
    plon, cx, cy = (np.array(nd.mean(v, lab, idx)) for v in (lon, x, y))
    keep = []   # largest patches east of REF_MIN_LON, >= REF_SEP_KM apart
    for k in np.argsort(-area):
        if plon[k] < REF_MIN_LON:
            continue
        if all(np.hypot(cx[k] - cx[j], cy[k] - cy[j]) >= REF_SEP_KM * 1000 for j in keep):
            keep.append(k)
        if len(keep) == N_REF:
            break
    keep = [k + 1 for k in keep]
    ref_patches = np.zeros(dem.shape, "uint8")
    for i, k in enumerate(keep, 1):
        ref_patches[lab == k] = i
    ref = ref_patches > 0

    # ---- tundra rule v2
    urban_core = min_patch(nd.binary_closing(lc == BUILT, iterations=3), 0.5) & (lc == BUILT)
    z_lo, z_hi = np.percentile(dem[urban_core], [5, 95])
    s_hi = np.percentile(slope[urban_core], 90)
    big_water = min_patch(lc == WATER, BIG_WATER_KM2)
    _, (iy, ix) = nd.distance_transform_edt(~big_water, return_indices=True)
    hand = dem - dem[iy, ix]
    tundra = (valid & np.isin(lc, TUNDRA_CLASSES) & (dem >= z_lo) & (dem <= z_hi)
              & (slope <= s_hi) & (hand >= HAND_MIN) & ~osm_buf
              & (km_from((lc == BUILT) | osm_buf) > TUNDRA_BUILT_DIST_KM))
    tundra = min_patch(tundra, 0.5)
    tl, tn = nd.label(tundra)
    plon = nd.mean(lon, tl, range(1, tn + 1))
    plat = nd.mean(lat, tl, range(1, tn + 1))
    nw = np.isin(tl, 1 + np.where(plon < NW_MAX_LON)[0])
    se = np.isin(tl, 1 + np.where((plon > SE_MIN_LON) & (plat < SE_MAX_LAT))[0])

    # ---- zone raster (later codes do not overwrite earlier ones)
    zones = np.zeros(dem.shape, "uint8")
    for code, m in [(3, accident & valid), (2, osm_core), (1, ref),
                    (4, lc == BUILT), (5, nw), (6, se)]:
        zones[(zones == 0) & m] = code

    ZDIR.mkdir(exist_ok=True)
    write(ZDIR / "zones_utm30.tif", zones, prof)
    write(ZDIR / "ref_patches_utm30.tif", ref_patches, prof)
    # all bedrock candidates (any longitude): fallback reference for products
    # that do not cover the chosen patches (e.g. NISAR track 6 DESC)
    write(ZDIR / "ref_candidates_utm30.tif", ref_all, prof)

    pts = {ZONES[c]: centre_point(zones == c, lon, lat) for c in ZONES if (zones == c).any()}
    for i in range(1, ref_patches.max() + 1):
        pts[f"ref_{i}"] = centre_point(ref_patches == i, lon, lat)
    (ZDIR / "zones_points.json").write_text(json.dumps(pts, indent=1))

    # ---- report
    print(f"urban terraces: elev {z_lo:.0f}-{z_hi:.0f} m, slope p90 {s_hi:.1f} deg")
    print(f"{len(hulls)} mine-shaft clusters from {len(shafts)} shafts, "
          f"buffer {MINE_BUFFER_KM} km")
    print(f"reference candidates: {n} patches, kept {len(keep)} largest east of {REF_MIN_LON}")
    for i in range(1, ref_patches.max() + 1):
        m = ref_patches == i
        print(f"  ref_{i}: {m.sum() * PX_KM2:.2f} km2, elev {dem[m].mean():.0f} m, "
              f"slope {slope[m].mean():.1f}, centre {pts[f'ref_{i}']}")
    for c, name in ZONES.items():
        m = zones == c
        print(f"{c} {name:11s} {m.sum() * PX_KM2:7.1f} km2  "
              f"elev {np.nanmean(np.where(m, dem, np.nan)):5.0f} m  "
              f"slope {np.nanmean(np.where(m, slope, np.nan)):4.1f}  pt {pts.get(name)}")

    # ---- figure
    OUT.mkdir(exist_ok=True)
    b = rasterio.transform.array_bounds(prof["height"], prof["width"], prof["transform"])
    ext = (b[0], b[2], b[1], b[3])
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.imshow(np.ma.masked_where(~valid, dem), extent=ext, cmap="Greys_r")
    for c, col in COLORS.items():
        m = zones == c
        ax.imshow(np.ma.masked_where(~m, m), extent=ext, cmap=ListedColormap([col]),
                  alpha=0.8, interpolation="nearest")
    for g in getattr(hull_buf, "geoms", [hull_buf]):
        ax.plot(*g.exterior.xy, "--", color="red", lw=1)
    for h in hulls:
        if h.geom_type == "Polygon":
            ax.plot(*h.exterior.xy, "-", color="darkred", lw=0.8)
    for f in shafts:
        px, py = to_utm.transform(*f["geometry"]["coordinates"])
        ax.plot(px, py, "x", color="darkred", ms=4)
    for name, (la, lo) in pts.items():
        px, py = to_utm.transform(lo, la)
        ax.plot(px, py, "o", ms=6, mfc="white", mec="black")
        ax.annotate(name, (px, py), xytext=(4, 4), textcoords="offset points", fontsize=8,
                    bbox=dict(fc="white", alpha=0.75, lw=0, pad=1))
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in ZONES]
    ax.legend(handles, [ZONES[c] for c in ZONES], loc="lower left", fontsize=8)
    ax.set_title("Zones v3. x = OSM mine shafts; solid = shaft hull; "
                 f"dashed = hull + {MINE_BUFFER_KM:.0f} km (reference excluded)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "masks_v3.png", dpi=130)


if __name__ == "__main__":
    main()
