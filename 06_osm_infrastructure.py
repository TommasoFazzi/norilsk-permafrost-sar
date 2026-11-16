"""Step 06: OpenStreetMap land use (industrial / quarry / landfill) -> mask.

WorldCover maps most of Norilsk's industry, pits and tailings as "bare",
not "built-up". OSM land-use polygons fill that gap: they are excluded
from both the InSAR reference and the tundra control masks (step 08).
Ways and multipolygon relations (outer minus inner rings) are used.

    python 06_osm_infrastructure.py            -> osm_landuse.geojson, masks,
                                                  osm_mineshafts.geojson
    python 06_osm_infrastructure.py --tanks    -> storage tanks / power plants
"""
import argparse
import json

import numpy as np
import rasterio
import requests
from pyproj import Transformer
from rasterio.features import rasterize
from shapely.geometry import LineString, Point, Polygon, mapping
from shapely.ops import polygonize, transform, unary_union  # noqa: F401

from config import BBOX, DIRS, UTM_EPSG

URL = "https://overpass-api.de/api/interpreter"
s, w, n, e = BBOX[1], BBOX[0], BBOX[3], BBOX[2]
LANDUSE = "industrial|quarry|landfill"
Q_LANDUSE = f"""
[out:json][timeout:300];
(way["landuse"~"^({LANDUSE})$"]({s},{w},{n},{e});
 relation["landuse"~"^({LANDUSE})$"]({s},{w},{n},{e}););
out geom tags;
"""
Q_TANKS = f"""
[out:json][timeout:300];
(nwr["man_made"="storage_tank"]({s},{w},{n},{e});
 nwr["power"="plant"]({s},{w},{n},{e}););
out center tags;
"""
Q_MINES = f"""
[out:json][timeout:300];
nwr["man_made"="mineshaft"]({s},{w},{n},{e});
out center tags;
"""
BUFFER_M = 100  # ~2-3 InSAR pixels: edges, haul roads, misregistration


def overpass(q):
    r = requests.post(URL, data={"data": q}, timeout=600,
                      headers={"User-Agent": "norilsk-insar-study/0.1"})
    r.raise_for_status()
    return r.json()["elements"]


def ring(coords):
    return [(p["lon"], p["lat"]) for p in coords]


def element_polygon(el):
    if el["type"] == "way":
        c = ring(el.get("geometry", []))
        return Polygon(c) if len(c) >= 4 and c[0] == c[-1] else None
    lines = {"outer": [], "inner": []}
    for m in el.get("members", []):
        if m["type"] == "way" and m.get("geometry"):
            lines.setdefault(m.get("role") or "outer", []).append(
                LineString(ring(m["geometry"])))
    outer = unary_union(list(polygonize(lines["outer"])))
    inner = unary_union(list(polygonize(lines["inner"])))
    return None if outer.is_empty else outer.difference(inner)


def landuse():
    feats, geoms = [], []
    for el in overpass(Q_LANDUSE):
        g = element_polygon(el)
        if g is None or g.is_empty:
            continue
        g = g.buffer(0)
        tags = el.get("tags", {})
        feats.append({"type": "Feature", "geometry": mapping(g),
                      "properties": {"osm_id": f"{el['type']}/{el['id']}",
                                     "landuse": tags.get("landuse"),
                                     "name": tags.get("name")}})
        geoms.append(g)
    out = DIRS["osm"] / "osm_landuse.geojson"
    out.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))

    # rasterize on the step-05 grid, buffered in metres
    with rasterio.open(DIRS["dem"] / "dem_utm30.tif") as src:
        prof, shape, tr_px = src.profile, src.shape, src.transform
    to_utm = Transformer.from_crs(4326, UTM_EPSG, always_xy=True).transform
    utm = [transform(to_utm, g).buffer(BUFFER_M) for g in geoms]
    mask = (rasterize([(g, 1) for g in utm], out_shape=shape, transform=tr_px,
                      dtype="uint8") if utm else np.zeros(shape, "uint8"))
    prof.update(dtype="uint8", nodata=None)
    with rasterio.open(DIRS["osm"] / "osm_landuse_mask.tif", "w", **prof) as d:
        d.write(mask, 1)
    # unbuffered version: zone statistics (industrial zone) use this one
    core = (rasterize([(transform(to_utm, g), 1) for g in geoms], out_shape=shape,
                      transform=tr_px, dtype="uint8")
            if geoms else np.zeros(shape, "uint8"))
    with rasterio.open(DIRS["osm"] / "osm_landuse_core.tif", "w", **prof) as d:
        d.write(core, 1)
    area = {}
    for f, g in zip(feats, utm):
        k = f["properties"]["landuse"]
        area[k] = area.get(k, 0) + g.area / 1e6
    print(f"{len(feats)} polygons -> {out.name}")
    print("area incl. buffer (km2):", {k: round(v, 1) for k, v in area.items()})
    print(f"mask: {mask.sum() * 0.03 ** 2:.1f} km2 (overlaps merged)")


def points(query, fname):
    feats = []
    for el in overpass(query):
        c = el.get("center") or {"lat": el.get("lat"), "lon": el.get("lon")}
        if c["lat"] is None:
            continue
        tags = el.get("tags", {})
        feats.append({"type": "Feature",
                      "geometry": mapping(Point(c["lon"], c["lat"])),
                      "properties": {"osm_id": f"{el['type']}/{el['id']}", **tags}})
    out = DIRS["osm"] / fname
    out.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    print(f"{len(feats)} features -> {out.name}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tanks", action="store_true")
    if ap.parse_args().tanks:
        points(Q_TANKS, "osm_tanks_plants.geojson")
    else:
        landuse()
        # underground mines (Talnakh): shafts outline the mined area
        points(Q_MINES, "osm_mineshafts.geojson")
