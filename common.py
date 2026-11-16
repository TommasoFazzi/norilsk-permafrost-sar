"""Shared helpers: the 30 m UTM analysis grid, zones, reprojection, dates."""
import json
import re

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject

from config import DIRS

ZONE_NAMES = {1: "reference", 2: "industrial", 3: "accident", 4: "urban",
              5: "control_NW", 6: "control_SE"}


def grid():
    """Profile of the step-05 grid (UTM 45N, 30 m): every layer is put on it."""
    with rasterio.open(DIRS["dem"] / "dem_utm30.tif") as s:
        return s.profile.copy()


def read(path, band=1):
    with rasterio.open(path) as s:
        return s.read(band).astype("float32"), s.profile


def zones():
    z, _ = read(DIRS["zones"] / "zones_utm30.tif")
    ref, _ = read(DIRS["zones"] / "ref_patches_utm30.tif")
    slope, _ = read(DIRS["dem"] / "slope_deg_utm30.tif")
    return z.astype("uint8"), ref.astype("uint8"), slope


def to_grid(arr, src_transform, src_crs, prof, resampling=Resampling.bilinear):
    """Reproject a float array onto the analysis grid (NaN = no data)."""
    out = np.full((prof["height"], prof["width"]), np.nan, "float32")
    reproject(arr.astype("float32"), out, src_transform=src_transform, src_crs=src_crs,
              dst_transform=prof["transform"], dst_crs=prof["crs"],
              src_nodata=np.nan, dst_nodata=np.nan, resampling=resampling)
    return out


def pair_dates(name):
    """First two YYYYMMDD groups in a product name -> (ref, sec) strings."""
    d = re.findall(r"(?<!\d)(20\d{6})(?:T\d+)?(?!\d)", name)
    if len(d) < 2:
        raise ValueError(f"no pair dates in {name}")
    return d[0], d[1]


def write_stack(path, bands, prof, descriptions=None, nodata=np.nan):
    p = dict(prof, count=len(bands), dtype="float32", nodata=nodata,
             compress="deflate", driver="GTiff")
    with rasterio.open(path, "w", **p) as d:
        for i, b in enumerate(bands, 1):
            d.write(np.asarray(b, "float32"), i)
            if descriptions:
                d.set_band_description(i, descriptions[i - 1])


def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=1, default=str))
