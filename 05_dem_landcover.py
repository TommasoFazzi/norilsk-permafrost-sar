"""Step 05: Copernicus DEM (30 m) + ESA WorldCover (10 m) on one UTM 30 m grid.

Public AWS buckets, no login. Also derives slope and two candidate masks
for the experimental design:
  ref_candidate    flat bare rock on high ground  -> InSAR reference
  tundra_candidate natural tundra/wetland, >2 km from built-up -> control
Thresholds are first guesses: inspect the masks in QGIS and adjust.

    python 05_dem_landcover.py
"""
import numpy as np
import rasterio
import requests
from rasterio.merge import merge
from rasterio.warp import Resampling, calculate_default_transform, reproject
from scipy.ndimage import distance_transform_edt

from config import BBOX, DIRS, UTM_EPSG

DEM_URL = ("https://copernicus-dem-30m.s3.amazonaws.com/"
           "Copernicus_DSM_COG_10_N{lat:02d}_00_E{lon:03d}_00_DEM/"
           "Copernicus_DSM_COG_10_N{lat:02d}_00_E{lon:03d}_00_DEM.tif")
WC_URL = ("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
          "ESA_WorldCover_10m_2021_v200_N{lat:02d}E{lon:03d}_Map.tif")
RES = 30.0
# WorldCover classes
BUILT, BARE, WATER = 50, 60, 80
NATURAL = [20, 30, 90, 95, 100]  # shrub, grass, wetland, mangrove, moss/lichen


def get(url, dst):
    if dst.exists():
        return dst
    print(f"  downloading {url.rsplit('/', 1)[-1]}")
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(dst, "wb") as f:
            for c in r.iter_content(1 << 20):
                f.write(c)
    return dst


def tiles_1deg(b):
    return [(la, lo) for la in range(int(np.floor(b[1])), int(np.floor(b[3])) + 1)
            for lo in range(int(np.floor(b[0])), int(np.floor(b[2])) + 1)]


def tiles_3deg(b):
    return sorted({(3 * int(la // 3), 3 * int(lo // 3)) for la, lo in tiles_1deg(b)})


def target_grid():
    """Fixed UTM 30 m grid covering BBOX, shared by every layer."""
    dtr, dw, dh = calculate_default_transform(
        "EPSG:4326", UTM_EPSG, 1000, 1000, *BBOX, resolution=RES)
    return dtr, dw, dh


def to_utm(srcs, resampling, dst_path, dtype):
    mosaic, tr = merge(srcs, bounds=BBOX)
    dtr, dw, dh = target_grid()
    out = np.zeros((dh, dw), dtype=dtype)
    reproject(mosaic[0], out, src_transform=tr, src_crs=srcs[0].crs,
              dst_transform=dtr, dst_crs=UTM_EPSG, resampling=resampling)
    prof = dict(driver="GTiff", height=dh, width=dw, count=1, dtype=dtype,
                crs=f"EPSG:{UTM_EPSG}", transform=dtr, compress="deflate")
    with rasterio.open(dst_path, "w", **prof) as d:
        d.write(out, 1)
    return out, prof


def write(path, arr, prof):
    p = dict(prof, dtype=arr.dtype.name)
    with rasterio.open(path, "w", **p) as d:
        d.write(arr, 1)


def main():
    dd, lc_dir = DIRS["dem"], DIRS["landcover"]
    dem_files = [get(DEM_URL.format(lat=la, lon=lo), dd / f"cop30_N{la}E{lo}.tif")
                 for la, lo in tiles_1deg(BBOX)]
    wc_files = [get(WC_URL.format(lat=la, lon=lo), lc_dir / f"wc2021_N{la}E{lo}.tif")
                for la, lo in tiles_3deg(BBOX)]

    dem, prof = to_utm([rasterio.open(f) for f in dem_files], Resampling.bilinear,
                       dd / "dem_utm30.tif", "float32")
    gy, gx = np.gradient(dem, RES)
    slope = np.degrees(np.arctan(np.hypot(gx, gy))).astype("float32")
    write(dd / "slope_deg_utm30.tif", slope, prof)

    # land cover: majority class of 10 m pixels -> 30 m
    lc, _ = to_utm([rasterio.open(f) for f in wc_files], Resampling.mode,
                   lc_dir / "worldcover_utm30.tif", "uint8")

    built_dist_km = distance_transform_edt(lc != BUILT) * RES / 1000
    ref = (lc == BARE) & (dem > 300) & (slope < 10)
    tundra = np.isin(lc, NATURAL) & (dem < 200) & (slope < 5) & (built_dist_km > 2)
    write(lc_dir / "built_distance_km.tif", built_dist_km.astype("float32"), prof)
    write(lc_dir / "ref_candidate.tif", ref.astype("uint8"), prof)
    write(lc_dir / "tundra_candidate.tif", tundra.astype("uint8"), prof)

    px_km2 = (RES / 1000) ** 2
    print(f"DEM {dem.min():.0f}-{dem.max():.0f} m; grid {dem.shape}")
    print(f"Built-up: {(lc == BUILT).sum() * px_km2:.1f} km2")
    print(f"Reference candidate: {ref.sum() * px_km2:.1f} km2")
    print(f"Tundra control candidate: {tundra.sum() * px_km2:.1f} km2")


if __name__ == "__main__":
    main()
