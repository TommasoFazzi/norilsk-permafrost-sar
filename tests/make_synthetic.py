"""Synthetic end-to-end test data for steps 09-12 (no login needed).

Builds a fake data folder with the real static layers (DEM, land cover,
zones, climate: symlinked) plus
  - HyP3-like multi-burst products (unw_phase, corr, lv_theta GeoTIFFs,
    UTM 45N, 40 m) for the planned 2020 + 2025 pairs
  - NISAR-GUNW-like HDF5 files (4 pairs: 111 A / 6 D, early / late 2026)
with a KNOWN signal: vertical(t) = B[zone] * sqrt(TDD_cum(t)) (mm), plus an
east motion EAST_MM in the industrial zone for NISAR, plus noise.
Then steps 09-12 must recover B (see tests/check_synthetic.py).

    python tests/make_synthetic.py OUT_DIR
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import ZONE_NAMES, to_grid  # noqa: E402
from config import S1_WAVELENGTH, ZONE_CLIMATE  # noqa: E402

REAL = Path(__file__).resolve().parents[1] / "data"
B = {"reference": 0.0, "industrial": -0.8, "accident": -1.2, "urban": -0.6,
     "control_NW": -0.4, "control_SE": -0.25}
B_OUTSIDE = -0.3
EAST_MM = 3.0            # industrial zone, NISAR windows only
INC_DEG = 39.0
ATM_MM = 2.0             # per-date smooth atmospheric noise (std)
PIX_MM = 1.0             # per-pair pixel noise (std)
NISAR_WL = 299_792_458.0 / 1.2575e9
rng = np.random.default_rng(42)


def tdd(point):
    df = pd.read_csv(REAL / "climate" / f"era5land_daily_{point}.csv",
                     index_col="date", parse_dates=["date"])
    return df["tdd_cum"]


def b_map(z):
    b = np.full(z.shape, B_OUTSIDE, "float32")
    for c, name in ZONE_NAMES.items():
        b[z == c] = B[name]
    return b


def sqrt_tdd_map(z, date):
    """sqrt(TDD) per pixel from the zone's climate point (norilsk_centre outside)."""
    t = pd.Timestamp(date)
    out = np.full(z.shape, np.sqrt(tdd("norilsk_centre")[t]), "float32")
    for c, name in ZONE_NAMES.items():
        out[z == c] = np.sqrt(tdd(ZONE_CLIMATE[name])[t])
    return out


def atm(shape):
    f = gaussian_filter(rng.standard_normal(shape), 30)
    return (f / f.std() * ATM_MM).astype("float32")


def main(out):
    out = Path(out)
    for k in ("dem", "landcover", "osm", "zones", "climate", "s1"):
        (out / k).parent.mkdir(parents=True, exist_ok=True)
        if (out / k).is_dir() and not (out / k).is_symlink() and not any((out / k).iterdir()):
            (out / k).rmdir()
        if not (out / k).exists():
            (out / k).symlink_to(REAL / k)
    with rasterio.open(REAL / "zones" / "zones_utm30.tif") as s:
        z30, p30 = s.read(1), s.profile
    with rasterio.open(REAL / "landcover" / "worldcover_utm30.tif") as s:
        lc30 = s.read(1)
    b30 = b_map(z30)

    # ---- S1 HyP3-like products on a 40 m UTM grid covering the 30 m grid
    tr30 = p30["transform"]
    W, N = tr30.c, tr30.f
    w40 = int(p30["width"] * 30 / 40) + 1
    h40 = int(p30["height"] * 30 / 40) + 1
    tr40 = from_origin(W, N, 40, 40)
    p40 = dict(driver="GTiff", width=w40, height=h40, count=1, dtype="float32",
               crs=p30["crs"], transform=tr40, nodata=0, compress="deflate")
    g40 = dict(p40, height=h40, width=w40)

    def on40(a):
        return to_grid(a, tr30, p30["crs"], g40)

    pairs = pd.read_csv(REAL / "hyp3" / "pairs_p136.csv", parse_dates=["ref_date", "sec_date"])
    dates = sorted(set(pairs.ref_date) | set(pairs.sec_date))
    vert = {d: b30 * sqrt_tdd_map(z30, d) + atm(z30.shape) for d in dates}
    water = on40((lc30 == 80).astype("float32")) > 0.5
    cos_inc = np.cos(np.radians(INC_DEG))
    pdir = out / "hyp3" / "norilsk_p136"
    for _, r in pairs.iterrows():
        d0, d1 = r.ref_date.strftime("%Y%m%d"), r.sec_date.strftime("%Y%m%d")
        name = f"S1_136_000000s1n00-136231s2n02-136231s3n02_IW_{d0}_{d1}_VV_INT40_SYNT"
        f = pdir / name
        f.mkdir(parents=True, exist_ok=True)
        dlos = on40((vert[r.sec_date] - vert[r.ref_date]) * cos_inc)
        dlos += rng.standard_normal(dlos.shape).astype("float32") * PIX_MM
        phase = -dlos / 1000 * 4 * np.pi / S1_WAVELENGTH
        coh = np.where(water, 0.1, 0.6).astype("float32")
        phase[~np.isfinite(phase)] = 0
        for key, arr in [("unw_phase", phase), ("corr", coh),
                         ("lv_theta", np.full(phase.shape, np.radians(90 - INC_DEG), "float32"))]:
            with rasterio.open(f / f"{name}_{key}.tif", "w", **p40) as d:
                d.write(arr.astype("float32"), 1)
    print(f"S1: {len(pairs)} synthetic pairs, {len(dates)} dates -> {pdir}")

    # ---- NISAR GUNW-like HDF5, 80 m grid, y decreasing
    ndir = out / "nisar"
    ndir.mkdir(parents=True, exist_ok=True)
    x80 = W + 40 + 80 * np.arange(int(p30["width"] * 30 / 80))
    y80 = N - 40 - 80 * np.arange(int(p30["height"] * 30 / 80))
    tr80 = from_origin(W, N, 80, 80)
    g80 = dict(p40, width=x80.size, height=y80.size, transform=tr80)
    ind = (z30 == 2).astype("float32")
    los_v = {"A": (-0.60, -0.12), "D": (0.60, -0.12)}
    pairs_n = [("111", "A", "20260621", "20260703"), ("111", "A", "20260820", "20260901"),
               ("006", "D", "20260626", "20260708"), ("006", "D", "20260825", "20260906")]
    cx = W + 2000 * np.arange(int(x80.size * 80 / 2000) + 2)
    cy = N - 2000 * np.arange(int(y80.size * 80 / 2000) + 2)
    for trk, dr, d0, d1 in pairs_n:
        up = b30 * (sqrt_tdd_map(z30, d1) - sqrt_tdd_map(z30, d0)) + atm(z30.shape) * 0.5
        east = ind * EAST_MM
        e, n = los_v[dr]
        u = np.sqrt(1 - e * e - n * n)
        los = to_grid(e * east + u * up, tr30, p30["crs"], g80)
        los += rng.standard_normal(los.shape).astype("float32") * PIX_MM
        phase = SIGN_NISAR * los / 1000 * 4 * np.pi / NISAR_WL
        phase[~np.isfinite(phase)] = 0
        name = (f"NISAR_L2_PR_GUNW_023_{trk}_{dr}_037_024_2000_SH_{d0}T212950_{d0}T213027_"
                f"{d1}T212950_{d1}T213027_P05023_N_F_J_001.h5")
        with h5py.File(ndir / name, "w") as h:
            g = "science/LSAR/GUNW/grids/frequencyA"
            h[f"{g}/centerFrequency"] = 1.2575e9
            ui = f"{g}/unwrappedInterferogram"
            h[f"{ui}/xCoordinates"] = x80
            h[f"{ui}/yCoordinates"] = y80
            pj = h.create_dataset(f"{ui}/projection", data=32645)
            pj.attrs["epsg_code"] = 32645
            h[f"{ui}/HH/unwrappedPhase"] = phase.astype("float32")
            h[f"{ui}/HH/coherenceMagnitude"] = np.full(phase.shape, 0.7, "float32")
            rg = "science/LSAR/GUNW/metadata/radarGrid"
            h[f"{rg}/xCoordinates"] = cx
            h[f"{rg}/yCoordinates"] = cy
            h[f"{rg}/heightAboveEllipsoid"] = np.array([-500., 0., 500., 1000.])
            h[f"{rg}/losUnitVectorX"] = np.full((4, cy.size, cx.size), e, "float32")
            h[f"{rg}/losUnitVectorY"] = np.full((4, cy.size, cx.size), n, "float32")
    print(f"NISAR: {len(pairs_n)} synthetic GUNW -> {ndir}")


SIGN_NISAR = -1.0   # phase = SIGN * LOS * 4pi/lambda, inverse of 10_nisar.py

if __name__ == "__main__":
    main(sys.argv[1])
