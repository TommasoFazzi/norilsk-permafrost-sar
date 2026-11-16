"""Quick-look PNGs of step 05 outputs with config.POINTS overlaid.

    python quicklook.py        -> outputs/*.png + point/mask diagnostics
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import BoundaryNorm, ListedColormap
from pyproj import Transformer

from config import DIRS, POINTS, UTM_EPSG

OUT = Path(__file__).parent / "outputs"
OUT.mkdir(exist_ok=True)

WC = {0: ("nodata", "#ffffff"), 10: ("tree", "#006400"), 20: ("shrub", "#ffbb22"),
      30: ("grass", "#ffff4c"), 40: ("crop", "#f096ff"), 50: ("built", "#fa0000"),
      60: ("bare", "#b4b4b4"), 70: ("snow/ice", "#f0f0f0"), 80: ("water", "#0064c8"),
      90: ("wetland", "#0096a0"), 95: ("mangrove", "#00cf75"),
      100: ("moss/lichen", "#fae6a0")}


def load(path):
    with rasterio.open(path) as s:
        return s.read(1), s.transform, s.bounds


def to_px(transform):
    tr = Transformer.from_crs(4326, UTM_EPSG, always_xy=True)
    out = {}
    for k, (lat, lon) in POINTS.items():
        x, y = tr.transform(lon, lat)
        c, r = ~transform * (x, y)
        out[k] = (x, y, int(r), int(c))
    return out


def plot(arr, bounds, pts, title, fname, **kw):
    fig, ax = plt.subplots(figsize=(9, 8))
    ext = (bounds.left, bounds.right, bounds.bottom, bounds.top)
    im = ax.imshow(arr, extent=ext, interpolation="nearest", **kw)
    for k, (x, y, *_ ) in pts.items():
        ax.plot(x, y, "o", ms=7, mfc="magenta", mec="black")
        ax.annotate(k, (x, y), xytext=(5, 5), textcoords="offset points",
                    fontsize=8, color="black",
                    bbox=dict(fc="white", alpha=0.7, lw=0, pad=1))
    ax.set_title(title)
    ax.set_xlabel(f"UTM E (EPSG:{UTM_EPSG}) [m]")
    ax.set_ylabel("UTM N [m]")
    fig.colorbar(im, ax=ax, shrink=0.7)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=130)
    plt.close(fig)


def main():
    dem, tr, b = load(DIRS["dem"] / "dem_utm30.tif")
    slope, *_ = load(DIRS["dem"] / "slope_deg_utm30.tif")
    lc, *_ = load(DIRS["landcover"] / "worldcover_utm30.tif")
    ref, *_ = load(DIRS["landcover"] / "ref_candidate.tif")
    tun, *_ = load(DIRS["landcover"] / "tundra_candidate.tif")
    bd, *_ = load(DIRS["landcover"] / "built_distance_km.tif")
    pts = to_px(tr)
    nod = dem == 0

    plot(np.ma.masked_where(nod, dem), b, pts, "Copernicus DEM 30 m [m]",
         "dem.png", cmap="terrain")
    plot(np.ma.masked_where(nod, slope), b, pts, "Slope [deg]", "slope.png",
         cmap="magma", vmax=30)
    codes = sorted(WC)
    cmap = ListedColormap([WC[c][1] for c in codes])
    norm = BoundaryNorm(codes + [101], cmap.N)
    plot(lc, b, pts, "ESA WorldCover 2021 (30 m mode)", "landcover.png",
         cmap=cmap, norm=norm)
    hill = np.ma.masked_where(nod, dem)
    for m, name, col in [(ref, "ref_candidate", "red"),
                         (tun, "tundra_candidate", "green")]:
        fig, ax = plt.subplots(figsize=(9, 8))
        ext = (b.left, b.right, b.bottom, b.top)
        ax.imshow(hill, extent=ext, cmap="Greys_r")
        ax.imshow(np.ma.masked_where(m == 0, m), extent=ext,
                  cmap=ListedColormap([col]), alpha=0.8, interpolation="nearest")
        ax.imshow(np.ma.masked_where(lc != 50, lc), extent=ext,
                  cmap=ListedColormap(["yellow"]), alpha=0.9,
                  interpolation="nearest")
        for k, (x, y, *_ ) in pts.items():
            ax.plot(x, y, "o", ms=7, mfc="magenta", mec="black")
            ax.annotate(k, (x, y), xytext=(5, 5), textcoords="offset points",
                        fontsize=8, bbox=dict(fc="white", alpha=0.7, lw=0, pad=1))
        ax.set_title(f"{name} ({col}) over DEM; WorldCover built-up in yellow")
        fig.tight_layout()
        fig.savefig(OUT / f"{name}.png", dpi=130)
        plt.close(fig)

    # ---- diagnostics
    px = 0.03 ** 2
    print(f"nodata (reprojection corners): {nod.mean() * 100:.1f}% of grid")
    print(f"built-up {np.sum(lc == 50) * px:.1f} km2, bare {np.sum(lc == 60) * px:.1f} km2")
    print(f"ref: {ref.sum() * px:.1f} km2; elev p5/p50/p95 = "
          f"{np.percentile(dem[ref == 1], [5, 50, 95]).round()}; "
          f"<2 km from built: {np.mean(bd[ref == 1] < 2) * 100:.0f}%")
    print(f"tundra: {tun.sum() * px:.1f} km2; slope p50 = "
          f"{np.median(slope[tun == 1]):.1f} deg")
    print(f"slope >45 deg: {np.mean(slope > 45) * 100:.2f}% of grid")
    print("\npoint            elev  slope  class        built_km  ref  tundra")
    for k, (_, _, r, c) in pts.items():
        if 0 <= r < dem.shape[0] and 0 <= c < dem.shape[1]:
            print(f"{k:16s} {dem[r, c]:5.0f} {slope[r, c]:6.1f}  "
                  f"{WC.get(int(lc[r, c]), ('?',))[0]:12s} {bd[r, c]:7.2f}  "
                  f"{ref[r, c]:3d}  {tun[r, c]:5d}")
        else:
            print(f"{k:16s} OUTSIDE GRID")


if __name__ == "__main__":
    main()
