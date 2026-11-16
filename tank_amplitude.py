"""Locate the HPP-3 (TPP-3) tank collapse of 2020-05-29 with S1 amplitude.

A large steel tank is a strong, stable radar reflector; its collapse should
show as a bright target disappearing between the last pre-event and first
post-event acquisitions on the same track (path 136 DESC, same geometry).
Reads small windows from Planetary Computer `sentinel-1-rtc` COGs
(geocoded, terrain-corrected gamma0; anonymous access, no Earthdata login).

    python tank_amplitude.py [summer|spring] -> outputs/tank_amplitude_*.png
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import planetary_computer as pc
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

from config import DIRS, POINTS

EVENT = "2020-05-29"
PATH = 136
HALF = 750  # m, default window half-size
OUT = Path(__file__).parent / "outputs"
OUT.mkdir(exist_ok=True)
# name: (lat, lon, half-window m). tpp3_osm = centroid of OSM way/525598517
# "Норильская ТЭЦ-3", included because candidate A is 0.8 km from it.
CANDS = {"hpp3_tank_a": (*POINTS["hpp3_tank_a"], HALF),
         "hpp3_tank_b": (*POINTS["hpp3_tank_b"], HALF),
         "tpp3_osm": (69.3227, 87.9499, 1500)}


# Spring pre/post is dominated by snowmelt (10 May snow, 3 June wet/bare).
# Default: same snow-free season one year apart, so only real changes remain.
PERIODS = {"spring": ["2020-04-01/2020-07-10"],
           "summer": ["2019-06-20/2019-08-31", "2020-06-20/2020-08-31"]}


def items(mode):
    cat = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=pc.sign_inplace)
    lat, lon = POINTS["hpp3_tank_a"]
    out = []
    for dt in PERIODS[mode]:
        s = cat.search(collections=["sentinel-1-rtc"], datetime=dt,
                       intersects={"type": "Point", "coordinates": [lon, lat]},
                       query={"sat:relative_orbit": {"eq": PATH}})
        out += list(s.items())
    return sorted(out, key=lambda i: i.datetime)


def window(item, lat, lon, half=HALF, pol="vv"):
    with rasterio.open(item.assets[pol].href) as src:
        x, y = Transformer.from_crs(4326, src.crs, always_xy=True).transform(lon, lat)
        w = from_bounds(x - half, y - half, x + half, y + half, src.transform)
        a = src.read(1, window=w, boundless=True, fill_value=0).astype("float32")
        tr = src.window_transform(w)
        crs = src.crs
    a[a <= 0] = np.nan
    return 10 * np.log10(a), tr, crs


def osm_points(crs, tr):
    path = DIRS["osm"] / "osm_tanks_plants.geojson"
    if not path.exists():
        return []
    t = Transformer.from_crs(4326, crs, always_xy=True)
    pts = []
    for f in json.loads(path.read_text())["features"]:
        lon, lat = f["geometry"]["coordinates"]
        c, r = ~tr * t.transform(lon, lat)
        kind = "plant" if "power" in f["properties"] else "tank"
        pts.append((c, r, kind))
    return pts


def main(mode="summer"):
    its = items(mode)
    dates = [i.datetime.strftime("%Y-%m-%d") for i in its]
    pre = [k for k, d in enumerate(dates) if d < EVENT]
    post = [k for k, d in enumerate(dates) if d > EVENT]
    print("pre:", [dates[k] for k in pre], "\npost:", [dates[k] for k in post])

    for name, (lat, lon, half) in CANDS.items():
        stack, tr, crs = [], None, None
        for it in its:
            a, tr, crs = window(it, lat, lon, half)
            stack.append(a)
        stack = np.array(stack)
        n = stack.shape[1]
        # last pre vs first post, and multi-date means (speckle-robust)
        d_pair = stack[post[0]] - stack[pre[-1]]
        d_mean = np.nanmean(stack[post], 0) - np.nanmean(stack[pre], 0)
        pts = osm_points(crs, tr)

        cols = len(dates) + 2
        fig, axs = plt.subplots(1, cols, figsize=(2.6 * cols, 3.2))
        for k, (ax, d) in enumerate(zip(axs, dates)):
            ax.imshow(stack[k], cmap="gray", vmin=-25, vmax=5)
            ax.set_title(("PRE " if k in pre else "POST ") + d[5:], fontsize=8)
        for ax, arr, t in [(axs[-2], d_pair, f"{dates[post[0]][5:]} - {dates[pre[-1]][5:]}"),
                           (axs[-1], d_mean, "mean post - mean pre")]:
            ax.imshow(arr, cmap="RdBu", vmin=-10, vmax=10)
            ax.set_title(t + " [dB]", fontsize=8)
        for ax in axs:
            ax.plot(n / 2, n / 2, "+", color="magenta", ms=10)
            for c, r, kind in pts:
                if 0 <= c < n and 0 <= r < n:
                    ax.plot(c, r, "o" if kind == "tank" else "s", mfc="none",
                            mec="yellow" if kind == "tank" else "cyan", ms=4, mew=0.6)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{name} ({lat}, {lon}) S1 path {PATH} RTC VV, "
                     f"{2 * half} m window; +=candidate, o=OSM tank, []=OSM plant",
                     fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / f"tank_amplitude_{mode}_{name}.png", dpi=130)
        plt.close(fig)

        # strongest losses in the window (bright pre target that vanished)
        bright_pre = np.nanmean(stack[pre], 0) > 0  # gamma0 > 0 dB: point-like
        loss = np.where(bright_pre, d_mean, np.nan)
        px = abs(tr.a)
        print(f"\n{name}: {np.sum(bright_pre)} bright pre-event px; "
              f"px with mean loss < -6 dB: {np.nansum(loss < -6)}")
        if np.isfinite(loss).any():
            r, c = np.unravel_index(np.nanargmin(loss), loss.shape)
            x, y = tr * (c + 0.5, r + 0.5)
            lo, la = Transformer.from_crs(crs, 4326, always_xy=True).transform(x, y)
            print(f"  largest loss {loss[r, c]:.1f} dB at {la:.4f}N {lo:.4f}E "
                  f"({np.hypot(c - n / 2, r - n / 2) * px:.0f} m from candidate)")


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "summer")
