"""Presentation asset: Sentinel-2 L2A true colour of the 2020 TPP-3 diesel spill.

Planetary Computer (anonymous). Mosaics the two MGRS tiles of each date over
TPP-3 -> Daldykan -> Ambarnaya. Also a 2 km chip of the tank farm in the
least cloudy summer scene before (2019) and after (2020) the collapse, with
its lon/lat bounds for map overlays. Credit: contains modified Copernicus
Sentinel data 2019-2020.

    python s2_spill.py   -> outputs/web/s2_spill_<date>.png (+ titled copy),
                            outputs/web/s2_tank_<before|after>.png + s2_tank_chips.json
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

BOX = (87.75, 69.30, 88.10, 69.55)   # lon/lat
SCENES = {"20200531": ["S2B_MSIL2A_20200531T060639_R134_T45WWS_20201027T072154",
                       "S2B_MSIL2A_20200531T060639_R134_T45WWT_20200919T124641"],
          "20200608": ["S2A_MSIL2A_20200608T061641_R034_T45WWS_20200919T124640",
                       "S2A_MSIL2A_20200608T061641_R034_T45WWT_20200919T124645"]}
TPP3 = (69.3227, 87.9499)
TANK5 = (69.3272, 87.9339)          # config.POINTS["tank5_zhang2022"]
CHIP_HALF_DEG = (0.009, 0.026)      # lat, lon half-size: ~2 km x 2 km
CHIP_WINDOWS = {"before": "2019-06-20/2019-09-10", "after": "2020-06-20/2020-09-10"}
OUT = Path(__file__).parent / "outputs"
STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"


def read_chip(item, box):
    with rasterio.open(item.assets["visual"].href) as s:
        t = Transformer.from_crs(4326, s.crs, always_xy=True)
        xs, ys = t.transform([box[0], box[2], box[0], box[2]], [box[1], box[1], box[3], box[3]])
        return s.read(window=from_bounds(min(xs), min(ys), max(xs), max(ys), s.transform),
                      boundless=True, fill_value=0)


def tank_chips(catalog, max_items=15):
    """Clearest summer scene before and after the collapse, 2 km around tank
    5 (score: fraction of white pixels - cloud or plant steam - plus no-data
    within the chip), saved with its lon/lat bounds."""
    lat, lon = TANK5
    hy, hx = CHIP_HALF_DEG
    box = (lon - hx, lat - hy, lon + hx, lat + hy)
    meta = {"bounds_lonlat": box, "tank5_latlon": TANK5,
            "credit": "Contains modified Copernicus Sentinel data 2019-2020"}
    for name, window in CHIP_WINDOWS.items():
        items = sorted(catalog.search(collections=["sentinel-2-l2a"], bbox=box, datetime=window,
                                      query={"eo:cloud_cover": {"lt": 30}}).item_collection(),
                       key=lambda i: i.properties["eo:cloud_cover"])[:max_items]
        scored = []
        for item in items:
            a = read_chip(item, box)
            score = float(np.mean(a.min(0) > 170) + np.mean(a.max(0) == 0))
            scored.append((score, item, a))
        score, item, a = min(scored, key=lambda x: x[0])
        plt.imsave(OUT / "web" / f"s2_tank_{name}.png", np.clip(np.moveaxis(a, 0, -1) / 255 * 1.3, 0, 1))
        meta[name] = {"file": f"s2_tank_{name}.png", "item": item.id,
                      "date": item.datetime.strftime("%Y-%m-%d"),
                      "chip_white_or_nodata_fraction": round(score, 3)}
        print(name, item.id, round(score, 3))
    (OUT / "web" / "s2_tank_chips.json").write_text(json.dumps(meta, indent=1))


def main():
    (OUT / "web").mkdir(parents=True, exist_ok=True)
    catalog = pystac_client.Client.open(STAC, modifier=pc.sign_inplace)
    tank_chips(catalog)
    if all((OUT / "web" / f"s2_spill_{d}.png").exists() for d in SCENES):
        return
    col = catalog.get_collection("sentinel-2-l2a")
    for d, ids in SCENES.items():
        out = None
        for iid in ids:
            with rasterio.open(col.get_item(iid).assets["visual"].href) as s:
                t = Transformer.from_crs(4326, s.crs, always_xy=True)
                xs, ys = t.transform([BOX[0], BOX[2], BOX[0], BOX[2]],
                                     [BOX[1], BOX[1], BOX[3], BOX[3]])
                w = from_bounds(min(xs), min(ys), max(xs), max(ys), s.transform)
                a = s.read(window=w, boundless=True, fill_value=0)
                tr = s.window_transform(w)
            out = a if out is None else np.where(out.sum(0) > 0, out, a)
        rgb = np.clip(np.moveaxis(out, 0, -1) / 255 * 1.3, 0, 1)
        plt.imsave(OUT / "web" / f"s2_spill_{d}.png", rgb)
        c, r = ~tr * t.transform(TPP3[1], TPP3[0])
        fig, ax = plt.subplots(figsize=(9, 10))
        ax.imshow(rgb)
        ax.plot(c, r, "o", mfc="none", mec="magenta", ms=14, mew=2)
        ax.annotate("TPP-3", (c, r), xytext=(10, -10), textcoords="offset points",
                    color="magenta", fontsize=10, weight="bold")
        ax.axis("off")
        ax.set_title(f"Sentinel-2 L2A true colour {d[:4]}-{d[4:6]}-{d[6:]}, "
                     "TPP-3 -> Daldykan -> Ambarnaya\n"
                     "Contains modified Copernicus Sentinel data 2020", fontsize=9)
        fig.tight_layout()
        fig.savefig(OUT / f"s2_spill_{d}_titled.png", dpi=130)
        plt.close(fig)
        print(d, rgb.shape, f"valid {(out.sum(0) > 0).mean():.2f}")


if __name__ == "__main__":
    main()
