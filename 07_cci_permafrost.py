"""Step 07: ESA CCI Permafrost v5 (ground temperature, active layer) cropped to AOI.

Annual, ~1 km, model output (CryoGrid driven by MODIS LST + ERA5), so
"medium" evidence: modelled, not measured. v5 ends in 2023.
Files are circumpolar (large): each is streamed, cropped, then deleted.

    python 07_cci_permafrost.py [--years 2015 2023]
"""
import argparse
import re
import tempfile
from pathlib import Path
from urllib.parse import urljoin

import requests

from config import BBOX, DIRS

BASE = "https://dap.ceda.ac.uk/neodc/esacci/permafrost/data/"
PRODUCTS = {
    "GTD": "ground_temperature/L4/area4/pp/v05.0/northern_hemisphere/",
    "ALT": "active_layer_thickness/L4/area4/pp/v05.0/northern_hemisphere/",
}
OUT = DIRS["cci"]


def list_nc(url, depth=2):
    r = requests.get(url, timeout=60)
    if r.status_code != 200:
        print(f"  cannot list {url} ({r.status_code})")
        return []
    links = set(re.findall(r'href="([^"?#]+)"', r.text))
    files = [urljoin(url, h) for h in links if h.endswith(".nc")]
    if not files and depth:
        subdirs = [urljoin(url, h) for h in links
                   if h.endswith("/") and urljoin(url, h).startswith(url)
                   and urljoin(url, h) != url]
        for d in subdirs:
            files += list_nc(d, depth - 1)
    return sorted(set(files))


def crop(src, dst):
    import xarray as xr
    from pyproj import CRS, Transformer
    ds = xr.open_dataset(src)
    gm = next((v for v in ds.data_vars if "grid_mapping_name" in ds[v].attrs), None)
    xname = next(c for c in ds.coords if c.lower() in ("x", "lon", "longitude"))
    yname = next(c for c in ds.coords if c.lower() in ("y", "lat", "latitude"))
    if gm:
        tr = Transformer.from_crs(4326, CRS.from_cf(ds[gm].attrs), always_xy=True)
        xs, ys = tr.transform([BBOX[0], BBOX[2], BBOX[0], BBOX[2]],
                              [BBOX[1], BBOX[1], BBOX[3], BBOX[3]])
    else:  # regular lat/lon grid
        xs, ys = [BBOX[0], BBOX[2]], [BBOX[1], BBOX[3]]
    x, y = ds[xname], ds[yname]
    sx = slice(min(xs), max(xs)) if x[0] < x[-1] else slice(max(xs), min(xs))
    sy = slice(min(ys), max(ys)) if y[0] < y[-1] else slice(max(ys), min(ys))
    sub = ds.sel({xname: sx, yname: sy})
    sub.load().to_netcdf(dst)
    print(f"  -> {dst.name} {dict(sub.sizes)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=int, nargs=2, default=[2015, 2023])
    a = ap.parse_args()
    years = range(a.years[0], a.years[1] + 1)
    for prod, rel in PRODUCTS.items():
        print(f"{prod}: listing {BASE + rel}")
        for url in list_nc(BASE + rel):
            name = url.rsplit("/", 1)[-1]
            m = re.search(r"-(\d{4})-fv", name)
            if not m or int(m.group(1)) not in years:
                continue
            dst = OUT / name.replace(".nc", ".aoi.nc")
            if dst.exists():
                continue
            print(f"  {name}")
            with tempfile.TemporaryDirectory() as td:
                tmp = Path(td) / name
                with requests.get(url + "?download=1", stream=True, timeout=900) as r:
                    r.raise_for_status()
                    with open(tmp, "wb") as f:
                        for c in r.iter_content(1 << 22):
                            f.write(c)
                crop(tmp, dst)


if __name__ == "__main__":
    main()
