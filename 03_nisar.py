"""Step 03: NISAR L-band products over the AOI (ASF, Earthdata login for download).

PROVISIONAL (calibrated) products exist for acquisitions from 2026-06-17,
so this gives summer 2026 only. L-band keeps coherence over tundra
better than Sentinel-1 C-band: the key cross-check for the InSAR result.

    python 03_nisar.py search                      # list GUNW + GCOV, no login
    python 03_nisar.py download --type GUNW --max 3
    python 03_nisar.py download --type GUNW --match _111_A_037_ --match 20260621T
    python 03_nisar.py subset                      # crop .h5 to AOI -> small .nc
Full frames are large (GB). 'subset --delete-full' removes the originals.
"""
import argparse

import pandas as pd

from config import BBOX, DIRS, NISAR_START, bbox_wkt

OUT = DIRS["nisar"]


def h5_bytes(b):
    """NISAR 'bytes' is {filename: {"bytes": n, "format": ...}} (h5 + browse PNGs)."""
    if isinstance(b, dict):
        return sum(v.get("bytes", 0) for k, v in b.items() if k.endswith(".h5"))
    return b or 0


def search():
    import asf_search as asf
    res = asf.geo_search(dataset="NISAR", processingLevel=["GUNW", "GCOV"],
                         dataMaturity="PROVISIONAL",
                         intersectsWith=bbox_wkt(BBOX), start=NISAR_START)
    rows = [{"file": p.properties["fileID"], "url": p.properties["url"],
             "type": p.properties.get("processingLevel"),
             "start": p.properties["startTime"],
             "path": p.properties.get("pathNumber"),
             "frame": p.properties.get("frameNumber"),
             "direction": p.properties.get("flightDirection"),
             "mb": round(h5_bytes(p.properties.get("bytes")) / 1e6)}
            for p in res]
    df = pd.DataFrame(rows).sort_values("start")
    df.to_csv(OUT / "nisar_products.csv", index=False)
    print(df.groupby(["type", "path", "direction"]).size().to_string())
    print(f"\n{len(df)} products -> {OUT / 'nisar_products.csv'}")


def download(ptype, n, match=None):
    import asf_search as asf
    from auth import asf_session
    df = pd.read_csv(OUT / "nisar_products.csv")
    df = df[df["type"].str.contains(ptype, na=False)]
    if match:  # every substring must appear in the file name
        for m in match:
            df = df[df["file"].str.contains(m, regex=False)]
    df = df.head(n)
    print(f"Downloading {len(df)} {ptype} files (~{df.mb.sum()} MB)")
    asf.download_urls(df.url.tolist(), str(OUT), session=asf_session(),
                      processes=2)


def crop_h5(path):
    """Crop every geocoded grid in a NISAR L2 .h5 to the AOI."""
    import h5py
    import numpy as np
    import xarray as xr
    from pyproj import Transformer

    out_vars, attrs = {}, {"source": path.name}
    with h5py.File(path, "r") as f:
        grids = []
        f.visit(lambda n: grids.append(n.rsplit("/", 1)[0])
                if n.endswith("xCoordinates") else None)
        for g in grids:
            grp = f[g]
            if "projection" not in grp:
                continue
            proj = grp["projection"]
            epsg = int(proj.attrs.get("epsg_code", proj[()]))
            x, y = grp["xCoordinates"][:], grp["yCoordinates"][:]
            tr = Transformer.from_crs(4326, epsg, always_xy=True)
            xs, ys = tr.transform([BBOX[0], BBOX[2], BBOX[0], BBOX[2]],
                                  [BBOX[1], BBOX[1], BBOX[3], BBOX[3]])
            ix = np.where((x >= min(xs)) & (x <= max(xs)))[0]
            iy = np.where((y >= min(ys)) & (y <= max(ys)))[0]
            if ix.size == 0 or iy.size == 0:
                continue
            sx, sy = slice(ix[0], ix[-1] + 1), slice(iy[0], iy[-1] + 1)
            tag = g.split("/grids/")[-1].replace("/", "_")
            dims = (f"y_{tag}", f"x_{tag}")
            coords = {dims[0]: y[sy], dims[1]: x[sx]}

            def grab(name, obj, grp=grp, shape=(y.size, x.size)):
                if isinstance(obj, h5py.Dataset) and obj.shape == shape:
                    key = f"{tag}_{name.replace('/', '_')}"
                    out_vars[key] = xr.DataArray(obj[sy, sx], dims=dims,
                                                 coords=coords)
            grp.visititems(grab)
            attrs[f"epsg_{tag}"] = epsg
    if not out_vars:
        print(f"  {path.name}: AOI not inside this frame")
        return None
    ds = xr.Dataset(out_vars, attrs=attrs)
    dst = path.with_suffix(".aoi.nc")
    ds.to_netcdf(dst)
    print(f"  {path.name} -> {dst.name} ({len(out_vars)} layers)")
    return dst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["search", "download", "subset"])
    ap.add_argument("--type", default="GUNW", choices=["GUNW", "GCOV"])
    ap.add_argument("--max", type=int, default=3)
    ap.add_argument("--match", action="append",
                    help="substring the file name must contain (repeatable)")
    ap.add_argument("--delete-full", action="store_true")
    a = ap.parse_args()
    if a.cmd == "search":
        search()
    elif a.cmd == "download":
        download(a.type, a.max, a.match)
    else:
        for h5 in sorted(OUT.glob("*.h5")):
            if crop_h5(h5) and a.delete_full:
                h5.unlink()


if __name__ == "__main__":
    main()
