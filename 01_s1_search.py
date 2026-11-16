"""Step 01: Sentinel-1 IW scenes and bursts over the AOI, summers 2017-2026.

No login needed. Outputs:
  data/s1/s1_scenes.csv          all IW SLC scenes (any frame)
  data/s1/s1_stacks_summary.csv  dates per summer per path, where the union of
                                 that day's frames covers >= 95 % of the AOI
  data/s1/s1_bursts.csv          SLC bursts (VV) on the chosen path whose
                                 footprints intersect the AOI
Frame numbers drift between years on the same path (e.g. 360 vs 358+363),
so stacks are defined by path + burst ID, not by frame.

    python 01_s1_search.py [--path P] [--min-burst-cov 0.02]

Single-burst time series over any window (e.g. winters before the TPP-3 tank
collapse); skips the summer scene search:
    python 01_s1_search.py --burst 136_290630_IW2 --window 2017-06-01 2020-10-15 \\
        --out s1_bursts_tank.csv
"""
import argparse

import asf_search as asf
import pandas as pd
from shapely import wkt
from shapely.geometry import shape
from shapely.ops import unary_union

from config import BBOX, DIRS, SUMMER, YEARS, bbox_wkt

aoi = wkt.loads(bbox_wkt(BBOX))


def cov(geom):
    return round(geom.intersection(aoi).area / aoi.area, 3)


def search_scenes():
    rows, geoms = [], {}
    for y in YEARS:
        res = asf.geo_search(
            dataset=asf.DATASET.SENTINEL1,
            processingLevel=asf.PRODUCT_TYPE.SLC,
            beamMode=asf.BEAMMODE.IW,
            intersectsWith=bbox_wkt(BBOX),
            start=f"{y}-{SUMMER[0]}",
            end=f"{y}-{SUMMER[1]}T23:59:59Z",
        )
        for p in res:
            pr = p.properties
            geoms[pr["sceneName"]] = shape(p.geometry)
            rows.append({
                "granule": pr["sceneName"], "time": pr["startTime"],
                "platform": pr["platform"], "path": pr["pathNumber"],
                "frame": pr["frameNumber"], "direction": pr["flightDirection"],
                "aoi_coverage": cov(geoms[pr["sceneName"]]),
            })
        print(f"{y}: {len(res)} scenes")
    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("No scenes found: check BBOX / dates.")
    df["time"] = pd.to_datetime(df["time"])
    df["date"] = df["time"].dt.date
    df["year"] = df["time"].dt.year
    df = df.sort_values("time")
    df.to_csv(DIRS["s1"] / "s1_scenes.csv", index=False)

    # union of all frames acquired on the same day on the same path
    day = (df.groupby(["path", "direction", "date", "year"])["granule"]
           .apply(lambda g: cov(unary_union([geoms[x] for x in g])))
           .rename("union_cov").reset_index())
    full = day[day.union_cov >= 0.95]
    summary = (full.groupby(["path", "direction", "year"]).size()
               .unstack("year", fill_value=0))
    summary["total"] = summary.sum(axis=1)
    summary = summary.sort_values("total", ascending=False)
    summary.to_csv(DIRS["s1"] / "s1_stacks_summary.csv")
    print(f"\nSaved {len(df)} scenes")
    print("\nDates per summer where the day's frames cover >=95% of the AOI:")
    print(summary.to_string())
    return summary


def search_bursts(path, min_cov):
    rows = []
    for y in YEARS:
        res = asf.geo_search(
            dataset=asf.DATASET.SLC_BURST, relativeOrbit=path,
            polarization="VV", intersectsWith=bbox_wkt(BBOX),
            start=f"{y}-{SUMMER[0]}", end=f"{y}-{SUMMER[1]}T23:59:59Z")
        for p in res:
            pr = p.properties
            rows.append({"burst_id": pr["burst"]["fullBurstID"],
                         "granule": pr["fileID"],
                         "date": pr["startTime"][:10],
                         "aoi_coverage": cov(shape(p.geometry))})
    b = pd.DataFrame(rows).drop_duplicates(["burst_id", "date"])
    keep = b.groupby("burst_id").aoi_coverage.max()
    keep = keep[keep >= min_cov].index
    b = b[b.burst_id.isin(keep)].sort_values(["date", "burst_id"])
    n_ids = len(keep)
    per_date = b.groupby("date").burst_id.nunique()
    b["complete"] = b.date.map(per_date == n_ids)
    b.to_csv(DIRS["s1"] / "s1_bursts.csv", index=False)
    print(f"\nPath {path}: {n_ids} bursts with >= {min_cov:.0%} AOI coverage: "
          f"{sorted(keep)}")
    yr = b.drop_duplicates("date").assign(year=lambda d: d.date.str[:4])
    print(yr.groupby("year").complete.agg(dates="size", complete="sum")
          .T.to_string())


def search_one_burst(full_id, start, end, out):
    """Every acquisition (any season) of one burst ID -> CSV in the format of
    s1_bursts.csv, so 02_s1_hyp3.py --bursts can use it."""
    path = int(full_id.split("_")[0])
    res = asf.search(dataset=asf.DATASET.SLC_BURST, relativeOrbit=path,
                     polarization="VV", fullBurstID=full_id,
                     start=start, end=f"{end}T23:59:59Z")
    rows = [{"burst_id": p.properties["burst"]["fullBurstID"], "granule": p.properties["fileID"],
             "date": p.properties["startTime"][:10], "aoi_coverage": cov(shape(p.geometry)),
             "complete": True} for p in res]
    b = pd.DataFrame(rows).drop_duplicates(["burst_id", "date"]).sort_values("date")
    b.to_csv(DIRS["s1"] / out, index=False)
    print(f"{full_id}: {len(b)} dates {b.date.min()} .. {b.date.max()} -> {out}")
    print(b.groupby(b.date.str[:7]).size().to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--path", type=int, help="default: path with most dates")
    ap.add_argument("--min-burst-cov", type=float, default=0.02)
    ap.add_argument("--burst", help="single full burst ID, e.g. 136_290630_IW2")
    ap.add_argument("--window", nargs=2, metavar=("START", "END"), help="with --burst")
    ap.add_argument("--out", default="s1_bursts_tank.csv", help="with --burst")
    a = ap.parse_args()
    if a.burst:
        if not a.window:
            raise SystemExit("--burst needs --window START END")
        search_one_burst(a.burst, *a.window, a.out)
        return
    summary = search_scenes()
    path = a.path or (int(summary.index[0][0]) if not summary.empty else None)
    if path is None:
        raise SystemExit("No path covers the AOI: pass --path.")
    search_bursts(path, a.min_burst_cov)
    print(f"\nNext: python 02_s1_hyp3.py plan --path {path} --years 2020 2025")


if __name__ == "__main__":
    main()
