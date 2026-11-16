"""Step 02: Sentinel-1 interferograms on demand via ASF HyP3 (free, Earthdata login).

Builds a small-baseline network inside each summer (each date paired with
the next 1-2 acquisitions), so MintPy can invert it to a time series.
Uses INSAR_ISCE_MULTI_BURST on the bursts listed by step 01 (s1_bursts.csv),
so frame-number drift between years does not matter. Only dates with the
full burst set are used. Cost comes from the public HyP3 /costs table
(10x2 looks = 40 m; 4-9 bursts = 5 credits/pair). Free quota 8,000/month.

    python 02_s1_hyp3.py plan     --path P [--years 2020 2025]
    python 02_s1_hyp3.py submit   --path P [--years ...]
    python 02_s1_hyp3.py status   --path P
    python 02_s1_hyp3.py download --path P
Start small (e.g. --years 2020 2025), check one product, then add years.

Tank time series (one burst, all seasons, 20 m, products cropped to
+-CROP km around config.POINTS[--crop-point] right after download):
    python 02_s1_hyp3.py plan --path 136 --bursts s1_bursts_tank.csv --name norilsk_tank \\
        --looks 5x1 --continuous --n-next 3 --years 2017 2018 2019 2020
    ... submit / status / download with the same options (+ --crop-km 6)
"""
import argparse
import zipfile
from itertools import islice

import pandas as pd
import rasterio
import requests
from pyproj import Transformer
from rasterio.windows import Window, from_bounds

from config import DIRS, POINTS, YEARS

LOOKS = "10x2"
COSTS_URL = "https://hyp3-api.asf.alaska.edu/costs"


def credits_per_pair(n_bursts, looks=LOOKS):
    table = requests.get(COSTS_URL, timeout=30).json()
    return table["INSAR_ISCE_MULTI_BURST"]["cost_table"][looks][str(n_bursts)]


def build_pairs(dates, n_next=2, annual=False, continuous=False):
    """dates: list of (date, granules) where granules is a tuple of burst
    names for that date. Returns list of (g1, g2, d1, d2).

    Inside each year: pair each scene with the next `n_next` scenes
    (continuous=True: one chain over the whole period, across New Year).
    annual=True also links consecutive summers (late-summer to late-summer)
    so the network is connected across years; C-band coherence over one year
    is often poor on tundra, decent on built-up/rocky ground.
    """
    df = pd.DataFrame(dates, columns=["date", "granule"])
    df["date"] = pd.to_datetime(df["date"])
    df = df.drop_duplicates("date").sort_values("date")
    pairs = []
    groups = [(None, df)] if continuous else df.groupby(df["date"].dt.year)
    for _, g in groups:
        g = g.reset_index(drop=True)
        for i in range(len(g)):
            for j in range(i + 1, min(i + 1 + n_next, len(g))):
                pairs.append((g.granule[i], g.granule[j], g.date[i], g.date[j]))
    if annual:
        years = sorted(df["date"].dt.year.unique())
        for y0, y1 in zip(years, years[1:]):
            a = df[df["date"].dt.year == y0].iloc[-1]
            target = a.date + pd.DateOffset(years=y1 - y0)
            b_all = df[df["date"].dt.year == y1]
            b = b_all.iloc[(b_all.date - target).abs().argsort().iloc[0]]
            if abs((b.date - target).days) <= 15:
                pairs.append((a.granule, b.granule, a.date, b.date))
    return pairs


def chunks(it, n):
    it = iter(it)
    while batch := list(islice(it, n)):
        yield batch


def job_name(path):
    return f"norilsk_p{path}"


def crop_product(folder, lat, lon, km):
    """Crop every GeoTIFF of an unzipped product to +-km around (lat, lon) and
    drop the quick-look files, to keep high-resolution products small."""
    for tif in sorted(folder.glob("*.tif")):
        with rasterio.open(tif) as s:
            x, y = Transformer.from_crs(4326, s.crs, always_xy=True).transform(lon, lat)
            w = from_bounds(x - km * 1000, y - km * 1000, x + km * 1000, y + km * 1000, s.transform)
            w = w.round_offsets().round_lengths().intersection(Window(0, 0, s.width, s.height))
            data = s.read(window=w)
            prof = dict(driver="GTiff", width=int(w.width), height=int(w.height), count=s.count,
                        dtype=s.dtypes[0], crs=s.crs, transform=s.window_transform(w),
                        nodata=s.nodata, compress="deflate")
        tmp = tif.with_suffix(".tmp")
        with rasterio.open(tmp, "w", **prof) as d:
            d.write(data)
        tmp.replace(tif)
    for junk in list(folder.glob("*.kmz")) + list(folder.glob("*.png")):
        junk.unlink()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "submit", "status", "download"])
    ap.add_argument("--path", type=int, required=True)
    ap.add_argument("--years", type=int, nargs="+", default=YEARS)
    ap.add_argument("--n-next", type=int, default=2)
    ap.add_argument("--annual", action="store_true")
    ap.add_argument("--continuous", action="store_true", help="pairs across years (all seasons)")
    ap.add_argument("--bursts", default="s1_bursts.csv", help="burst list from step 01 (data/s1)")
    ap.add_argument("--name", help="HyP3 job name and download folder (default norilsk_p<path>)")
    ap.add_argument("--looks", default=LOOKS, choices=["20x4", "10x2", "5x1"])
    ap.add_argument("--crop-km", type=float, default=0, help="download: crop products to +-km")
    ap.add_argument("--crop-point", default="tank5_zhang2022", help="config.POINTS key for --crop-km")
    ap.add_argument("--out", help="download folder under data/hyp3 (default: the job name), "
                                  "e.g. a second, wider crop of the same jobs")
    ap.add_argument("--no-water-mask", action="store_true",
                    help="do not let HyP3 zero the phase over its water mask (it covers the "
                         "TPP-3 tank farm next to the Daldykan river)")
    a = ap.parse_args()
    name = a.name or job_name(a.path)

    if a.cmd in ("plan", "submit"):
        b = pd.read_csv(DIRS["s1"] / a.bursts)
        b = b[b.burst_id.str.startswith(f"{a.path}_") & b.complete
              & b.date.str[:4].astype(int).isin(a.years)]
        if b.empty:
            raise SystemExit("No complete burst dates: run 01 with --path.")
        n_bursts = b.burst_id.nunique()
        by_date = (b.sort_values("burst_id").groupby("date").granule
                   .apply(tuple).reset_index().values.tolist())
        pairs = build_pairs(by_date, a.n_next, a.annual, a.continuous)
        plan = pd.DataFrame(pairs, columns=["ref", "sec", "ref_date", "sec_date"])
        plan["dt_days"] = (plan.sec_date - plan.ref_date).dt.days
        plan_file = f"pairs_{name}.csv" if a.name else f"pairs_p{a.path}.csv"
        plan[["ref_date", "sec_date", "dt_days"]].to_csv(DIRS["hyp3"] / plan_file, index=False)
        cpp = credits_per_pair(n_bursts, a.looks)
        print(plan.groupby(plan.ref_date.dt.year).size().rename("pairs").to_string())
        print(f"\n{n_bursts} bursts/pair, {a.looks} looks, {cpp} credits/pair")
        print(f"{len(plan)} pairs, {len(plan) * cpp} credits")
        if a.cmd == "plan":
            return

    import hyp3_sdk
    from auth import earthdata_credentials
    user, pwd = earthdata_credentials()
    hyp3 = hyp3_sdk.HyP3(username=user, password=pwd)

    if a.cmd == "submit":
        print(f"Credits available: {hyp3.check_credits()}")
        done = {(tuple(sorted(j.job_parameters["reference"])),
                 tuple(sorted(j.job_parameters["secondary"])))
                for j in hyp3.find_jobs(name=name)
                if j.job_type == "INSAR_ISCE_MULTI_BURST"}
        todo = [p for p in pairs
                if (tuple(sorted(p[0])), tuple(sorted(p[1]))) not in done]
        print(f"Already submitted: {len(pairs) - len(todo)}; new: {len(todo)}")
        if input("Submit? [y/N] ").lower() != "y":
            return
        prepared = [hyp3_sdk.HyP3.prepare_insar_isce_multi_burst_job(
            list(g1), list(g2), name=name, looks=a.looks,
            apply_water_mask=not a.no_water_mask) for g1, g2, *_ in todo]
        for batch in chunks(prepared, 50):
            hyp3.submit_prepared_jobs(batch)
        print("Submitted. Processing takes ~1 h per batch; run 'status'.")

    elif a.cmd == "status":
        jobs = hyp3.find_jobs(name=name)
        print(pd.Series([j.status_code for j in jobs]).value_counts().to_string())

    elif a.cmd == "download":
        out = DIRS["hyp3"] / (a.out or name)
        out.mkdir(exist_ok=True)
        jobs = hyp3.find_jobs(name=name).filter_jobs(
            succeeded=True, pending=False, running=False, failed=False,
            include_expired=False)
        have = {p.stem for p in out.iterdir()}
        new = [j for j in jobs
               if j.files[0]["filename"].removesuffix(".zip") not in have]
        print(f"{len(jobs)} succeeded, {len(new)} to download -> {out}")
        for j in new:
            for z in j.download_files(out):
                with zipfile.ZipFile(z) as zf:
                    zf.extractall(out)
                z.unlink()
                if a.crop_km:
                    crop_product(out / z.stem, *POINTS[a.crop_point], a.crop_km)


if __name__ == "__main__":
    main()
