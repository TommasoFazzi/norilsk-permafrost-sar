"""Step 11: zone x date table and the sqrt(TDD) test, per zone and per summer.

For each summer stack from step 09 (vertical, mm, relative to first date and
referenced to the plateau patches) and each zone:
  disp(t) = mean vertical displacement of the zone's pixels at date t
  fit     disp(t) = a + b * sqrt(TDD_cum(t))      (Stefan-type thaw model)
TDD_cum is the cumulative thawing degree days (ERA5-Land, from 1 Jan) at the
zone's climate point (config.ZONE_CLIMATE). b < 0 means subsidence that
grows with thaw; b is in mm per sqrt(degC day).

Uncertainty of b: OLS standard error (dates are few, so this is a lower
bound: spatially correlated atmospheric noise is not in it). The per-date
zone uncertainty is the spatial std / sqrt(n_eff), n_eff = n_px / PX_PER_IND
(neighbouring pixels are not independent).

Comparisons: PAIRED fit of the per-date difference series (zone_a - zone_b,
or zone - mean of the two controls) against sqrt(TDD): the common-mode
term shared by all lowland zones (troposphere relative to the high
reference) cancels. Also writes disp_rel_mm (zone minus mean of controls).

Gap model (PREDICTION, flagged): per zone, b averaged over summers,
predicted end-of-summer subsidence for 2022-2024 from ERA5 TDD;
leave-one-summer-out check of b.

Outputs: outputs/web/zone_timeseries.csv, outputs/zone_fits.csv,
outputs/zone_comparisons.csv, outputs/gap_prediction.csv,
outputs/fig_disp_vs_sqrt_tdd.png, outputs/fig_subsidence_<year>.png.

    python 11_analysis.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio

from common import ZONE_NAMES, zones
from config import DIRS, OUTPUTS, SUMMER, ZONE_CLIMATE

PX_PER_IND = 25          # ~150 m correlation length on a 30 m grid
GAP_YEARS = [2022, 2023, 2024]
COLORS = {"reference": "black", "industrial": "tab:orange", "accident": "tab:red",
          "urban": "gold", "control_NW": "tab:green", "control_SE": "tab:cyan"}


def tdd_series(point):
    df = pd.read_csv(DIRS["climate"] / f"era5land_daily_{point}.csv",
                     index_col="date", parse_dates=["date"])
    return df["tdd_cum"]


def fit(x, y):
    """OLS y = a + b x. Returns b, se_b, a, r2, n."""
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = len(x)
    if n < 3 or np.ptp(x) == 0:
        return dict(b=np.nan, se_b=np.nan, a=np.nan, r2=np.nan, n=n)
    X = np.c_[np.ones(n), x]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    s2 = res @ res / (n - 2)
    cov = s2 * np.linalg.inv(X.T @ X)
    r2 = 1 - (res @ res) / np.sum((y - y.mean()) ** 2) if np.ptp(y) > 0 else np.nan
    return dict(b=coef[1], se_b=np.sqrt(cov[1, 1]), a=coef[0], r2=r2, n=n)


def zone_table(z):
    rows = []
    for f in sorted(DIRS["ts"].glob("s1_*_vert.tif")):
        year = int(f.stem.split("_")[1])
        with rasterio.open(f) as s:
            stack = s.read().astype("float32")
            dates = [s.descriptions[i] for i in range(s.count)]
        for c, name in ZONE_NAMES.items():
            m = z == c
            tdd = tdd_series(ZONE_CLIMATE[name])
            for i, d in enumerate(dates):
                v = stack[i][m]
                v = v[np.isfinite(v)]
                if v.size < 10:
                    continue
                t = pd.Timestamp(d)
                rows.append({"zone": name, "summer": year, "date": t.date(),
                             "disp_mm": float(np.mean(v)), "disp_median_mm": float(np.median(v)),
                             "se_mm": float(np.std(v) / np.sqrt(max(v.size / PX_PER_IND, 1))),
                             "n_px": int(v.size), "tdd_cum": float(tdd.get(t, np.nan))})
    return pd.DataFrame(rows)


def main():
    z, _, _ = zones()
    (OUTPUTS / "web").mkdir(parents=True, exist_ok=True)
    tab = zone_table(z)
    if tab.empty:
        raise SystemExit("No s1_*_vert.tif in data/ts: run 09 first.")
    tab["sqrt_tdd"] = np.sqrt(tab["tdd_cum"])
    tab[["zone", "summer", "date", "disp_mm", "disp_median_mm", "se_mm", "n_px",
         "tdd_cum"]].to_csv(OUTPUTS / "web" / "zone_timeseries.csv", index=False)

    fits = []
    for (zone, yr), g in tab.groupby(["zone", "summer"]):
        r = fit(g["sqrt_tdd"].values, g["disp_mm"].values)
        end = g.sort_values("date").iloc[-1]
        fits.append({"zone": zone, "summer": yr, **r,
                     "seasonal_disp_mm": float(end["disp_mm"]),
                     "tdd_end": float(end["tdd_cum"])})
    fits = pd.DataFrame(fits)
    fits.round(4).to_csv(OUTPUTS / "zone_fits.csv", index=False)

    # Zone comparisons: PAIRED fit on the per-date difference series. The
    # lowland zones share a large common-mode term (troposphere between the
    # ~80 m lowland and the ~650 m reference), which cancels in a difference
    # but not in two independently fitted slopes.
    comps = []
    for yr, g in tab.groupby("summer"):
        piv = g.pivot_table(index="date", columns="zone", values="disp_mm")
        sq = g.pivot_table(index="date", columns="zone", values="sqrt_tdd")
        for a, b in [("industrial", "control_NW"), ("industrial", "control_SE"),
                     ("industrial", "controls"),
                     ("accident", "controls"), ("urban", "controls"),
                     ("control_NW", "control_SE")]:
            def series(z):
                if z == "controls":
                    if {"control_NW", "control_SE"} <= set(piv):
                        return (piv[["control_NW", "control_SE"]].mean(axis=1),
                                sq[["control_NW", "control_SE"]].mean(axis=1))
                    return None, None
                return (piv[z], sq[z]) if z in piv else (None, None)
            (da, xa), (db, xb) = series(a), series(b)
            if da is None or db is None:
                continue
            r = fit(((xa + xb) / 2).values, (da - db).values)
            end = (da - db).dropna()
            comps.append({"summer": yr, "zone_a": a, "zone_b": b, "diff_b": r["b"],
                          "se": r["se_b"], "z": r["b"] / r["se_b"] if r["se_b"] > 0 else np.nan,
                          "r2": r["r2"], "n": r["n"],
                          "seasonal_diff_mm": float(end.iloc[-1]) if len(end) else np.nan})
    # zone series relative to the mean of the two tundra controls (web)
    rel = []
    for (yr, d), g in tab.groupby(["summer", "date"]):
        g = g.set_index("zone")
        if {"control_NW", "control_SE"} <= set(g.index):
            c = g.loc[["control_NW", "control_SE"], "disp_mm"].mean()
            for zn in g.index:
                rel.append({"zone": zn, "summer": yr, "date": d, "disp_rel_mm": g.at[zn, "disp_mm"] - c})
    if rel:
        tab = tab.merge(pd.DataFrame(rel), on=["zone", "summer", "date"], how="left")
        tab[["zone", "summer", "date", "disp_mm", "disp_rel_mm", "disp_median_mm", "se_mm",
             "n_px", "tdd_cum"]].to_csv(OUTPUTS / "web" / "zone_timeseries.csv", index=False)
    comps = pd.DataFrame(comps)
    comps.round(4).to_csv(OUTPUTS / "zone_comparisons.csv", index=False)

    # gap model (prediction): b averaged over summers; LOSO check
    gap = []
    for zone, g in fits.groupby("zone"):
        b_mean = g["b"].mean()
        tdd = tdd_series(ZONE_CLIMATE[zone])
        for y in GAP_YEARS:
            s0, s1 = (pd.Timestamp(f"{y}-{SUMMER[0]}"), pd.Timestamp(f"{y}-{SUMMER[1]}"))
            pred = b_mean * (np.sqrt(tdd.get(s1, np.nan)) - np.sqrt(tdd.get(s0, np.nan)))
            gap.append({"zone": zone, "year": y, "predicted_seasonal_mm": pred,
                        "b_used": b_mean, "type": "PREDICTION (not measured)"})
        for _, row in g.iterrows():
            others = g[g.summer != row.summer]["b"]
            if len(others):
                gap.append({"zone": zone, "year": int(row.summer),
                            "loso_b_pred": others.mean(), "b_obs": row.b,
                            "type": "leave-one-summer-out check"})
    pd.DataFrame(gap).round(3).to_csv(OUTPUTS / "gap_prediction.csv", index=False)

    # figure: disp vs sqrt(TDD)
    years = sorted(tab.summer.unique())
    fig, axs = plt.subplots(1, len(years), figsize=(6 * len(years), 4.5), squeeze=False)
    for ax, yr in zip(axs[0], years):
        for zone, g in tab[tab.summer == yr].groupby("zone"):
            ax.errorbar(g.sqrt_tdd, g.disp_mm, yerr=g.se_mm, fmt="o", ms=4,
                        color=COLORS.get(zone), label=zone, capsize=2)
            f = fits[(fits.zone == zone) & (fits.summer == yr)].iloc[0]
            if np.isfinite(f.b):
                xx = np.linspace(g.sqrt_tdd.min(), g.sqrt_tdd.max(), 10)
                ax.plot(xx, f.a + f.b * xx, "-", color=COLORS.get(zone), lw=1)
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_title(f"Summer {yr}: S1 vertical vs sqrt(TDD)")
        ax.set_xlabel("sqrt(cumulative TDD) [sqrt(degC day)]")
        ax.set_ylabel("vertical displacement since first date [mm] (+ up)")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_disp_vs_sqrt_tdd.png", dpi=130)
    plt.close(fig)

    # figure: seasonal subsidence map per summer (last date)
    for f in sorted(DIRS["ts"].glob("s1_*_vert.tif")):
        with rasterio.open(f) as s:
            last = s.read(s.count)
            d0, d1 = s.descriptions[0], s.descriptions[-1]
        fig, ax = plt.subplots(figsize=(10, 7.5))
        lim = np.nanpercentile(np.abs(last), 98) if np.isfinite(last).any() else 1
        im = ax.imshow(last, cmap="RdBu", vmin=-lim, vmax=lim)
        ax.contour(z == 2, levels=[0.5], colors="orange", linewidths=0.5)
        ax.contour(z == 3, levels=[0.5], colors="magenta", linewidths=0.8)
        ax.contour(z == 1, levels=[0.5], colors="black", linewidths=0.8)
        ax.set_title(f"S1 vertical displacement {d0} -> {d1} [mm], + up "
                     "(slope < 5 deg only). Orange industrial, magenta TPP-3, black reference",
                     fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, shrink=0.7)
        fig.tight_layout()
        fig.savefig(OUTPUTS / f"fig_subsidence_{f.stem.split('_')[1]}.png", dpi=130)
        plt.close(fig)

    pd.set_option("display.width", 140)
    print(fits[["zone", "summer", "b", "se_b", "r2", "n", "seasonal_disp_mm"]].round(3).to_string(index=False))
    print()
    print(comps.round(3).to_string(index=False))
    (OUTPUTS / "analysis_summary.json").write_text(json.dumps(
        {"fits": fits.round(4).to_dict("records"),
         "comparisons": comps.round(4).to_dict("records")}, default=str, indent=1))


if __name__ == "__main__":
    main()
