"""Step 20: GNSS station NRIL (Norilsk) - does the bedrock move?

NRIL is an IGS station on the roof of the Seismology Building in Norilsk
(stainless-steel rod, 2.6 m monument; bedrock: flood basalts, basement
outcrop; IGS site log nril00rus_20250108). It therefore measures the
building on bedrock, not the permafrost ground. Daily positions from the
Nevada Geodetic Laboratory (IGS20 frame, free), 2000-09 .. 2023-01.

Used for two assumptions:
- the InSAR reference (plateau bedrock) is stable: bedrock vertical trend
  and seasonal amplitude at NRIL;
- Zhang et al. 2022 report -30.8 mm/yr over their whole study area relative
  to SE hills: a regional (crustal) subsidence of that size would show here.

Outputs: data/gnss/NRIL.tenv3 (cache), outputs/gnss_nril.json,
outputs/gnss_nril_up.csv (weekly median, mm), outputs/fig_gnss_nril.png

    python 20_gnss_nril.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

from common import save_json
from config import DATA, OUTPUTS

URL = "https://geodesy.unr.edu/gps_timeseries/IGS20/tenv3/IGS20/NRIL.tenv3"
CACHE = DATA / "gnss" / "NRIL.tenv3"
COLS = ("site ymd dec mjd week dow reflon e0 east n0 north u0 up ant sig_e sig_n sig_u "
        "c_en c_eu c_nu lat lon h").split()
WINDOW = (2017.4, 2020.5)          # Sentinel-1 tank series period
SUMMER = ("06-15", "09-30")


def load():
    if not CACHE.exists():
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(URL, timeout=120)
        r.raise_for_status()
        CACHE.write_text(r.text)
    df = pd.read_csv(CACHE, sep=r"\s+", skiprows=1, names=COLS)
    df["date"] = pd.to_datetime(df.mjd - 40587, unit="D")
    df["up_mm"] = (df.u0 + df.up) * 1000
    return df


def fit(df):
    """Trend (mm/yr) with annual + semi-annual terms; formal se inflated by
    the lag-1 autocorrelation of the residuals (daily GNSS is not white)."""
    t = df.dec.values - df.dec.mean()
    w = 2 * np.pi * df.dec.values
    X = np.c_[np.ones_like(t), t, np.cos(w), np.sin(w), np.cos(2 * w), np.sin(2 * w)]
    coef, *_ = np.linalg.lstsq(X, df.up_mm.values, rcond=None)
    res = df.up_mm.values - X @ coef
    rho = np.corrcoef(res[:-1], res[1:])[0, 1]
    n_eff = len(res) * (1 - rho) / (1 + rho)
    se = np.sqrt(res @ res / (len(res) - X.shape[1]) * np.linalg.inv(X.T @ X)[1, 1] * len(res) / n_eff)
    peak_day = int((np.degrees(np.arctan2(coef[3], coef[2])) % 360) / 360 * 365.25)
    return {"trend_mm_yr": round(float(coef[1]), 2), "trend_se": round(float(se), 2),
            "annual_amp_mm": round(float(np.hypot(coef[2], coef[3])), 1),
            "annual_peak_day_of_year": peak_day,
            "residual_sd_mm": round(float(res.std()), 1), "n_days": int(len(res))}, X @ coef


def main():
    df = load()
    rep = {"station": "NRIL", "lat": round(float(df.lat.iloc[-1]), 4),
           "lon": round(float(df.lon.iloc[-1] % 360), 4),
           "monument": "roof of Seismology Building, 2.6 m steel rod; bedrock flood basalts",
           "source": URL, "first": str(df.date.min().date()), "last": str(df.date.max().date())}
    rep["all"], _ = fit(df)
    win = df[(df.dec >= WINDOW[0]) & (df.dec <= WINDOW[1])]
    rep["2017_2020"], model = fit(win)
    summers = {}
    for y in range(2017, 2021):
        s = df[(df.date >= f"{y}-{SUMMER[0]}") & (df.date <= f"{y}-{SUMMER[1]}")]
        if len(s) > 60:
            a = s.up_mm.iloc[:10].median()
            b = s.up_mm.iloc[-10:].median()
            summers[str(y)] = round(float(b - a), 1)
    rep["summer_change_mm"] = summers
    rep["zhang2022_area_mean_mm_yr"] = -30.8
    save_json(OUTPUTS / "gnss_nril.json", rep)
    wk = (df.set_index("date").up_mm.resample("W").median().dropna())
    (wk - wk.loc["2017"].mean()).round(1).rename("up_mm_vs_2017_mean").to_csv(OUTPUTS / "gnss_nril_up.csv")
    print(rep)

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(wk.index, wk - wk.loc["2017"].mean(), ".", ms=2, color="grey", label="NRIL weekly median")
    ax.plot(win.date, model - wk.loc["2017"].mean(), color="C3", lw=1,
            label=f"fit 2017-2020: {rep['2017_2020']['trend_mm_yr']:+.1f} ± {rep['2017_2020']['trend_se']:.1f} mm/yr")
    t0 = pd.Timestamp("2018-08-25")
    ax.plot([t0, pd.Timestamp("2020-05-29")], [0, -30.8 * 1.76], "k:", label="Zhang 2022 area mean, -30.8 mm/yr")
    ax.axvline(pd.Timestamp("2020-05-29"), color="red", lw=0.8)
    ax.set_ylabel("vertical [mm], + up (vs 2017 mean)")
    ax.set_title("GNSS NRIL (Norilsk, on bedrock): the bedrock is not sinking", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUTS / "fig_gnss_nril.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
