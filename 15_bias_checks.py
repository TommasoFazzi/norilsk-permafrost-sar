"""Step 15: is the local wetland signal of step 14 a radar artefact?

The wetland result (step 14) could come from the processing instead of the
ground. Four checks on the HyP3 pairs of each summer:

  1. Closure phase. For every triplet of consecutive dates (i, i+1, i+2):
         c = d(i,i+1) + d(i+1,i+2) - d(i,i+2)      (LOS mm)
     on the raw unwrapped phase, each pair re-referenced to the plateau.
     Ground motion and atmosphere belong to single dates and cancel
     exactly. What is left:
       - whole 2-pi cycles -> unwrapping error
       - the remainder     -> phase bias ("fading signal": soil moisture,
         vegetation, filtering), the C-band artefact that can look like a
         repeatable subsidence of wet ground.
  2. Network test. Displacement over the same span (date 0 -> date 2m)
     from short pairs only (consecutive dates) and from long pairs only
     (dates i, i+2). Their difference is the sum of the closures. A phase
     bias changes with the pair length; ground motion does not.
  3. DEM error. Per-date perpendicular baseline B_t (from the pair
     baselines); a DEM error dz adds dz * B_t / (R sin(inc)) to the series.
     Per-pixel joint fit disp = a + b*sqrt(TDD) + c*B_t/(R sin inc), and the
     DEM error that would be needed to fake the wetland slope.
  4. Selection. Wetland signal on pixels solved with the FULL pair network
     in both summers only (no pixel-dependent sub-networks).

All "local" values are relative to the ~2 km neighbourhood (step 14 high
pass), by ESA WorldCover class; se = std / sqrt(n_px / 25) as in step 14.
The pairs are processed as in v0 (conncomp 'any', height + ramp): the
closure test needs every unwrapped component, and local (2 km) values are
not affected by whole-cycle offsets between components except at their
edges (step 16 compares v0 and v1).

Outputs: outputs/bias_checks.json, outputs/fig_bias_checks.png

    python 15_bias_checks.py
"""
import importlib
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import read, save_json
from config import DIRS, OUTPUTS, S1_WAVELENGTH

ts = importlib.import_module("09_timeseries")
loc = importlib.import_module("14_local_signal")

YEARS = (2020, 2025)
WETLAND = 90
CLASSES = {c: n for c, n in loc.WC.items() if c in (20, 30, 50, 60, 90, 100)}
SIGMA_PX = loc.HP_KM * 1000 / loc.PX


def to_mm(phase_rad):
    return -phase_rad * S1_WAVELENGTH / (4 * np.pi) * 1000


def to_rad(mm):
    return -mm / 1000 * 4 * np.pi / S1_WAVELENGTH


def local(a2d):
    return a2d - loc.nangauss(a2d, SIGMA_PX)


def by_class(a2d, lc, extra=None):
    """Mean +- se of a map per WorldCover class (optionally within a mask)."""
    out = {}
    for c, name in CLASSES.items():
        m = (lc == c) & np.isfinite(a2d) & (True if extra is None else extra)
        if m.sum() >= 500:
            out[name] = {"mean": round(float(np.mean(a2d[m])), 2),
                         "se": round(loc.se_mean(a2d[m]), 2), "n_px": int(m.sum())}
    return out


def pair_baseline(unw_path):
    """Perpendicular baseline (m) and centre slant range (m) of a HyP3 pair."""
    stem = str(unw_path).removesuffix("_unw_phase.tif")
    txt = open(f"{stem}.txt").read()
    bperp = float(re.search(r"Baseline:\s*([-\d.eE]+)", txt).group(1))
    rng = float(re.search(r"Slant range center:\s*([-\d.eE]+)", txt).group(1))
    return bperp, rng


def ols_stack(Y, X):
    """Per-pixel OLS on pixels with all dates finite. Y (n_t, n_pix), X (n_t, k).
    Returns coefficients (k, n_pix), NaN where any date is missing."""
    full = np.all(np.isfinite(Y), axis=0)
    coef = np.full((X.shape[1], Y.shape[1]), np.nan, "float32")
    sol, *_ = np.linalg.lstsq(X, Y[:, full], rcond=None)
    coef[:, full] = sol
    return coef


def closure(b0, ok, pairs, dates, lc, shape):
    """Closure phase of every consecutive-date triplet."""
    idx = {p: k for k, p in enumerate(pairs)}
    rows, maps = [], []
    for i in range(len(dates) - 2):
        d0, d1, d2 = dates[i:i + 3]
        ks = [idx.get(p) for p in ((d0, d1), (d1, d2), (d0, d2))]
        if None in ks:
            continue
        m = ok[ks[0]] & ok[ks[1]] & ok[ks[2]]
        c = np.where(m, b0[ks[0]] + b0[ks[1]] - b0[ks[2]], np.nan)
        phi = to_rad(c)
        wrapped = np.angle(np.exp(1j * phi))
        cycles = np.round((phi - wrapped) / (2 * np.pi))
        bias = to_mm(wrapped).reshape(shape)
        cyc = cycles.reshape(shape)
        maps.append(bias)
        jumps = {}
        for cl, name in CLASSES.items():
            mm = (lc == cl) & np.isfinite(cyc)
            if mm.sum() >= 500:
                jumps[name] = round(float(np.mean(cyc[mm] != 0)), 4)
        rows.append({"triplet": [d0, d1, d2],
                     "valid_frac": round(float(m.mean()), 3),
                     "unwrap_error_frac": round(float(np.mean(cycles[m] != 0)), 4),
                     "unwrap_error_frac_by_class": jumps,
                     "bias_mm_by_class": by_class(bias, lc),
                     "bias_local_mm_by_class": by_class(local(bias), lc)})
    return rows, maps


def network_test(b1, ok, vert, pairs, dates, inc, ctx, lc):
    """Same span from short pairs, long pairs and the SBAS solution (vert)."""
    idx = {p: k for k, p in enumerate(pairs)}
    m = (len(dates) - 1) // 2
    short = [idx.get((dates[i], dates[i + 1])) for i in range(2 * m)]
    long_ = [idx.get((dates[2 * j], dates[2 * j + 2])) for j in range(m)]
    if None in short or None in long_:
        return None
    cos_inc = np.cos(np.radians(inc)).ravel()
    flat =(ctx.slope < ts.SLOPE_MAX).ravel()
    need = np.all(ok[short + long_], axis=0) & flat
    res = {"span": [dates[0], dates[2 * m]], "n_short": len(short), "n_long": len(long_),
           "px_all_pairs": int(need.sum())}
    vmaps = {}
    for name, ks in (("short", short), ("long", long_)):
        v = np.where(need, b1[ks].sum(axis=0) / cos_inc, np.nan).reshape(ctx.shape)
        vmaps[name] = local(v)
    sb = np.where(need.reshape(ctx.shape), vert[2 * m], np.nan)
    vmaps["sbas"] = local(sb)
    for name, v in vmaps.items():
        res[f"{name}_local_mm"] = by_class(v, lc)
    d = vmaps["short"] - vmaps["long"]
    res["short_minus_long_local_mm"] = by_class(d, lc)
    return res, vmaps


def dem_error(items, pairs, dates, inc_deg, stack, ctx, lc):
    """Per-date baselines, sensitivity bound and joint per-pixel fit on the
    vertical SBAS stack (n_dates, H, W)."""
    pb = {}
    rng = []
    for r, s, unw in items:
        bp, rr = pair_baseline(unw)
        pb[(r, s)] = bp
        rng.append(rr)
    A = ts.design(dates, pairs)
    bt = np.r_[0.0, np.linalg.lstsq(A, np.array([pb[p] for p in pairs]), rcond=None)[0]]
    inc = np.radians(inc_deg)
    sens = bt / (np.mean(rng) * np.sin(inc)) * 1000 / np.cos(inc)   # vertical mm per m DEM error
    x = loc.sqrt_tdd(dates)
    beta = np.polyfit(x, sens, 1)[0]                     # mm/sqrt(TDD) per m of DEM error
    Y = stack.reshape(len(dates), -1)
    X2 = np.c_[np.ones_like(x), x]
    X3 = np.c_[np.ones_like(x), x, sens]
    span = x[-1] - x[0]
    out = {"baselines_m": np.round(bt, 1).tolist(),
           "r_baseline_vs_sqrt_tdd": round(float(np.corrcoef(bt, x)[0, 1]), 3),
           "slope_bias_mm_per_sqrt_tdd_per_m_dem_error": round(float(beta), 4)}
    for name, X in (("without_dem_term", X2), ("with_dem_term", X3)):
        coef = ols_stack(Y, X)
        seas = (coef[1] * span).reshape(ctx.shape)
        out[f"seasonal_local_{name}"] = by_class(local(seas), lc)
        if name == "with_dem_term":
            dz = coef[2].reshape(ctx.shape)
            out["dem_error_m_by_class"] = by_class(dz, lc)
    wet = out["seasonal_local_without_dem_term"].get("herbaceous wetland", {}).get("mean")
    if wet is not None and beta != 0:
        out["dem_error_needed_to_fake_wetland_m"] = round(float(wet / span / beta), 1)
    return out


def main():
    ctx = ts.Context()
    lc = read(DIRS["landcover"] / "worldcover_utm30.tif")[0]
    pairs_by_year = ts.find_pairs()
    report, figmaps, full_net, local_seas = {}, {}, {}, {}
    for year in YEARS:
        items = sorted(pairs_by_year[year])
        b_raw, coh, inc = ts.load_year(items, ctx)
        b0, coh0, pairs, dates, keep, _ = ts.preprocess(b_raw, coh, items, ctx, "none", False)
        b1, coh1, pairs1, dates1, _, _ = ts.preprocess(b_raw, coh, items, ctx, "height_ramp", False)
        del b_raw
        assert pairs == pairs1 and dates == dates1
        ok = np.isfinite(b0) & (coh0 >= ts.COH_MIN)
        items = [items[k] for k in keep]
        r = {"dates": dates, "pairs": pairs}
        r["closure"], cmaps = closure(b0, ok, pairs, dates, lc, ctx.shape)
        r["closure_mean_bias_local_by_class"] = by_class(
            local(np.nanmean(np.array(cmaps), axis=0)), lc) if cmaps else None
        _, vert, *_ = ts.solve(dates, pairs, b1, coh1, inc, ctx)
        nt = network_test(b1, ok, vert, pairs, dates, inc, ctx, lc)
        if nt:
            r["network"], figmaps[year] = nt
        r["dem_error"] = dem_error(items, pairs, dates, float(np.nanmedian(inc)), vert, ctx, lc)
        full_net[year] = np.all(ok, axis=0).reshape(ctx.shape)
        local_seas[year] = local(loc.seasonal_map(vert, loc.sqrt_tdd(dates)))
        r["full_network_frac"] = round(float(full_net[year].mean()), 3)
        report[str(year)] = r
        print(f"{year}: {len(pairs)} pairs, {len(r['closure'])} triplets")

    # 4. selection: local seasonal maps (step 14 method) on full-network pixels only
    both = full_net[YEARS[0]] & full_net[YEARS[1]]
    sel = {}
    for year in YEARS:
        sel[str(year)] = {"all_pixels": by_class(local_seas[year], lc),
                          "full_network_both_summers": by_class(local_seas[year], lc, both)}
    a, b = local_seas[YEARS[0]], local_seas[YEARS[1]]
    sel["r_local_full_network"] = round(loc.corr(np.where(both, a, np.nan),
                                                 np.where(both, b, np.nan))[0], 3)
    report["selection"] = sel

    # network replication: are the long-pair-only local maps reproducible?
    if len(figmaps) == 2:
        report["network_replication_r"] = {
            k: round(loc.corr(figmaps[YEARS[0]][k], figmaps[YEARS[1]][k])[0], 3)
            for k in ("short", "long", "sbas")}
    save_json(OUTPUTS / "bias_checks.json", report)

    # ---- printout
    W = "herbaceous wetland"
    g = lambda d, k=W: d.get(k, {}).get("mean")
    for year in YEARS:
        r = report[str(year)]
        print(f"\n== {year} ==")
        for t in r["closure"]:
            print(f"  closure {t['triplet']}: unwrap errors {t['unwrap_error_frac']:.2%} "
                  f"(wetland {t['unwrap_error_frac_by_class'].get(W)}), bias local wetland "
                  f"{g(t['bias_local_mm_by_class'])} mm")
        if "network" in r:
            n = r["network"]
            print(f"  span {n['span']}: wetland local short {g(n['short_local_mm'])}, "
                  f"long {g(n['long_local_mm'])}, SBAS {g(n['sbas_local_mm'])} mm "
                  f"(short-long {n['short_minus_long_local_mm'].get(W)})")
        d = r["dem_error"]
        print(f"  DEM: r(B, sqrtTDD) {d['r_baseline_vs_sqrt_tdd']}, wetland local seasonal "
              f"{g(d['seasonal_local_without_dem_term'])} -> {g(d['seasonal_local_with_dem_term'])} mm "
              f"with DEM term; DEM error needed to fake it {d.get('dem_error_needed_to_fake_wetland_m')} m")
        s = report["selection"][str(year)]
        print(f"  selection: wetland local {g(s['all_pixels'])} (all) vs "
              f"{g(s['full_network_both_summers'])} (full network both summers)")
    print(f"\nreplication r, full network: {report['selection']['r_local_full_network']}")
    print(f"replication r by network: {report.get('network_replication_r')}")

    # ---- figure
    if len(figmaps) == 2:
        fig, axs = plt.subplots(2, 3, figsize=(16, 9))
        for row, year in enumerate(YEARS):
            for ax, k in zip(axs[row], ("short", "long", "sbas")):
                im = ax.imshow(figmaps[year][k], cmap="RdBu", vmin=-15, vmax=15)
                ax.set_title(f"{year} local, {k} pairs [mm]", fontsize=9)
                ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=axs, shrink=0.6)
        fig.savefig(OUTPUTS / "fig_bias_checks.png", dpi=100)
        plt.close(fig)


if __name__ == "__main__":
    main()
