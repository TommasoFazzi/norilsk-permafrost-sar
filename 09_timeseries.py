"""Step 09: Sentinel-1 burst interferograms -> vertical displacement time series per summer.

Input: unzipped HyP3 INSAR_ISCE_MULTI_BURST products in data/hyp3/norilsk_p136/
(each folder: *_unw_phase.tif, *_corr.tif, *_lv_theta.tif).
Small least-squares SBAS in numpy (MintPy not used for v0):

  1. every pair is reprojected onto the 30 m analysis grid (step 05)
  2. each pair is re-referenced to the mean of the plateau patches (HyP3
     pairs each have their own arbitrary reference point); pairs where the
     reference is not measurable (< REF_COV_MIN usable, e.g. snow on the
     plateau in late September) are dropped, with any date left isolated.
  3. per summer, per pixel: solve A x = b, x = cumulative LOS displacement at
     each date relative to the first date. A pair is usable at a pixel if
     finite and coherence >= COH_MIN; pixels are solved on their own usable
     sub-network if it still connects all dates.
     The reference is re-applied per date (mean of the patches).
  1b. (--tropo, default height_ramp) empirical stratified troposphere: per
     pair, fit  b = c + k*h (+ planar ramp)  on coherent pixels between
     TROPO_HMIN and TROPO_HMAX m (plateau flanks and hills: bedrock, little
     thaw), subtract k*h (+ ramp) everywhere. The lowland (50-150 m) is NOT
     used in the fit, so lowland deformation cannot leak into k; the cost is
     a short extrapolation in height (250 -> 80 m).
     Each patch is also compared with the others (relative stability).
  4. LOS -> vertical: d_v = d_los / cos(incidence), only where slope < 5 deg;
     steeper pixels are kept in the LOS stack but NaN in the vertical one.

Sign convention (as MintPy): d_los = -phase * wavelength / (4 pi), positive =
motion TOWARD the satellite. Vertical: positive = uplift, negative = subsidence.
Incidence = 90 deg - lv_theta (lv_theta = look-vector elevation angle, rad).

Outputs (data/ts/):
  s1_<year>_vert.tif   one band per date, mm, band description = date
  s1_<year>_los.tif    same, LOS mm (no slope mask)
  s1_<year>_meta.json  dates, pairs, residual RMS, reference-patch stability

    python 09_timeseries.py [--years 2020 2025]
"""
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import Resampling

from common import grid, pair_dates, read, save_json, to_grid, write_stack, zones
from config import DIRS, S1_WAVELENGTH

COH_MIN = 0.3
REF_COV_MIN = 0.5   # min fraction of reference pixels usable in a pair
SLOPE_MAX = 5.0
PRODUCTS = DIRS["hyp3"] / "norilsk_p136"
TROPO_HMIN, TROPO_HMAX = 250.0, 1100.0
TIED_MIN = 0.5      # conncomp 'reference': min fraction of coherent lowland in the reference component


def find_pairs():
    """{year: [(ref, sec, folder)]} from the unzipped HyP3 product folders."""
    out = defaultdict(list)
    for unw in sorted(PRODUCTS.rglob("*_unw_phase.tif")):
        ref, sec = pair_dates(unw.name)
        out[int(ref[:4])].append((ref, sec, unw))
    return out


def load_pair(unw, prof, ref_flat=None):
    """unw_phase, corr, lv_theta of one pair on the analysis grid.
    ref_flat (conncomp 'reference'): also mask every pixel whose SNAPHU
    connected component is not the one of the plateau reference patches.
    Phase in another component is unwrapped up to an unknown whole number
    of cycles (2.8 cm LOS) relative to the reference."""
    stem = str(unw).removesuffix("_unw_phase.tif")
    layers, ccg = {}, None
    for key in ("unw_phase", "corr", "lv_theta"):
        with rasterio.open(f"{stem}_{key}.tif") as s:
            a = s.read(1).astype("float32")
            if s.nodata is not None:
                a[a == s.nodata] = np.nan
            if key == "unw_phase":
                a[a == 0] = np.nan          # HyP3 uses 0 outside valid data
            if key == "unw_phase" and Path(f"{stem}_conncomp.tif").exists():
                with rasterio.open(f"{stem}_conncomp.tif") as cc:
                    ccn = cc.read(1)
                    a[ccn == 0] = np.nan     # not unwrapped
                    if ref_flat is not None:
                        ccg = to_grid(ccn.astype("float32"), cc.transform, cc.crs, prof,
                                      Resampling.nearest)
            layers[key] = to_grid(a, s.transform, s.crs, prof)
    if ccg is not None:
        c = ccg.ravel()
        refc = c[ref_flat]
        refc = refc[np.isfinite(refc) & (refc > 0)]
        label = np.bincount(refc.astype(int)).argmax() if refc.size else -1
        layers["unw_phase"] = np.where((c == label).reshape(ccg.shape), layers["unw_phase"], np.nan)
    return layers


def tropo_fit(b, ok, dem, xx, yy, mode):
    """Per-pair empirical troposphere (see docstring). Returns the model to
    subtract (n_pairs, n_pix) and the height coefficients (mm/km)."""
    fit_px = (dem >= TROPO_HMIN) & (dem <= TROPO_HMAX)
    model, ks = np.zeros_like(b), []
    for k in range(b.shape[0]):
        m = ok[k] & fit_px
        cols = [np.ones(m.sum()), dem[m]] + ([xx[m], yy[m]] if mode == "height_ramp" else [])
        coef, *_ = np.linalg.lstsq(np.c_[tuple(cols)], b[k, m], rcond=None)
        model[k] = coef[1] * dem + ((coef[2] * xx + coef[3] * yy) if mode == "height_ramp" else 0)
        ks.append(round(float(coef[1]) * 1000, 1))
    return model, ks


def design(dates, pairs):
    idx = {d: i for i, d in enumerate(dates)}
    A = np.zeros((len(pairs), len(dates) - 1))
    for k, (r, s) in enumerate(pairs):
        if idx[r] > 0:
            A[k, idx[r] - 1] = -1
        A[k, idx[s] - 1] = 1
    return A


def date_components(pairs):
    """Connected components of the pair graph, as sorted date lists."""
    parent = {}
    find = lambda d: d if parent.setdefault(d, d) == d else find(parent[d])
    for r, s in pairs:
        parent[find(r)] = find(s)
    groups = defaultdict(set)
    for d in list(parent):
        groups[find(d)].add(d)
    return [sorted(g) for g in groups.values()]


def connected_dates(pairs):
    """Dates of the largest connected component of the pair graph."""
    return max(date_components(pairs), key=len)


def invert(dates, pairs, b, ok):
    """Per-pixel least-squares SBAS with missing pairs.

    b: (n_pairs, n_pix) LOS mm; ok: same shape, True where the pair is usable
    at that pixel. Pixels are grouped by their pattern of usable pairs; each
    group is solved if its sub-network still connects all dates.
    Returns x (n_dates, n_pix) relative to the first date, rms (mm)."""
    A = design(dates, pairs)
    n = len(dates) - 1
    x = np.full((len(dates), b.shape[1]), np.nan, "float32")
    weights = (1 << np.arange(len(pairs), dtype=np.int64))[:, None]
    code = (ok.astype(np.int64) * weights).sum(axis=0)
    res2, nres = 0.0, 0
    for c in np.unique(code[code > 0]):
        sel = np.flatnonzero(code == c)
        rows = np.flatnonzero((c >> np.arange(len(pairs))) & 1)
        As = A[rows]
        if len(rows) < n or np.linalg.matrix_rank(As) < n:
            continue
        sol, *_ = np.linalg.lstsq(As, b[np.ix_(rows, sel)], rcond=None)
        x[0, sel] = 0
        x[1:, sel] = sol
        r = As @ sol - b[np.ix_(rows, sel)]
        res2 += float(np.sum(r ** 2)); nres += r.size
    return x, (np.sqrt(res2 / nres) if nres else float("nan"))


def reference_stats(stack, ref_patches):
    """Per-patch mean series, each minus the mean of the other patches (mm)."""
    ids = [i for i in np.unique(ref_patches) if i > 0]
    means = {i: np.nanmean(stack[:, ref_patches == i], axis=1) for i in ids}
    stats = {}
    for i in ids:
        others = np.nanmean([means[j] for j in ids if j != i], axis=0) if len(ids) > 1 else 0
        dev = means[i] - others
        stats[f"ref_{i}"] = {"n_px": int(np.sum(np.isfinite(stack[-1][ref_patches == i]))),
                             "dev_vs_others_mm": np.round(dev, 2).tolist(),
                             "max_abs_dev_mm": float(np.nanmax(np.abs(dev)))}
    return stats


class Context:
    """Static layers shared by every summer (also used by steps 15-17)."""

    def __init__(self):
        self.prof = grid()
        self.z, self.ref_patches, self.slope = zones()
        self.shape = (self.prof["height"], self.prof["width"])
        self.dem = read(DIRS["dem"] / "dem_utm30.tif")[0].ravel()
        self.yy, self.xx = [v.ravel() / 1000.0 for v in
                            np.mgrid[0:self.shape[0], 0:self.shape[1]].astype("float32") * 30]
        self.ref_flat = (self.ref_patches > 0).ravel()
        self.lowland = np.isin(self.z, [2, 3, 4, 5, 6]).ravel()   # analysis zones off the plateau


def load_year(items, ctx, conncomp="any"):
    """Raw pairs of one summer on the grid: b (n_pairs, n_pix) LOS mm (each
    pair with its own HyP3 reference), coherence, median incidence (deg).
    conncomp: 'any' = every unwrapped component kept (v0); 'reference' =
    only the component that contains the plateau reference (load_pair)."""
    b, coh, inc = [], [], []
    for r, s, unw in items:
        L = load_pair(unw, ctx.prof, ctx.ref_flat if conncomp == "reference" else None)
        b.append(-L["unw_phase"] * S1_WAVELENGTH / (4 * np.pi) * 1000)  # mm LOS
        coh.append(L["corr"])
        inc.append(90 - np.degrees(L["lv_theta"]))
    return (np.array(b).reshape(len(items), -1), np.array(coh).reshape(len(items), -1),
            np.nanmedian(np.array(inc), axis=0))


def preprocess(b, coh, items, ctx, tropo="height_ramp", verbose=True, conncomp="any"):
    """Steps 1b + 1: troposphere model, re-reference to the plateau, drop pairs
    without a measurable reference (and, with conncomp 'reference', pairs in
    which < TIED_MIN of the coherent lowland is connected to the reference).
    Returns b, coh, pairs, dates, keep, tropo_k."""
    tropo_k = None
    if tropo != "none":
        model, tropo_k = tropo_fit(b, np.isfinite(b) & (coh >= COH_MIN), ctx.dem,
                                   ctx.xx, ctx.yy, tropo)
        b = b - model
    else:
        b = b.copy()
    # each HyP3 pair has its own arbitrary reference point: re-reference
    # every pair to the plateau patches BEFORE inversion; drop pairs where
    # the reference is not measurable (e.g. early snow on the plateau)
    pairs_all = [(r, s) for r, s, _ in items]
    keep = []
    for k in range(len(items)):
        refv = b[k, ctx.ref_flat]
        cov = np.mean(np.isfinite(refv) & (coh[k, ctx.ref_flat] >= COH_MIN))
        if cov < REF_COV_MIN:
            if verbose:
                print(f"  drop {pairs_all[k][0]}-{pairs_all[k][1]}: reference coverage {cov:.0%}")
            continue
        if conncomp == "reference":
            cohl = ctx.lowland & (coh[k] >= COH_MIN)
            tied = np.mean(np.isfinite(b[k, cohl])) if cohl.any() else 0.0
            if tied < TIED_MIN:
                if verbose:
                    print(f"  drop {pairs_all[k][0]}-{pairs_all[k][1]}: lowland tied to "
                          f"reference {tied:.0%}")
                continue
        b[k] -= np.nanmean(refv)
        keep.append(k)
    dates = connected_dates([pairs_all[k] for k in keep])
    keep = [k for k in keep if pairs_all[k][0] in dates and pairs_all[k][1] in dates]
    return b[keep], coh[keep], [pairs_all[k] for k in keep], dates, keep, tropo_k


def solve(dates, pairs, b, coh, inc, ctx):
    """Step 2-4: per-pixel SBAS, per-date plateau reference, LOS -> vertical.
    Returns los, vert (n_dates, H, W), rms, valid (flat), n_ref."""
    ok = np.isfinite(b) & (coh >= COH_MIN)
    los, rms = invert(dates, pairs, b, ok)
    valid = np.isfinite(los[-1])
    los = los.reshape(len(dates), *ctx.shape)
    ref_mask = ctx.ref_patches > 0
    n_ref = int(np.sum(np.isfinite(los[-1][ref_mask])))
    if n_ref == 0:
        raise SystemExit("no coherent pixel in reference patches")
    los -= np.nanmean(los[:, ref_mask], axis=1)[:, None, None]
    cos_inc = np.cos(np.radians(inc)).reshape(ctx.shape)
    vert = np.where(ctx.slope < SLOPE_MAX, los / cos_inc, np.nan)
    return los, vert, rms, valid, n_ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", type=int, nargs="+", default=[2020, 2025])
    ap.add_argument("--tropo", choices=["none", "height", "height_ramp"], default="height_ramp")
    ap.add_argument("--conncomp", choices=["any", "reference"], default="any",
                    help="'reference': keep only phase connected to the plateau reference")
    a = ap.parse_args()
    ctx = Context()
    DIRS["ts"].mkdir(parents=True, exist_ok=True)
    pairs_by_year = find_pairs()
    if not pairs_by_year:
        raise SystemExit(f"No *_unw_phase.tif under {PRODUCTS}: run 02 download.")

    for year in a.years:
        items = sorted(pairs_by_year.get(year, []))
        if not items:
            print(f"{year}: no pairs, skipped")
            continue
        b, coh, inc = load_year(items, ctx, a.conncomp)
        b, coh, pairs, dates, keep, tropo_k = preprocess(b, coh, items, ctx, a.tropo,
                                                         conncomp=a.conncomp)
        los, vert, rms, valid, n_ref = solve(dates, pairs, b, coh, inc, ctx)
        stab = reference_stats(los, ctx.ref_patches)

        write_stack(DIRS["ts"] / f"s1_{year}_los.tif", list(los), ctx.prof, dates)
        write_stack(DIRS["ts"] / f"s1_{year}_vert.tif", list(vert), ctx.prof, dates)
        meta = {"year": year, "dates": dates, "pairs": pairs,
                "coh_min": COH_MIN, "slope_max_deg": SLOPE_MAX,
                "valid_px": int(valid.sum()), "ref_px": n_ref,
                "residual_rms_mm": round(rms, 2), "tropo": a.tropo, "conncomp": a.conncomp,
                "tropo_k_mm_per_km": None if tropo_k is None else [tropo_k[k] for k in keep],
                "incidence_median_deg": float(np.nanmedian(inc)),
                "reference_patch_stability": stab,
                "sign": "LOS + toward satellite; vertical + up (subsidence < 0)"}
        save_json(DIRS["ts"] / f"s1_{year}_meta.json", meta)
        print(f"{year}: {len(dates)} dates, {len(pairs)} pairs, valid px "
              f"{valid.mean() * 100:.0f}%, residual RMS {rms:.1f} mm")
        for k, v in stab.items():
            print(f"  {k}: {v['n_px']} px, max |dev vs other patches| {v['max_abs_dev_mm']:.1f} mm")


if __name__ == "__main__":
    main()
