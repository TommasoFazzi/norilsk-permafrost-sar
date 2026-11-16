"""Compare what steps 09-11 recovered with the synthetic truth.

    NORILSK_OUT=... python tests/check_synthetic.py
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config import DIRS, OUTPUTS  # noqa: E402
from tests.make_synthetic import B, EAST_MM  # noqa: E402

fits = pd.read_csv(OUTPUTS / "zone_fits.csv")
fits["b_true"] = fits.zone.map(B)
fits["err"] = fits.b - fits.b_true
print(fits[["zone", "summer", "b_true", "b", "se_b", "err", "r2"]].round(3).to_string(index=False))
# Absolute b carries a common-mode error (reference-patch atmosphere is
# subtracted from every zone on each date): loose tolerance. Zone
# differences (what the hypothesis test uses) cancel it: tight tolerance.
ok = (fits.err.abs() < 0.35).all()
comp = pd.read_csv(OUTPUTS / "zone_comparisons.csv")
truth = dict(B, controls=(B["control_NW"] + B["control_SE"]) / 2)
comp["diff_true"] = comp.zone_a.map(truth) - comp.zone_b.map(truth)
comp["err"] = comp.diff_b - comp.diff_true
print(comp[["summer", "zone_a", "zone_b", "diff_true", "diff_b", "err"]].round(3).to_string(index=False))
big = ~comp.zone_a.isin(["urban", "accident"])   # small/fragmented zones blur at 40 m
ok &= (comp[big].err.abs() < 0.1).all()
meta = json.loads((DIRS["ts"] / "nisar_meta.json").read_text())
for w in ("early", "late"):
    pz = meta[w]["per_zone"]
    print(f"NISAR {w}: industrial east {pz['industrial']['east_mean_mm']} mm (true {EAST_MM}), "
          f"control_NW east {pz['control_NW']['east_mean_mm']} mm (true 0)")
    ok &= abs(pz["industrial"]["east_mean_mm"] - EAST_MM) < 1.0
    ok &= abs(pz["control_NW"]["east_mean_mm"]) < 1.0
print("SYNTHETIC TEST", "PASSED" if ok else "FAILED")
sys.exit(0 if ok else 1)
