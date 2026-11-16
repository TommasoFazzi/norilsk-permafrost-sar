#!/usr/bin/env bash
# Rerun the final (v1) analysis and rebuild the frontend package.
# Needs the downloads: 02 download for norilsk_p136, norilsk_tank_nowm and
# norilsk_tank_wide (see README), 03 download (NISAR). Network: ERA5 and
# Sentinel-2 (cached after the first run), GNSS NRIL (cached).
#   ./run_v1.sh            -> outputs/*.json|csv|png, outputs/frontend/, outputs/frontend_package.zip
set -euo pipefail
cd "$(dirname "$0")"
PY=".venv/bin/python -W ignore"
$PY 09_timeseries.py --conncomp reference --years 2017 2018 2019 2020 2021 2025 2026
$PY 10_nisar.py
$PY 11_analysis.py
$PY 14_local_signal.py
$PY 15_bias_checks.py
$PY 16_robustness.py
$PY 16_robustness.py --all-summers-only
$PY 17_era5_tropo.py
$PY 18_tank_series.py
$PY 19_tank_checks.py
$PY 20_gnss_nril.py
$PY s2_spill.py
$PY 12_export_web.py      # frontend package (last: reads all of the above)
# 13_build_site.py rebuilds the original website in docs/ (v0 layout): run it
# only if you mean to republish that site.
