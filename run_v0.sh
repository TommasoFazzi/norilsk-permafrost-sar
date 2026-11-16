#!/usr/bin/env bash
# Rerun the v0 analysis (needs downloaded data: 02 download, 03 download).
# NOTE: this reproduces the friend's v0 and OVERWRITES data/ts with v0 stacks
# (unwrapping components not checked). For the final results use run_v1.sh.
#   ./run_v0.sh            real data
#   ./run_v0.sh --synthetic  end-to-end test on synthetic data (no login)
set -euo pipefail
cd "$(dirname "$0")"
PY=.venv/bin/python
if [[ "${1:-}" == "--synthetic" ]]; then
  SYN="${TMPDIR:-/tmp}/norilsk_synthetic"
  rm -rf "$SYN"; mkdir -p "$SYN/data" "$SYN/outputs"
  $PY tests/make_synthetic.py "$SYN/data"
  export NORILSK_DATA="$SYN/data" NORILSK_OUT="$SYN/outputs"
fi
$PY 09_timeseries.py --years 2020 2025
$PY 10_nisar.py
$PY 11_analysis.py
$PY 14_local_signal.py
$PY 12_export_web.py
if [[ "${1:-}" == "--synthetic" ]]; then
  $PY tests/check_synthetic.py
  echo "synthetic outputs: $NORILSK_OUT"
else
  $PY 13_build_site.py      # website data (docs/), real data only
fi
