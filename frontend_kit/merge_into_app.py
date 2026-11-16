"""Merge the Norilsk SAR datasets of this package into the Dancing-with-the-SARs app.

    python merge_into_app.py --app /path/to/NASA-Space-Apps-Challenge-2026-...-main [--dry-run]

What it does (standard library only, like the app's own importer):
1. checks every app/sar/NOR_*.json against app/catalog_entries.json with the
   rules backend/sar_data.py applies at startup (count, unique ids, quality
   count, first/last date, point ownership, finite values);
2. copies the observation files to <app>/backend/data/sar/;
3. replaces the NOR_* entries of <app>/backend/data/sar/catalog.json and
   writes the identical catalogue to <app>/frontend/src/data/sarCatalog.json
   (the app compares both by content/hash).

Re-run it after any backend/import_finalized_sar.py run: that importer
rebuilds the catalogue from its four workbooks and drops other entries.
Afterwards: cd frontend && npm run export:ai-context && npm run build, then
restart FastAPI. Read INTEGRATION.md first: the Norilsk points have no mean
velocity and no forecast, and need the small frontend changes listed there.
"""
import argparse
import json
import math
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREFIX = "NOR_"


def check(meta, rows):
    ids = {p["id"] for p in meta["points"]}
    problems = []
    if len(rows) != meta["observationCount"]:
        problems.append("observationCount")
    if len({r["id"] for r in rows}) != len(rows):
        problems.append("duplicate observation ids")
    if sum(r["qualityPass"] for r in rows) != meta["qualityPassedCount"]:
        problems.append("qualityPassedCount")
    if min(r["date"] for r in rows) != meta["startDate"] or max(r["date"] for r in rows) != meta["endDate"]:
        problems.append("start/end date")
    if any(r.get("pointId") is not None and r["pointId"] not in ids for r in rows):
        problems.append("observation point not in dataset")
    if any(not math.isfinite(r["value"]) for r in rows):
        problems.append("non-finite value")
    if problems:
        raise SystemExit(f"{meta['id']}: {', '.join(problems)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--app", type=Path, required=True, help="root of the app repository")
    ap.add_argument("--dry-run", action="store_true", help="check only, write nothing")
    a = ap.parse_args()
    sar_dir = a.app / "backend" / "data" / "sar"
    front = a.app / "frontend" / "src" / "data" / "sarCatalog.json"
    if not (sar_dir / "catalog.json").exists() or not front.exists():
        raise SystemExit(f"Not the app repository: {a.app}")

    ours = json.loads((HERE / "catalog_entries.json").read_text(encoding="utf-8"))
    for meta in ours:
        if not meta["id"].startswith(PREFIX):
            raise SystemExit(f"Unexpected dataset id {meta['id']}")
        check(meta, json.loads((HERE / "sar" / f"{meta['id']}.json").read_text(encoding="utf-8")))

    catalogue = [m for m in json.loads((sar_dir / "catalog.json").read_text(encoding="utf-8"))
                 if not m["id"].startswith(PREFIX)] + ours
    content = json.dumps(catalogue, ensure_ascii=False, allow_nan=False, indent=2) + "\n"
    print(f"{len(ours)} Norilsk dataset(s) checked; catalogue will hold {len(catalogue)} datasets")
    if a.dry_run:
        return
    for meta in ours:
        shutil.copyfile(HERE / "sar" / f"{meta['id']}.json", sar_dir / f"{meta['id']}.json")
    (sar_dir / "catalog.json").write_text(content, encoding="utf-8")
    front.write_text(content, encoding="utf-8")
    print("written: backend/data/sar/catalog.json, frontend/src/data/sarCatalog.json, "
          + ", ".join(f"backend/data/sar/{m['id']}.json" for m in ours))
    print("next: cd frontend && npm run export:ai-context && npm run build; restart FastAPI")


if __name__ == "__main__":
    main()
