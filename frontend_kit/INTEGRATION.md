# Norilsk results in the Dancing-with-the-SARs app

For the team serving the app (React + Cesium globe + FastAPI). The region
`norilsk` already exists with two sourced spill stories; this package adds
the measured results. What to say comes from `../story.json`; every file is
described in `../manifest.json` and `../README.md`.

## 0. Short on time? Do it in this order

| Step | Effort | What the user gets |
|---|---|---|
| 1. Paste `app/norilsk_stories.ts` into `stories.ts` | minutes, no logic | two result stories on the globe (thaw; tank 5), with sources; Earth AI learns them after `npm run export:ai-context` |
| 2. Show `images/` (before/after chips, figures) in those stories or a methods panel | small | the evidence as pictures |
| 3. Add 1-2 `overlays/` (e.g. `overlay_local_persistent.webp`) as image layers | ~20 lines Cesium | the "where the ground sinks" map |
| 4. Merge the SAR dataset (section 2 changes) | 4 small code changes | tank 5 readings in Analyze |

Steps 1-2 alone tell the whole story correctly.

## 1. What goes where

| Package file | App destination | How |
|---|---|---|
| `app/sar/NOR_TPP3_TANK_VERTICAL.json`, `app/catalog_entries.json` | `backend/data/sar/`, `catalog.json`, `frontend/src/data/sarCatalog.json` | `python app/merge_into_app.py --app <repo>` (checks like `sar_data.py`, keeps both catalogues identical) |
| `overlays/*.webp` + corners in `manifest.json` | globe overlays for region `norilsk` (`frontend/src/globe/overlays.ts`) | Cesium `SingleTileImageryProvider` with `Rectangle.fromDegrees(west, south, east, north)` from `overlay_corners_lonlat` |
| `app/norilsk_stories.ts` (from `story.json`) | `frontend/src/data/stories.ts` (`RegionStory`) and Earth AI context (`npm run export:ai-context`) | paste; category is `water` (no ground-motion category exists) |
| `tables/*.csv` | charts in Analyze (optional) | units in column names |
| `images/s2_tank_*.png` | Before/After for the tank story | `bounds_lonlat` in manifest |

Re-run `merge_into_app.py` after any `import_finalized_sar.py`: that importer
rebuilds the catalogue from its workbooks and drops other datasets. Then
`npm run export:ai-context`, `npm run build`, restart FastAPI.

## 2. The dataset is not EGMS: four small frontend changes are required

`NOR_TPP3_TANK_VERTICAL` has 3-6 points (tank 5 and the bright structures of
its row), vertical mm relative to stable targets 0.3-1.5 km away. Unlike the
Venice/London EGMS points:

1. **No mean velocity** (`velocity: null`). C-band loses the tanks at snowmelt
   and snow onset, so the series is split into segments that cannot be joined;
   a velocity would be invented. `ObservationPanel.tsx` calls
   `v.toLocaleString()` and would throw:
   ```ts
   const numeric = (v: number | null) =>
     v == null ? "—" : v.toLocaleString(undefined, { maximumFractionDigits: 2 });
   ```
   and in `sar.ts` type `velocity: number | null`.
2. **Segments.** Each reading has `details.segment`; each segment starts at 0.
   In `sar.ts` `observationPlot`, also start a new path when the segment
   changes (two segments are only 12 days apart, under the 45-day gap rule):
   ```ts
   rows[i].details.segment !== rows[i - 1].details.segment ||
   ```
   and in the "Change from" sentence of `ObservationPanel.tsx`, show the
   difference only when both readings have the same `details.segment`
   (otherwise: "Different segments: not comparable").
3. **No forecast.** The catalogue entry has `forecastAvailable: false`.
   Do not open the Forecast tab/panel for such datasets (App.tsx, the
   `mode === "forecast" && activeDataset?.kind === "subsidence"` branch);
   `/api/datasets/NOR_.../forecast` has nothing to return. Never fill it.
4. **Badge.** The panel prints "EGMS POINTS" for every subsidence dataset; use
   `dataset.pointBadge ?? "EGMS POINTS"` ("SENTINEL-1 TARGETS" here).

The backend needs no change: `Dataset` ignores the extra fields
(`forecastAvailable`, `pointBadge`), `points` is free-form.

## 3. Reading the numbers correctly (put this in the UI text)

- Negative = down. Values are tank minus surroundings, not absolute motion.
- 20 m pixels: a point is a small area, not a building part; tilt within the
  45 m tank cannot be measured.
- The key reading: winter 2019-20 (segment starting 2019-12-06), tank 5 goes
  to about -14 mm by March-April while its neighbours stay within a few mm.
  Snow-free summers show nothing unusual. It is a hint of a failing
  foundation, not proof (see `story.json` caveats).
- Maps (`overlays/`): red = down, blue = up, transparent near 0. The
  `local_*` maps are relative to ~2 km surroundings (where the reproducible
  signal is); the `s1_seasonal_*` maps are relative to plateau bedrock and
  include atmosphere at 10-50 km scale.

## 4. Sources to cite in stories

- Zhang et al. 2022, Remote Sensing 14(19), 5036, doi:10.3390/rs14195036
  (claims permafrost-driven subsidence at the tank; not confirmed here).
- Rostekhnadzor investigation, Nov 2020 (reported by Interfax,
  https://www.interfax.ru/russia/737323): weak pile foundation, no permafrost
  degradation found.
- Data: Copernicus Sentinel-1 (ASF HyP3), NISAR (NASA/ISRO, ASF DAAC), ERA5
  (Copernicus C3S via Open-Meteo), GNSS NRIL (Nevada Geodetic Laboratory),
  Sentinel-2 (Microsoft Planetary Computer), OpenStreetMap.
