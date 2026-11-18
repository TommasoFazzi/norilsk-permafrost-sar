# Results: Norilsk ground motion and the TPP-3 tank collapse

NASA Space Apps 2026, *Dancing with the SARs*. Final results are reproduced by
`./run_v1.sh`; the frontend package (`outputs/frontend/`, step 12) carries the
same conclusions in `story.json`. Below the final conclusions comes the full
record: v0 (the friend's first analysis, `./run_v0.sh`), the v1 robustness
checks (steps 15-17) that supersede it where they disagree, and the tank
analysis (steps 18-20).

## Final conclusions

**Question.** Did the ground under fuel tank no. 5 of TPP-3 sink because the
permafrost was degrading, and did that bring the tank down on 29 May 2020
(Zhang et al. 2022), or did the tank's pile foundation fail (the federal
inspectorate Rostekhnadzor, Nov 2020)?

**Terms.** Permafrost stays frozen all year. Above it, the active layer (the
top metre or two) thaws every summer and refreezes every winter: the surface
settles in summer and rises again in winter. *Thawing ground* below means
permafrost degradation: the active layer thawing deeper than usual into the
ice-rich top of the permafrost, which shows as extra settlement late in the
thaw season (Jul-Oct) that does not recover. All InSAR numbers are relative
to the surroundings of each target.

1. **Seasonal thaw of the active layer is real, local and steady.**
   Herbaceous wetlands sink a few mm more than their ~2 km surroundings every
   summer, following √TDD (Stefan) with the same coefficient in all six
   summers 2017-2025 (−0.195 ± 0.009 mm/√(°C·day), constant p = 0.07, no
   acceleration); the pattern repeats between any two summers (r 0.15-0.43,
   z ≥ 6.3) and NISAR L-band sees it in 2026. Lowland motion relative to
   bedrock is not measurable (atmosphere).
2. **Industry is not sinking more than other ground.** Industrial zones are
   indistinguishable from random lowland (z 0.24 / 0.07); TPP-3 never sinks
   more than its surroundings in the summers 2017-2026.
3. **No extra settlement of tank 5 in the thaw seasons.** In the snow-free
   summers 2017-2019 and the last thaw season (Jun-Nov 2019), when the active
   layer thaws, the tank moved like comparable targets (|z| ≤ 0.7; detectable:
   ~1 cm per season); Zhang's steady −35 mm/yr is not seen relative to the
   surroundings. The GNSS station NRIL shows that the bedrock under Norilsk
   did not sink (+3.0 ± 0.6 mm/yr 2017-2020); it says nothing about the
   lowland ground above the bedrock, and a uniform settlement of the whole
   area around the tank would be invisible to these relative measures.
4. **In its last winter the tank alone settled, by a small amount.** From
   6 Dec 2019 to 4 Apr 2020, with the active layer frozen, tank 5 settled
   14.4 mm relative to its surroundings, step by step from January. Only 1 of
   121 comparable targets (bright, same winter coherence, 1.5-5.5 km away)
   deviated as much (z −2.75, p ≈ 0.01); none of 20 other tanks and tank-row
   structures did, including the one 117 m away. It passes the checks we could
   run on the existing data (no unwrapping errors; median −13.8,
   coherence-weighted −15.1, height-error-corrected −13.9 mm; any single pair
   left out −7.8 to −20.5 mm) and does not follow blowing-snow episodes. But
   it rests on a few low-coherence winter pairs, single tanks elsewhere show
   ~1 cm excursions in other periods, and a snow drift cannot be fully
   excluded: a hint, not proof.
5. **Reading: this fits a failing foundation better than thawing ground.**
   Before the collapse the radar shows no extra thaw-season settlement at the
   tank or the plant, which thawing ground (as defined above) would produce;
   it shows a small settlement of the tank itself in its final months, while
   the active layer was frozen. This is consistent with the inspection (pile
   cap and piles too weak, 33 piles failed in a chain, ~30 % not on rock,
   tank returned to service in 2019 without notifying the regulator, no
   permafrost degradation found), which remains the main evidence for the
   cause.

**What remains open.** Warming (not thawing) of the frozen ground around the
piles cannot be excluded: it weakens the piles without any thaw settlement,
and at 6-12 m depth the yearly temperature peak arrives in winter. Also open:
a snow drift at the tank (point 4); tilt within the tank, below the 20 m
resolution; a uniform settlement of the whole area, invisible to relative
measures; one satellite track, with no other radar archive for 2019-20 (ASF
has no ALOS-2 there).

---

# Full record

> **v1 in one paragraph.** The local wetland subsidence passes every test
> we could run on the existing data (processing variants, spatial-shift null
> z = −6 to −9, phase closure, short- vs long-pair networks, DEM error, soil
> moisture), follows √TDD with the **same coefficient in both summers**
> (−0.205 ± 0.011 and −0.208 ± 0.045 mm/√(°C·day)), and is seen
> **independently by NISAR L-band in 2026** (z = −2.7 and −4.9). The industry
> results do not survive honest error bars: industrial ground is
> indistinguishable from random lowland ground, and the number of
> "persistent" OSM sites is at chance level (p = 0.17). Lowland motion
> relative to bedrock remains unmeasurable. A processing bug (unwrapping
> components, below) affected v0 zone-scale numbers, not the local result.
>
> **TPP-3 tank no. 5 (step 18).** No anomalous ground subsidence at the tank
> or the plant in the thaw seasons 2017-2019, so a deepening thaw of the
> ground under the tank (Zhang et al. 2022) is not supported. In the last
> winter before the collapse (Dec 2019 - Apr 2020) the tank settled ~1.5 cm
> relative to the plant (z = −2.75, p ≈ 0.01), progressively, while the
> active layer was frozen: consistent with a failing foundation (Rostekhnadzor),
> but a hint, not proof (few noisy pairs, possible snow drift; step 19 rules
> out unwrapping errors; see "Weak points"). Warming of the frozen ground around
> the piles and uniform settlement of the whole area cannot be ruled out.

## One-paragraph answer

**Thaw subsidence is real, reproducible, and controlled by the ground, not
by the industry.** Comparing every 30 m pixel with its own ~2 km
neighbourhood, the map of seasonal motion in summer 2020 (Sentinel-1A)
correlates with the one in 2025 (Sentinel-1C) at r = 0.36, against
0.007 ± 0.023 for chance (2025 map shifted by 6 km), i.e. ~15σ. The
atmosphere is independent between summers, so the repeating part is the ground.
Wetland ground (WorldCover 90, 48 km²) sinks 4.7 mm (2020) and 3.7 mm (2025)
more than its surroundings, following √TDD at every date (r = −0.99 and
−0.96): this is how a deepening thaw front behaves (Stefan's law). Industrial
sites do **not** sink more than their surroundings: the big plants
(Nadezhda smelter, TPP-3, ventilation shafts) sink *less*, in both summers.
At zone scale, industry is not detectably different from tundra (opposite
signs in the two summers). Whether the whole lowland is sinking relative to
bedrock is still not measurable: at 10-50 km scale the atmosphere dominates.

## What we measured

| | |
|---|---|
| Sentinel-1 | path 136 descending, 6 IW bursts, HyP3 ISCE2 multi-burst, 40 m, VV |
| Summer 2020 | 8 dates (15 Jun to 19 Sep), 13 pairs after QC, SBAS residual 2.0 mm |
| Summer 2025 | 7 dates (13 Jun to 5 Sep), 11 pairs after QC, SBAS residual 1.9 mm |
| Troposphere | per pair, phase vs height (+ planar ramp) fitted on coherent pixels at 250-1100 m only, subtracted everywhere |
| Coverage | 55 % (2020) and 46 % (2025) of the 30 m grid solved; incidence ~41° |
| NISAR 2026 | 4 GUNW pairs (12 days): asc 111 21 Jun–3 Jul, 20 Aug–1 Sep; desc 6 26 Jun–8 Jul, 25 Aug–6 Sep |
| Climate | ERA5-Land daily air temperature → cumulative thawing degree-days (TDD) |
| Zones | industrial (OSM, 69 km²), TPP-3 site (r = 1 km), city (WorldCover), tundra NW (64 km²), tundra SE (83 km²), plateau bedrock reference (3 patches, 3.9 km²) |

Model tested per zone and summer: `disp(t) = a + b·√TDD(t)` (Stefan-type thaw).

### Local, reproducible signal (step 14, the strongest result)

Per pixel, the seasonal motion is the √TDD fit over the summer; the local map
is that minus its 2 km Gaussian average.

| test | value |
|---|---|
| r(local 2020, local 2025), 23,095 samples 150 m apart | **0.363** |
| same with the 2025 map shifted 6 km (24 directions) | 0.007 ± 0.023 (max 0.052) |
| r of the raw (not high-passed) seasonal maps | 0.013 (atmosphere) |
| r of local map with local relief (DEM error check) | −0.09 / +0.14; partial r after removing relief 0.38 |
| pixels < −10 mm in both summers | 21,808 (19.6 km²) vs 6,699 expected if unrelated (3.3×) |

By land cover, end of summer, vs ~2 km surroundings (mm, ±1σ):

| class | 2020 | 2025 |
|---|---|---|
| herbaceous wetland (48 km²) | **−4.7 ± 0.25** | **−3.7 ± 0.24** |
| moss and lichen (213 km²) | +0.7 | +0.7 |
| grassland (115 km²) | −0.2 | +0.9 |
| tree cover | +2.7 | +2.4 |
| bare / sparse | +3.2 | +1.5 |
| built-up | +3.5 | +0.1 |

Per OSM site (inside minus a 0.3-1.5 km ring): 87 sites, site-level r between
summers 0.41. 10 sites have the same sign and > 2σ in both summers: 7 sink
*less* than their surroundings (Nadezhda plant +16.5 / +7.1 mm, TPP-3
+14.0 / +5.4, ventilation shafts VS-1,2,4 +12.8 / +9.2, Norilskgeologiya,
Kayerkan coal mine, two unnamed), 3 sink more (an unnamed site at 69.211 N
88.342 E −14.6 / −3.3, Mayak ventilation shaft, a quarry). Landfills and
tailings on average −2.8 mm in both summers. Full list:
`outputs/local_facilities.csv`; hotspot clusters: `outputs/local_hotspots.csv`.

Reading: piles, gravel pads and rock fill keep built ground stable; the
ground that gives way on thaw is the wet, ice-rich natural ground. Caveat:
C-band phase is also sensitive to soil moisture, which could add an apparent
signal of a few mm in wetlands; the monotonic √TDD shape over the whole
summer is what thaw predicts, but the two effects are not separated here.

### Paired zone comparisons (after tropospheric correction)

Slope of the per-date difference series against √TDD, in mm per √(°C·day).
Negative = the first zone sinks more per unit of thaw.

| summer | comparison | Δb | σ | z | end-of-summer Δ (mm) |
|---|---|---|---|---|---|
| 2020 | industrial − tundra mean | +0.53 | 0.23 | +2.3 | +13.7 |
| 2020 | TPP-3 site − tundra mean | +0.81 | 0.28 | +2.8 | +24.4 |
| 2020 | tundra NW − tundra SE | +1.38 | 0.76 | +1.8 | +40.3 |
| 2025 | industrial − tundra mean | −0.30 | 0.21 | −1.4 | −10.6 |
| 2025 | TPP-3 site − tundra mean | −0.19 | 0.31 | −0.6 | −6.3 |
| 2025 | tundra NW − tundra SE | −0.33 | 0.51 | −0.6 | −6.2 |

Industry sinks significantly *less* than tundra in 2020 and not
significantly more in 2025; the two tundra controls differ by up to 40 mm.
At 20-30 km separation, zone comparisons still carry turbulent atmosphere;
the local test above does not.

### Reference stability

Each plateau patch compared with the mean of the other two, maximum
deviation over the summer: 2020: 7.2 / 6.6 / 4.8 mm; 2025: 2.2 / 6.9 / 5.4 mm.
No patch moves systematically relative to the others, so all three are kept.
In late September the plateau (~650 m) decorrelates (early snow): pairs with
less than 50 % usable reference pixels were dropped (2 pairs in 2020, 4 in 2025).

### NISAR cross-check

After ionosphere and troposphere screens from the product, out-of-swath
masking and a planar deramp, the median |east-west| motion on flat ground is
7.2 mm (early) and 9.0 mm (late) per 12 days, comparable to the median
|vertical| of 8.9 and 7.5 mm. With single pairs this is at the noise level, so
**NISAR cannot yet confirm the vertical-motion assumption** used for
Sentinel-1. Track 6 descending covers only the south-west half of the AOI.

## What it supports

| Hypothesis | v0 status | Why |
|---|---|---|
| A. Climatic thaw (natural ground subsides with √TDD) | **local evidence** | wetlands sink vs surroundings along √TDD in both summers; absolute motion vs bedrock still atmosphere-dominated (σ 0.25-0.6) |
| B. Anthropogenic thaw (industry > tundra) | **not supported** | zone pairs: opposite signs; site by site, the large plants sink less than their surroundings in both summers |
| C. Stable ground | **not testable yet** | same reason as A: "not significant" here means noisy, not stable |

## What it cannot say

1. **No multi-year trend.** Snow decorrelates C-band every winter, and
   Sentinel-1 acquired nothing over Norilsk in 2022–2024. Summers are compared,
   not linked.
2. **The atmosphere dominates absolute motion.** All lowland zones move
   together by tens of mm from one pass to the next (inter-zone correlation
   0.44–1.0, mostly above 0.85); a full reversal on 19 Sep 2020 shows it is not ground. The
   reference is ~550 m above the lowland zones, so stratified water vapour maps into
   apparent lowland motion. Paired tests cancel it; absolute rates do not.
3. **Error bars are statistical only.** They come from the scatter of 7–8
   dates. Zone-local atmosphere (km scale) is not in them; small zones (TPP-3,
   3 km²) are the most affected. Our synthetic test shows the method
   recovers zone differences to ±0.05 mm/√(°C·day) for large zones and
   ±0.2 for the TPP-3 zone.
4. **Vertical is assumed** (LOS / cos 41°, slopes < 5° only). NISAR does not
   yet confirm this (above).
5. **One tank is one pixel.** At 40 m we see the ground around TPP-3, not
   the foundation of tank no. 5. The two literature coordinates for the tank
   could not be resolved with Sentinel-1 amplitude (snowmelt in spring,
   no clear disappearing target in summer); we use the OSM TPP-3 site.
   (v0. Superseded: at 20 m, step 18 locates the tank by its amplitude loss
   after the collapse, confirmed on Google Earth.)
6. **Two summers.** 2020 was exceptionally warm (MAAT −3.7 °C, TDD 1752),
   2025 colder (−7.4 °C, 1435). Two summers can show a zone difference, not a
   climate trend.

## Method checks

- End-to-end synthetic test (`./run_v0.sh --synthetic`): known per-zone slopes
  and an east motion in the industrial zone are recovered (zone differences
  within 0.05 for large zones; NISAR east 2.8 mm vs 3.0 true).
- Sign conventions checked against product documentation: HyP3 positive phase =
  motion away from the sensor; NISAR LOS unit vectors target-to-sensor (ENU);
  ISCE3 interferogram = reference × conj(secondary).
- Each HyP3 pair carries its own reference point; pairs are re-referenced to
  the plateau patches before inversion (without this the SBAS residual was
  11 mm instead of 3–6 mm; with the empirical troposphere it is 2 mm).

## v1 robustness (steps 15-17)

Numbers from `outputs/bias_checks.json`, `outputs/robustness.json`,
`outputs/era5_tropo.json`.

### A processing bug in v0: unwrapping components

SNAPHU unwraps each connected component separately. In the June pairs the
plateau reference (still snowy) is one component and the lowland another, so
lowland phase relative to the reference is off by an unknown whole number of
cycles (2.8 cm LOS, 3.7 cm vertical). v0 masked only component 0. Fraction of
each zone in the reference component:

| pairs | industrial | tundra NW | tundra SE |
|---|---|---|---|
| 2020: 06-15/06-27, 06-15/07-21, 06-27/07-21 | 0.08 | 0.00 | 0.00-0.22 |
| 2025: all pairs with 06-13, plus 06-25/07-19, 07-07/07-31, 07-31/09-05 | 0.08 | 0.00 | 0.00-0.10 |
| all other pairs | 0.99-1.00 | 0.93-1.00 | 0.97-1.00 |

Triplet closures confirm it: 59-65 % of pixels close with a whole cycle in
the triplets that contain those pairs, < 2 % elsewhere (mid-summer 0.1 %).
Fix: `09_timeseries.py --conncomp reference` (v1): 2020 → 7 dates / 10 pairs,
2025 → 6 dates / 6 pairs; reference-patch stability improves (2020 max
deviation 4.5 / 4.1 / 3.3 mm instead of 7.2 / 6.6 / 4.8). Zone-scale numbers
move by tens of mm with the processing choice (tundra NW − SE, end of 2020:
40 mm v0, 72 mm v1 height+ramp, 36 mm v1 without troposphere model): zone
comparisons 20-40 km apart are dominated by processing and atmosphere.

### Local wetland signal

"Local" = relative to the 2 km neighbourhood; contrast = herbaceous wetland
minus moss/lichen (the dominant natural class), end of summer.

| test | result | reading |
|---|---|---|
| processing variants v0 / v1 / v1 no troposphere model | contrast −5.3 / −4.3 / −4.2 mm (2020), −5.1 / −4.8 / −4.1 mm (2025); r(2020, 2025) 0.36 / 0.37 / 0.37 | robust |
| v1 with height-only troposphere | r −0.07 | unsafe model (height coefficient absorbs the ramp and imprints relief); not used |
| spatial-shift null, wetland mask moved 3-20 km (≈350 shifts) | wetland local −4.0 mm, z −6.3 (2020); −4.0 mm, z −8.9 (2025); p < 0.003 | significant with honest error bars |
| high-pass scale 1 / 2 / 4 km | r 0.34 / 0.37 / 0.39 (null ≤ 0.05) | robust |
| phase closure (fading-signal bias) | wetland local bias per triplet −0.19 to +0.14 mm | small |
| same span from short pairs only vs long pairs only | −4.9 vs −4.3 mm (2020), −3.2 vs −2.6 mm (2025); short − long −0.6 ± 0.4 mm | ≤ 15-20 % of the signal can be phase bias |
| replication with long pairs only | r 0.16 (short 0.37, null ≈ 0.01) | weaker, still far above chance |
| DEM error | would need 59 m (2020) and −15 m (2025), opposite signs; joint fit changes the wetland value by ≤ 0.1 mm | excluded |
| Stefan consistency (contrast per date vs √TDD) | b = −0.205 ± 0.011 (2020, R² 0.99), −0.208 ± 0.045 (2025, R² 0.84); difference z = 0.05 | same ground property in a warm and a cold summer |
| ERA5 0-7 cm soil moisture vs √TDD (pooled, one intercept per summer) | √TDD alone t = −9.9 (AIC −13.7); moisture alone t = 3.0 (AIC +8.7); both: √TDD t = −7.2, moisture t = 1.3 | thaw explains the series, regional moisture does not |
| NISAR L-band 2026, asc + desc vertical, 12-day pairs | contrast −1.5 mm (early window, z −2.7), −1.2 mm (late, z −4.9); S1 Stefan slope predicts −1.1 / −0.4 mm | independent sensor, wavelength, year and geometry agree in sign; late window larger than predicted |
| pixels coherent in every pair of both summers | contrast −0.3 mm (2020), −1.4 mm (2025) | the signal sits mostly in pixels that decorrelate in early-June / late-September pairs |

Open caveats. (1) The coverage result: the wettest ground decorrelates first,
so this is expected for real thaw, but it is also where unwrapping is least
reliable. (2) In flooded herbaceous wetland, C-band phase can follow the
water level (double bounce); a seasonal drawdown would also look like
subsidence growing through the summer. Regional soil moisture does not
explain the series, but water level is not measured. Ground truth (thaw
depth, water level) is needed to tell them apart.

### Six summers (2017-2021, 2025; `16_robustness.py --all-summers-only`)

All summers processed with `--conncomp reference` (the June component
problem is present every year). 2026 has only 3 usable dates and is left out.

| summer | dates | wetland − moss/lichen (mm) | Stefan b (mm/√(°C·day)) | R² | annual TDD |
|---|---|---|---|---|---|
| 2017 | 8 | −3.1 | −0.191 ± 0.039 | 0.80 | 1133 |
| 2018 | 9 | −3.2 | −0.190 ± 0.023 | 0.91 | 1279 |
| 2019 | 6 | −4.0 | −0.209 ± 0.041 | 0.87 | 1379 |
| 2020 | 7 | −4.3 | −0.205 ± 0.011 | 0.99 | 1752 |
| 2021 | 7 | −1.3 | −0.089 ± 0.035 | 0.56 | 1264 |
| 2025 | 6 | −4.8 | −0.208 ± 0.045 | 0.84 | 1435 |

- Replication: the local maps of **all 15 pairs of summers** correlate,
  r = 0.15-0.43 (median 0.36), each 6-18 σ above the shift null.
- Stefan coefficient: weighted mean −0.195 ± 0.009; constant across summers
  within the errors except 2021 (χ² = 10.2, 5 dof, p = 0.07). No trend from
  2017 to 2025 and no clear link with summer warmth (r = −0.37): the
  wetland contrast behaves as a stable property of the ground, repeated
  every summer, not as an accelerating degradation.

### Industry with honest error bars

- Zones: the industrial mask shifted over random lowland ground gives the
  error bar of the paired slope. 2020: industrial − controls +0.71, random
  lowland − controls +0.62 ± 0.37 → z = 0.24; 2025: z = 0.07. Industrial
  ground is indistinguishable from random lowland.
- OSM sites (69 with coverage): 11 "persistent" (same sign, > 2 se in both
  summers) against 6.5 expected at shifted positions (95th percentile 13),
  p = 0.17; mean inside-minus-ring +0.2 mm (null std 0.9), p = 0.75. The
  per-site list in `outputs/local_facilities.csv` is not evidence on its
  own: "the big plants sink less" is not supported as a population result.

### Absolute motion (hypotheses A and C at zone scale)

- ERA5 layer delay between lowland and plateau (surface temperature, dew
  point, pressure at the pass time): per-pair height coefficient vs the
  empirical one of step 09: r = 0.53, slope 1.0 (16 pairs). Tundra
  common-mode vs ERA5 apparent motion: r = 0.86 (2025), 0.39 (2020), but
  ERA5 explains only about half of the amplitude.
- Tundra vs bedrock slopes change sign between summers with every
  correction (e.g. tundra SE −2.0 ± 0.7 in 2020, p = 0.03; +1.1 ± 0.7 in
  2025). Still atmosphere-dominated: not testable.

### Hypotheses after v1

| Hypothesis | v1 status |
|---|---|
| A. Climatic thaw | **local evidence, strengthened**: wetlands sink with √TDD in all six summers 2017-2025 with a stable coefficient (−0.195 ± 0.009), confirmed by NISAR 2026; no acceleration over the years; zone-vs-bedrock motion still not measurable |
| B. Anthropogenic thaw | **not supported**: zone null z < 0.3; OSM sites at chance level |
| C. Stable ground | not testable at zone scale |

## TPP-3 tank no. 5: ground or foundation? (step 18)

Question: was the ground under the tank sinking before the collapse
(2020-05-29), as Zhang et al. (2022) claim (60-70 mm from Aug 2018, ~−35
mm/yr, relative to hills to the SE), or did the tank's foundation fail, as
Rostekhnadzor concluded (Nov 2020)? Damage comes from *differential*
settlement, so every number below is the tank **minus winter-coherent
targets 0.3-1.5 km around it**, compared with 200 bright targets 1.5-5.5 km
away with the same winter coherence (the null). Data: Sentinel-1 burst
136_290630_IW2, 20 m, 255 pairs 2017-06 .. 2020-10, winters included, HyP3
without water mask (with it, the tank row is masked as river).

**Where the tank is.** Among bright targets within 300 m of the published
point, three separate places lose radar amplitude from summer 2019 to summer
2020 (−6.5 to −7.7 dB, each beyond the 0.11th percentile of bright pixels;
bright pixels normally change −2.0 .. +2.8 dB). One is 40 m from the
published point (`loss_3`), the others 206 m NW and 297 m SW. The published
point itself (−4.1 dB) and `loss_3` are the same tank.

**The series cannot be continuous.** The tank decorrelates at snowmelt
(Apr-Jun) and at snow onset (Nov-Dec); nothing bridges those gaps at C-band,
so a multi-year trend is not observable. Five intervals are (vertical mm,
+ up; Zhang's rate × the interval's length for comparison):

| Interval | Zhang rate | published | loss_3 | null mean ± sd (n) | z (published) |
|---|---|---|---|---|---|
| 2017-06-19 → 10-05 | −10.3 | −13.8 | −14.9 | −3.2 ± 5.7 (171) | −1.9 ¹ |
| 2017-10-29 → 12-28 | −5.7 | −7.0 | −1.9 | −0.1 ± 5.6 (93) | −1.3 |
| 2018-06-14 → 2019-03-29 | −27.6 | **+19.1** | +23.9 | +6.8 ± 9.1 (56) | +1.4 |
| 2019-06-09 → 11-24 | −16.1 | −1.2 | −0.9 | +3.2 ± 8.9 (89) | −0.5 |
| **2019-12-06 → 2020-04-04** | −11.5 | **−14.4** | −11.6 | +0.9 ± 5.6 (121) | **−2.75 (p = 0.008)** |

¹ All three candidates, 300 m apart, drop by the same −15 mm at snow onset
(Sep 23 → Oct 5): a tank-farm-wide snow effect, not the tank. Snow-free
summers alone (published point): 2017 −5.1, 2018 −1.4, 2019 +5.3 mm, all
within the null (|z| ≤ 0.6).

The +19 mm of 2018-19 is not a rise of the tank but small snow-related
steps, each normal except the first: end of snowmelt (Jun 14 → 26, ERA5-Land
snow depth 0.43 m → 0) +8.1 (null −0.8 ± 4.3), snow-free summer −2.3
(+1.5 ± 4.7), snow onset (Sep 30 → Nov 5) +4.6 (+0.1 ± 4.3), Nov 5 → Dec 11
+3.7 (+3.2 ± 5.2, i.e. common to all structures), Dec 11 → Mar 29 +5.0
(+1.2 ± 7.1). Snow on the ground delays the radar signal of the reference
targets; steel walls standing out of the snow are less affected, so
structures look like they rise. Frost heave of the reference ground would
do the opposite.

**The last winter is different, and only at the tank.**
- Progressive, not one bad date: Dec 6 0, Jan 11 −5.2, Jan 23 −8.1,
  Feb 28 −6.8, Mar 11 −13.9, Apr 4 −14.4 mm.
- Same season one year earlier (2018-12-11 → 2019-03-29): +5.0 mm (null
  +1.2 ± 7.1, z +0.5). Not a seasonal effect.
- Localised: of 290 bright targets within 400 m, median −0.8 mm (10-90 %:
  −7.3 / +4.4); the five most negative (−12.3 to −21.3 mm) are all within
  30 m of the published point (overlapping 3x3 windows: one place, not five
  independent confirmations), four of them among the pixels that lose
  amplitude after the collapse (−3.0 to −6.0 dB); next come scattered
  targets 130-330 m away at −10.7 to −11.9 mm. The other two
  amplitude-loss places do not move (−1.0 and −4.1 mm).
- Other bright structures of the same row (probably tanks; OSM does not map
  them), Dec 6 → Apr 4: 117 m SE stays within −2.5 / +1.8 mm, 272 m SE
  within −1.5 / +5.3, 206 m NW within −5.3 / −0.1. But the 272 m SE one
  dropped −11.3 mm in winter 2018-19 (z ≈ −1.8) and did not collapse:
  single structures do move ~1 cm in a winter. Tank 5's −14.4 mm is in the
  1 % tail, not unique in kind. (Per-date table and row structures: ad-hoc
  check, not yet in step 18.)
- Robust to the reference: rings 0.3-1.5, 0.5-2.5 and 1-3 km with
  stable-target thresholds 0.5-0.7 (and 0.2-0.8 km at 0.5; too few targets
  above that) give −13.5 to −16.0 mm (published), −10.5 to −13.0 (loss_3).
- No coherence precursor: relative winter coherence 2019-20 is average (z 0.0).
- TPP-3 as a whole never sinks more than its surroundings (40 m summer
  stacks, 2017-2026: +0.9 to +12.0 mm, part B).

**Against Zhang et al.** No steady sinking: over 2018-19 the relative series
went up 19 mm (snow steps, above) where their rate predicts −28, and in the last
thaw season (Jun-Nov 2019) it moved −1 mm where it predicts −16. Their −35
mm/yr at the site is also close to their own area-wide mean (−30.8 mm/yr),
i.e. even in their data the site is not special relative to the area; at
C-band over winter tundra, a regional trend relative to distant hills is
dominated by snow and atmosphere (our steps 15-17).

**Interpretation.**
- *Ground thaw* (a deeper active layer thawing the ice-rich top of the
  permafrost): not supported. That would show as settlement in Jul-Oct, when
  the thaw depth is greatest; there is none at the tank or the plant in the
  thaw seasons 2017-2019 (Jun-Nov 2019: −1 mm). The one anomalous movement
  happened in Dec-Apr, when the active layer is frozen, and was confined to
  one tank, not the ground of the tank farm.
- *Foundation of tank 5*: consistent with the data. A ~1.5 cm (up to ~2 cm in the
  strongest pixel) settlement of that tank relative to the plant, building up
  over the 4-5 months before the collapse, matches Rostekhnadzor's finding
  (insufficient bearing capacity of the pile cap and piles, chain failure of
  33 piles, ~1.5 m settlement at failure, ~30 % of piles not on rock, tank
  returned to service in 2019 without notifying the regulator, no permafrost
  degradation found). Note: this is settlement of the whole tank relative to
  its surroundings; tank damage thresholds (API 653, GOST) refer to
  differences *within* one tank (tilt, edge settlement), which 20 m pixels
  cannot measure, so the 1-2 cm cannot be compared with them.
- Not excluded: permafrost *warming* without thawing, which weakens the
  frozen ground around the piles and lets it creep (at 6-12 m the seasonal
  temperature peak arrives months late, in winter, so it would act exactly
  then; Rostekhnadzor found no degradation, and found defective piles);
  tilt or wall deformation instead of vertical settlement (the
  radar sees any motion away from it); load changes (filling history
  unknown).
- Limits: one track, 20 m pixels (the tank is ~45 m across: no tilt), five
  intervals × four candidate points tested (Bonferroni over intervals for
  the published point: p ≈ 0.04; p = 0.008 is 1 of 121 null targets, the
  resolution limit), no independent sensor for 2019-20.
- Weak points of the winter signal (assumption review; tested in step 19,
  outcomes in the next section):
  (1) it rests on a few low-coherence pairs: Mar 11 on ~4 pairs with tank
  coherence 0.30-0.45, Apr 4 on two pairs that disagree by ~13 mm;
  unweighted SBAS. (2) A 3x3 mean with a one-cycle unwrapping error (28 mm
  LOS) in 3 of 9 pixels gives about −12 mm vertical, the size of the signal;
  needs triplet closure and per-pixel checks. (3) Dry snow delays the signal
  by about 1 mm LOS per mm of snow water equivalent; a wind drift at the base
  of an 18 m cylinder would look like a local, progressive winter
  settlement, and the null of generic bright targets does not represent
  tall tanks (the 272 m SE structure: −11.3 mm in 2018-19). (4) The tank's
  position rests on Zhang's coordinate and the amplitude loss; the 206 m NW
  place also vanished and does not move. (5) Uniform settlement of the whole
  area (Zhang's −30.8 mm/yr) is invisible to a relative measure.
- What holds: no anomalous settlement of the tank or the plant relative to
  their surroundings in the thaw seasons 2017-2019 (detectable: ~1 cm per
  season), so Zhang's "tank sinking ~35 mm/yr" is not supported locally.
  The foundation verdict rests mainly on Rostekhnadzor's inspection; the
  InSAR is consistent with it, and the winter signal is a hint, not proof.

### Checks on the winter signal (steps 19-20)

Step 19 reruns step 18 on a wider crop of the same 255 interferograms
(`norilsk_tank_wide`, ±8 km, holding 19 OSM storage tanks of the Norilsk
depots 9-13 km east); tank 5's numbers are unchanged.

| Weak point | Test | Outcome |
|---|---|---|
| Unwrapping error in the 3x3 mean | triplet closures, winter 2019-20 | none over half a cycle (tank window 0/12, its pixels 0/43; null windows 0.2 %); no pixel off by a cycle |
| Few noisy pairs | median window; coherence-weighted SBAS; leave one pair out | Dec 6 → Apr 4: −13.8 / −15.1 mm; any single pair out: −7.8 to −20.5 (→ Mar 11: −12.8 to −15.3) |
| DEM (height) error | linear trend + per-date baseline term | Δh 9.8 m, corrected −13.9 mm |
| Null not tank-like | 19 depot tanks + tank-row structures | same winter: depot tanks −2.9 ± 4.8 mm (n 19), neighbour 117 m SE +1.6; 0 of 20 as low as tank 5. The wider crop also gives a second random draw of 200 comparable targets (same rule as step 18): +2.1 ± 3.9 mm (n 165), none as low, z −4.2. The z depends on the draw, so the headline keeps the more conservative step-18 value (z −2.75, 1 of 121) |
| Snow drift | ERA5 snowfall and blowing-snow hours between acquisitions | 2019-20 was the snowiest (196 cm) and windiest (1096 h, from S) of the three winters, but the tank's steps do not follow blowing snow (r −0.09; 2018-19: r 0.56 with no drop) and the neighbour in the same row did not move |
| Position | Google Earth (user) + amplitude loss | confirmed |
| Bedrock motion | GNSS NRIL on bedrock (step 20) | +1.65 mm/yr 2000-2023, +3.0 ± 0.6 mm/yr 2017-2020, seasonal ±8 mm: the bedrock is not sinking; says nothing about the lowland ground above it |

Not exclusive to tank 5 in other periods: at snow onset 2017 the whole TPP-3
row dropped 11-20 mm (tank 5 −13.7), and the depot tanks moved −12.7 ± 4.1 mm
in Jun-Nov 2019 while tank 5 moved −1.2. Single tanks can show ~1 cm
excursions; tank 5's is the only one in that winter, at the end of its life.

## Next steps (in order of impact)

1. ~~Other Sentinel-1 summers~~ done (six summers above); 2026 needs more
   dates.
2. **3-D weather-model troposphere** (GACOS, or ERA5 pressure levels via
   PyAPS with a Copernicus CDS account) for the absolute test; the surface
   ERA5 layer model of step 17 explains only part of it.
3. **Ground truth for the wetland signal**: thaw depth vs water level
   (CALM / field data, Sentinel-2 water extent through the summer).
4. Stack all NISAR 2026 pairs per track; repeat the wetland contrast and the
   east/vertical test.
5. Annual late-summer pairs on built-up and rock to link summers.
6. Tank 5: full-resolution single-reference Sentinel-1 stack (ISCE2 topsStack
   + MiaplPy/StaMPS) to resolve wall vs roof and remove the pixel hopping;
   X-band archives (COSMO-SkyMed, TerraSAR-X) or ALOS-2 for winter 2019-20 via
   ASI/DLR/JAXA catalogues; the full Rostekhnadzor report and the tank's
   foundation levelling records, the direct test of the winter settlement.
