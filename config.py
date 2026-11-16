"""Shared settings: study area, points of interest, periods, output folders."""
import os
from pathlib import Path

# ---------------------------------------------------------------- study area
# Bounding box (lon_min, lat_min, lon_max, lat_max), WGS84.
# Covers Norilsk city, Talnakh (N), Kayerkan + HPP-3 (W), Nadezhda plant,
# surrounding tundra and the Kharaelakh plateau NE of Talnakh (bedrock
# reference; BBOX extended east to 88.90 E for it).
BBOX = (87.55, 69.20, 88.90, 69.58)

# Metric CRS for derived rasters (slope etc.): UTM 45N covers 84-90 E.
UTM_EPSG = 32645

# Points of interest (lat, lon). APPROXIMATE: check every point on
# Google Earth / OSM before using it in analysis.
POINTS = {
    "norilsk_centre": (69.349, 88.201),
    # Early guesses of the tank-5 location, superseded by tank5_zhang2022.
    "hpp3_tank_a": (69.3278, 87.9644),
    "hpp3_tank_b": (69.3794, 87.7444),
    "talnakh": (69.49, 88.39),
    # TPP-3 (OSM way/525598517 centroid): centre of the "accident" zone.
    "tpp3_accident": (69.3227, 87.9499),
    # Collapse site given by Zhang et al. 2022 (Remote Sens. 14, 5036):
    # 69 19'37.74"N 87 56'1.86"E, 0.35 km from the TPP-3 polygon. Confirmed
    # on Google Earth historical imagery (2026-11-15) and by the Sentinel-1
    # amplitude loss after the collapse (step 18). (hpp3_tank_b equals the
    # Wikipedia page coordinate, 9 km away: a generic location, not the tank.)
    "tank5_zhang2022": (69.3272, 87.9339),
    # Centre of a +-8 km crop holding tank 5 and the OSM storage tanks of the
    # Norilsk depots 9-13 km east (tank-only null, step 19).
    "burst_tanks_centre": (69.34, 88.093),
    # Zone points from step 08 (zones v3): centre of the largest patch.
    "control_NW": (69.4467, 87.6366),
    "control_SE": (69.2301, 88.5044),
    # Kharaelakh plateau bedrock reference patches (step 08), east of 88.60 E,
    # outside the mine-shaft hulls + 3 km. plateau_reference = ref_1.
    "plateau_reference": (69.4955, 88.7356),
    "ref_1": (69.4955, 88.7356),
    "ref_2": (69.4849, 88.6588),
    "ref_3": (69.5122, 88.6786),
}

# ------------------------------------------------------------------ periods
YEARS = list(range(2017, 2027))       # summers 2017-2026
# Snow-free InSAR window (MM-DD). Snow decorrelates C-band; keep summer only.
SUMMER = ("06-10", "10-05")
NISAR_START = "2026-06-17"            # first PROVISIONAL (calibrated) products

# ------------------------------------------------------------------ outputs
# NORILSK_DATA / NORILSK_OUT override the folders (used by the synthetic test).
DATA = Path(os.environ.get("NORILSK_DATA", Path(__file__).parent / "data"))
OUTPUTS = Path(os.environ.get("NORILSK_OUT", Path(__file__).parent / "outputs"))
DIRS = {k: DATA / k for k in
        ["s1", "hyp3", "nisar", "climate", "dem", "landcover", "osm", "cci",
         "zones", "ts"]}

# Radar wavelengths (m). NISAR is read from the file when available.
S1_WAVELENGTH = 0.05546576
# Zone -> ERA5-Land point used for its TDD.
ZONE_CLIMATE = {"reference": "plateau_reference", "industrial": "norilsk_centre",
                "accident": "norilsk_centre", "urban": "norilsk_centre",
                "control_NW": "control_NW", "control_SE": "control_SE"}
for d in DIRS.values():
    d.mkdir(parents=True, exist_ok=True)


def bbox_wkt(b=BBOX) -> str:
    x0, y0, x1, y1 = b
    return f"POLYGON(({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))"
